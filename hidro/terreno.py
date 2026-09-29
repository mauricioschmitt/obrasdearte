"""
Análise de terreno para delimitação de microbacias.

Implementação em NumPy puro, sem dependência de GRASS/TauDEM, para que o
aluno consiga ler e modificar cada etapa. As três funções pesadas rodam
UMA VEZ no pré-processamento (por estado / por folha), nunca durante a
requisição do usuário.

Etapas:
  1. preencher_depressoes  -> Priority-Flood + epsilon (Barnes et al., 2014)
  2. direcao_d8            -> vetor de "receptores" (para onde cada célula drena)
  3. acumulacao_fluxo      -> nº de células a montante
  4. construir_doadores    -> grafo invertido em formato CSR, usado no traçado

Convenção: o terreno é representado como um grafo de células. Cada célula
tem no máximo um receptor (D8) e N doadores. Delimitar uma bacia é
percorrer os doadores a partir do exutório.
"""

from __future__ import annotations

import heapq

import numpy as np

# Deslocamentos dos 8 vizinhos, em ordem fixa (usada em todo o módulo)
DESLOC = np.array(
    [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)],
    dtype=np.int64,
)


def preencher_depressoes(
    mde: np.ndarray,
    mascara_valida: np.ndarray,
    epsilon: float = 1e-4,
) -> np.ndarray:
    """Remove depressões espúrias e áreas planas do MDE.

    Usa Priority-Flood com incremento epsilon: além de encher as depressões,
    impõe um gradiente mínimo nas áreas planas, garantindo que toda célula
    interna tenha ao menos um vizinho mais baixo. Sem isso, planícies e
    reservatórios viram buracos negros no traçado de fluxo.

    Parâmetros
    ----------
    mde : matriz float32 de elevações (m)
    mascara_valida : booleano, True onde há dado
    epsilon : incremento mínimo entre células vizinhas (m)

    Retorna
    -------
    Matriz de elevações condicionada hidrologicamente.
    """
    ny, nx = mde.shape
    preenchido = mde.astype(np.float64, copy=True)
    preenchido[~mascara_valida] = np.inf

    fechada = ~mascara_valida.copy()
    fila: list[tuple[float, int, int]] = []

    # Semeia a fila com a borda do domínio válido: são os pontos por onde a
    # água sai do recorte.
    for i in range(ny):
        for j in range(nx):
            if not mascara_valida[i, j]:
                continue
            na_borda = i == 0 or j == 0 or i == ny - 1 or j == nx - 1
            if not na_borda:
                # também é borda se faz divisa com nodata
                for di, dj in DESLOC:
                    vi, vj = i + di, j + dj
                    if not mascara_valida[vi, vj]:
                        na_borda = True
                        break
            if na_borda:
                heapq.heappush(fila, (float(preenchido[i, j]), i, j))
                fechada[i, j] = True

    while fila:
        z, i, j = heapq.heappop(fila)
        for di, dj in DESLOC:
            vi, vj = i + di, j + dj
            if vi < 0 or vj < 0 or vi >= ny or vj >= nx:
                continue
            if fechada[vi, vj]:
                continue
            zv = preenchido[vi, vj]
            if zv <= z:
                zv = z + epsilon
                preenchido[vi, vj] = zv
            fechada[vi, vj] = True
            heapq.heappush(fila, (float(zv), vi, vj))

    preenchido[~mascara_valida] = np.nan
    return preenchido.astype(np.float32)


def direcao_d8(
    mde: np.ndarray,
    mascara_valida: np.ndarray,
    tam_celula: float,
) -> np.ndarray:
    """Calcula o receptor D8 de cada célula.

    Retorna um vetor achatado onde receptor[k] é o índice da célula para
    onde k drena, ou -1 se k é exutório (sai do recorte) ou nodata.

    A escolha é pelo maior gradiente, corrigindo a distância nas diagonais.
    """
    ny, nx = mde.shape
    z = np.where(mascara_valida, mde, np.inf).astype(np.float64)

    melhor_grad = np.full((ny, nx), -np.inf)
    receptor = np.full((ny, nx), -1, dtype=np.int64)

    for di, dj in DESLOC:
        dist = tam_celula * np.hypot(di, dj)

        # z deslocado, com inf fora do domínio (impede drenar para fora)
        zv = np.full((ny, nx), np.inf)
        oi_dst = slice(max(0, di), ny + min(0, di))
        oj_dst = slice(max(0, dj), nx + min(0, dj))
        oi_src = slice(max(0, -di), ny + min(0, -di))
        oj_src = slice(max(0, -dj), nx + min(0, -dj))
        zv[oi_src, oj_src] = z[oi_dst, oj_dst]

        grad = (z - zv) / dist

        idx_vizinho = np.full((ny, nx), -1, dtype=np.int64)
        idx_plano = np.arange(ny * nx, dtype=np.int64).reshape(ny, nx)
        idx_vizinho[oi_src, oj_src] = idx_plano[oi_dst, oj_dst]

        melhora = (grad > melhor_grad) & (grad > 0) & np.isfinite(grad)
        melhor_grad = np.where(melhora, grad, melhor_grad)
        receptor = np.where(melhora, idx_vizinho, receptor)

    receptor[~mascara_valida] = -1
    return receptor.ravel()


def acumulacao_fluxo(
    receptor: np.ndarray,
    forma: tuple[int, int],
    mascara_valida: np.ndarray,
) -> np.ndarray:
    """Conta quantas células drenam para cada célula (inclusive ela mesma).

    Processa as células em ordem topológica usando o algoritmo de Kahn sobre
    o grafo de fluxo, o que é O(n) e dispensa recursão.
    """
    n = receptor.size
    valido = mascara_valida.ravel()

    grau_entrada = np.zeros(n, dtype=np.int32)
    tem_receptor = receptor >= 0
    np.add.at(grau_entrada, receptor[tem_receptor], 1)

    acumulacao = np.where(valido, 1.0, 0.0)

    pilha = list(np.flatnonzero((grau_entrada == 0) & valido))
    while pilha:
        k = pilha.pop()
        r = receptor[k]
        if r < 0:
            continue
        acumulacao[r] += acumulacao[k]
        grau_entrada[r] -= 1
        if grau_entrada[r] == 0:
            pilha.append(r)

    return acumulacao.reshape(forma)


def construir_doadores(receptor: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Inverte o grafo de fluxo para o formato CSR.

    Esta é a estrutura que torna a delimitação instantânea em tempo de
    requisição: a partir do exutório, basta percorrer os doadores.

    Retorna (inicio, doadores) tal que os doadores da célula k são
    doadores[inicio[k]:inicio[k+1]].
    """
    n = receptor.size
    tem_receptor = receptor >= 0
    destinos = receptor[tem_receptor]
    origens = np.flatnonzero(tem_receptor)

    contagem = np.bincount(destinos, minlength=n)
    inicio = np.zeros(n + 1, dtype=np.int64)
    np.cumsum(contagem, out=inicio[1:])

    ordem = np.argsort(destinos, kind="stable")
    doadores = origens[ordem].astype(np.int64)

    return inicio, doadores
