"""
Delimitação de microbacia a partir de um ponto e cálculo dos parâmetros
morfométricos que alimentam o módulo de vazão.

Tudo aqui roda em tempo de requisição e precisa ser rápido. O traçado a
montante visita apenas as células da bacia, então o custo é proporcional ao
tamanho da bacia, não ao tamanho do MDE do estado inteiro.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict

import numpy as np
from rasterio import features
from rasterio.transform import rowcol, xy
from shapely.geometry import shape, mapping
from shapely.ops import unary_union


@dataclass
class Morfometria:
    """Parâmetros da bacia. São exatamente as entradas do cálculo de Tc."""

    area_km2: float
    perimetro_km: float
    comprimento_talvegue_km: float
    cota_maxima_m: float
    cota_exutorio_m: float
    desnivel_m: float
    declividade_talvegue_m_m: float
    declividade_media_bacia_pct: float
    coef_compacidade: float
    fator_forma: float
    n_celulas: int

    def dict(self) -> dict:
        return {k: (round(v, 4) if isinstance(v, float) else v) for k, v in asdict(self).items()}


def ancorar_exutorio(
    linha: int,
    coluna: int,
    acumulacao: np.ndarray,
    raio: int = 5,
) -> tuple[int, int]:
    """Move o ponto clicado para a célula de maior acumulação na vizinhança.

    O usuário clica em cima do rio no mapa base, mas o rio do MDE está
    deslocado alguns pixels. Sem esta ancoragem, metade dos cliques delimita
    uma bacia de 3 pixels na encosta ao lado. É o erro mais comum em
    ferramentas desse tipo.
    """
    ny, nx = acumulacao.shape
    i0, i1 = max(0, linha - raio), min(ny, linha + raio + 1)
    j0, j1 = max(0, coluna - raio), min(nx, coluna + raio + 1)

    janela = acumulacao[i0:i1, j0:j1]
    if janela.size == 0:
        return linha, coluna

    di, dj = np.unravel_index(np.nanargmax(janela), janela.shape)
    return i0 + int(di), j0 + int(dj)


def celulas_a_montante(
    exutorio: int,
    inicio: np.ndarray,
    doadores: np.ndarray,
) -> np.ndarray:
    """Percorre o grafo invertido e devolve todas as células da bacia."""
    bacia: list[int] = []
    pilha = [exutorio]
    while pilha:
        k = pilha.pop()
        bacia.append(k)
        pilha.extend(doadores[inicio[k] : inicio[k + 1]].tolist())
    return np.asarray(bacia, dtype=np.int64)


def caminho_de_fluxo_mais_longo(
    celulas: np.ndarray,
    receptor: np.ndarray,
    exutorio: int,
    forma: tuple[int, int],
    tam_celula: float,
) -> tuple[float, int]:
    """Comprimento do talvegue principal, pelo maior caminho de fluxo.

    Calcula a distância de cada célula até o exutório seguindo o próprio
    grafo de fluxo. O maior valor é o comprimento do talvegue, definição
    mais defensável que "maior eixo da bacia".

    Retorna (comprimento em metros, índice da célula de cabeceira).
    """
    ny, nx = forma
    distancia = {int(exutorio): 0.0}

    # Ordena as células por profundidade no grafo garantindo que o receptor
    # já tenha distância conhecida quando a célula for processada.
    ordem = _ordem_de_jusante_para_montante(celulas, receptor, exutorio)

    for k in ordem:
        r = receptor[k]
        if r < 0 or int(r) not in distancia:
            continue
        i, j = divmod(int(k), nx)
        ri, rj = divmod(int(r), nx)
        passo = tam_celula * np.hypot(i - ri, j - rj)
        distancia[int(k)] = distancia[int(r)] + passo

    cabeceira = max(distancia, key=distancia.get)
    return distancia[cabeceira], cabeceira


def _ordem_de_jusante_para_montante(
    celulas: np.ndarray, receptor: np.ndarray, exutorio: int
) -> list[int]:
    """Ordena as células da bacia de jusante para montante (BFS)."""
    conjunto = set(celulas.tolist())
    doadores_locais: dict[int, list[int]] = {}
    for k in celulas.tolist():
        r = int(receptor[k])
        if r in conjunto:
            doadores_locais.setdefault(r, []).append(int(k))

    ordem: list[int] = []
    fila = [int(exutorio)]
    while fila:
        k = fila.pop(0)
        filhos = doadores_locais.get(k, [])
        ordem.extend(filhos)
        fila.extend(filhos)
    return ordem


def poligono_da_bacia(
    celulas: np.ndarray,
    forma: tuple[int, int],
    transform,
):
    """Converte as células da bacia em um polígono, suavizando os degraus."""
    mascara = np.zeros(forma, dtype=np.uint8)
    mascara.ravel()[celulas] = 1

    geoms = [
        shape(geom)
        for geom, valor in features.shapes(mascara, mask=mascara.astype(bool), transform=transform)
        if valor == 1
    ]
    if not geoms:
        raise ValueError("Nenhuma célula na bacia delimitada.")

    poligono = unary_union(geoms)
    if poligono.geom_type == "MultiPolygon":
        poligono = max(poligono.geoms, key=lambda g: g.area)

    # Remove os degraus de pixel sem deslocar o contorno de forma relevante
    return poligono.simplify(transform.a * 0.5, preserve_topology=True)


def calcular_morfometria(
    celulas: np.ndarray,
    poligono,
    mde: np.ndarray,
    receptor: np.ndarray,
    exutorio: int,
    forma: tuple[int, int],
    tam_celula: float,
) -> Morfometria:
    """Reúne os parâmetros da bacia a partir das células e do MDE."""
    ny, nx = forma
    area_m2 = len(celulas) * tam_celula**2
    area_km2 = area_m2 / 1e6
    perimetro_km = poligono.length / 1000.0

    z = mde.ravel()
    cotas = z[celulas]
    cotas = cotas[np.isfinite(cotas)]

    comprimento_m, cabeceira = caminho_de_fluxo_mais_longo(
        celulas, receptor, exutorio, forma, tam_celula
    )
    comprimento_km = comprimento_m / 1000.0

    cota_exutorio = float(z[exutorio])
    cota_cabeceira = float(z[cabeceira])
    desnivel = max(cota_cabeceira - cota_exutorio, 0.1)
    decl_talvegue = desnivel / max(comprimento_m, 1.0)

    decl_media_pct = _declividade_media(celulas, mde, forma, tam_celula)

    # Kc = 0,28 P / sqrt(A), com P em km e A em km2
    kc = 0.28 * perimetro_km / np.sqrt(area_km2) if area_km2 > 0 else float("nan")
    # Kf = A / L2
    kf = area_km2 / comprimento_km**2 if comprimento_km > 0 else float("nan")

    return Morfometria(
        area_km2=area_km2,
        perimetro_km=perimetro_km,
        comprimento_talvegue_km=comprimento_km,
        cota_maxima_m=float(cotas.max()) if cotas.size else float("nan"),
        cota_exutorio_m=cota_exutorio,
        desnivel_m=desnivel,
        declividade_talvegue_m_m=decl_talvegue,
        declividade_media_bacia_pct=decl_media_pct,
        coef_compacidade=float(kc),
        fator_forma=float(kf),
        n_celulas=int(len(celulas)),
    )


def _declividade_media(
    celulas: np.ndarray, mde: np.ndarray, forma: tuple[int, int], tam_celula: float
) -> float:
    """Declividade média das encostas pelo gradiente de Horn simplificado."""
    gy, gx = np.gradient(np.nan_to_num(mde, nan=0.0), tam_celula)
    declividade = np.hypot(gx, gy).ravel()
    valores = declividade[celulas]
    valores = valores[np.isfinite(valores)]
    return float(np.mean(valores) * 100.0) if valores.size else float("nan")


def tempo_concentracao_previa(m: Morfometria) -> dict:
    """Prévia do próximo módulo: Tc por três fórmulas, em minutos.

    Não é o módulo de vazão, é só para mostrar que os parâmetros da
    delimitação já estão prontos para alimentá-lo, e para deixar visível a
    dispersão entre fórmulas.
    """
    L_km = max(m.comprimento_talvegue_km, 1e-3)
    L_m = L_km * 1000.0
    S = max(m.declividade_talvegue_m_m, 1e-5)

    kirpich = 57.0 * (L_km**3 / m.desnivel_m) ** 0.385
    dooge = 21.88 * (m.area_km2**0.41) * (S**-0.17)
    corpo_eng = 191.5 * (L_km**0.76) * (m.desnivel_m / L_m) ** -0.19 / 60.0

    return {
        "kirpich_min": round(kirpich, 1),
        "dooge_min": round(dooge, 1),
        "corps_of_engineers_min": round(corpo_eng, 1),
        "observacao": "Três fórmulas, três respostas. A escolha precisa ser justificada na dissertação.",
    }


def ponto_para_celula(lon_x: float, lat_y: float, transform) -> tuple[int, int]:
    return rowcol(transform, lon_x, lat_y)


def celula_para_ponto(linha: int, coluna: int, transform) -> tuple[float, float]:
    return xy(transform, linha, coluna)


def para_geojson(poligono, propriedades: dict) -> dict:
    return {
        "type": "Feature",
        "geometry": mapping(poligono),
        "properties": propriedades,
    }
