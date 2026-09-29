#!/usr/bin/env python3
"""
Converte dados reais (MDE + uso do solo) no pacote que o site estático lê.

Roda UMA VEZ por região, na máquina de vocês. O resultado são dois arquivos
pequenos que você copia para a pasta do site. Depois disso não existe mais
servidor, nem processamento, nem custo.

    python preparar_regiao.py \
        --mde recortes/rh08_mde.tif \
        --uso recortes/rh08_mapbiomas.tif \
        --id rh08 --nome "RH8 · Vale do Itajaí" \
        --idf 1132,0.159,12.0,0.771 \
        --saida docs/dados

IMPORTANTE — recorte por região hidrográfica, não por município.
Uma bacia nunca cruza um divisor de águas, então recortando por região
hidrográfica toda microbacia cai inteira dentro de um arquivo. Recortando
por município, metade das bacias fica cortada ao meio e o cálculo sai errado
sem avisar. SC tem 10 regiões hidrográficas oficiais (RH1 a RH10).
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
import time
from pathlib import Path

import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.warp import reproject

sys.path.insert(0, str(Path(__file__).resolve().parent))
from hidro.terreno import DESLOC, direcao_d8, preencher_depressoes

# ---------------------------------------------------------------------------
# De MapBiomas para as 8 classes usadas no CN.
#
# Confira estes códigos contra a legenda da coleção que você baixar: eles
# mudam entre coleções, e um código não mapeado vira pastagem silenciosamente
# (o script avisa quando isso acontece).
# ---------------------------------------------------------------------------
MAPA_MAPBIOMAS = {
    1: [1, 3, 4, 5, 49],                                   # Floresta
    2: [9],                                                # Silvicultura
    3: [10, 11, 12, 23, 29, 32, 50],                       # Campo natural
    4: [15, 21],                                           # Pastagem
    5: [14, 18, 19, 20, 35, 36, 39, 40, 41, 46, 47, 48, 62],  # Agricultura
    6: [22, 25, 30],                                       # Solo exposto
    7: [24],                                               # Área urbanizada
    8: [26, 31, 33],                                       # Corpo d'água
}
CLASSE_PADRAO = 4
NOMES = {1: "Floresta", 2: "Silvicultura", 3: "Campo natural", 4: "Pastagem",
         5: "Agricultura", 6: "Solo exposto", 7: "Área urbanizada", 8: "Corpo d'água"}

LIMITE_AVISO_CELULAS = 6_000_000


def remapear_uso(bruto: np.ndarray) -> tuple[np.ndarray, dict]:
    """Traduz os códigos do MapBiomas para as 8 classes hidrológicas."""
    saida = np.full(bruto.shape, CLASSE_PADRAO, dtype=np.uint8)
    conhecidos = set()
    for destino, origens in MAPA_MAPBIOMAS.items():
        for cod in origens:
            saida[bruto == cod] = destino
            conhecidos.add(cod)

    presentes = set(np.unique(bruto).tolist())
    nao_mapeados = {
        int(c): int((bruto == c).sum())
        for c in presentes - conhecidos - {0, 27, 255}
        if (bruto == c).sum() > 0
    }
    return saida, nao_mapeados


def alinhar(caminho_uso: str, perfil_mde: dict) -> np.ndarray:
    """Reamostra o uso do solo exatamente sobre a grade do MDE."""
    destino = np.zeros((perfil_mde["height"], perfil_mde["width"]), dtype=np.int32)
    with rasterio.open(caminho_uso) as src:
        reproject(
            source=rasterio.band(src, 1),
            destination=destino,
            src_transform=src.transform, src_crs=src.crs,
            dst_transform=perfil_mde["transform"], dst_crs=perfil_mde["crs"],
            resampling=Resampling.nearest,   # classe temática: nunca interpolar
        )
    return destino


def preencher(mde: np.ndarray, mascara: np.ndarray) -> np.ndarray:
    """Usa richdem se estiver instalado, senão a implementação própria."""
    try:
        import richdem as rd
        arr = rd.rdarray(np.where(mascara, mde, -9999).astype(np.float64), no_data=-9999)
        return np.asarray(rd.FillDepressions(arr, epsilon=True)).astype(np.float32)
    except ImportError:
        print("      (richdem não instalado, usando a versão em NumPy — mais lenta)")
        return preencher_depressoes(mde, mascara)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mde", required=True, help="GeoTIFF de elevação, CRS projetado em metros")
    ap.add_argument("--uso", required=True, help="GeoTIFF do MapBiomas recortado na mesma área")
    ap.add_argument("--id", required=True, help="identificador curto, ex.: rh08")
    ap.add_argument("--nome", required=True, help="nome exibido no seletor do site")
    ap.add_argument("--saida", default="docs/dados")
    ap.add_argument("--idf", default="1000,0.15,10,0.75",
                    help="coeficientes K,a,b,c da curva i = K·TR^a/(t+b)^c")
    ap.add_argument("--reamostrar", type=float, default=0,
                    help="tamanho de célula alvo em metros (0 = manter o original)")
    args = ap.parse_args()

    t0 = time.time()
    saida = Path(args.saida)
    saida.mkdir(parents=True, exist_ok=True)

    # --- 1. MDE -------------------------------------------------------------
    print("[1/6] lendo o MDE")
    with rasterio.open(args.mde) as src:
        if src.crs is None or src.crs.is_geographic:
            raise SystemExit(
                f"O MDE precisa estar em CRS projetado em metros (CRS atual: {src.crs}).\n"
                "Em SC use SIRGAS 2000 / UTM 22S (EPSG:31982), ou 23S (EPSG:31983) no litoral norte.\n"
                "  gdalwarp -t_srs EPSG:31982 -r bilinear entrada.tif saida.tif"
            )

        fator = 1.0
        if args.reamostrar > 0:
            fator = abs(src.transform.a) / args.reamostrar
        altura = max(1, int(src.height * fator))
        largura = max(1, int(src.width * fator))

        mde = src.read(1, out_shape=(altura, largura),
                       resampling=Resampling.bilinear).astype(np.float32)
        transform = src.transform * src.transform.scale(src.width / largura,
                                                        src.height / altura)
        crs = src.crs
        nodata = src.nodata
        perfil = {"height": altura, "width": largura, "transform": transform, "crs": crs}

    celula = abs(transform.a)
    if abs(abs(transform.e) - celula) > 1e-6:
        raise SystemExit("As células não são quadradas. Reamostre o MDE antes.")

    lin, col = mde.shape
    ncel = lin * col
    print(f"      {lin} × {col} = {ncel:,} células de {celula:.1f} m  ·  {crs}")
    print(f"      área do recorte: {ncel * celula**2 / 1e6:,.0f} km²")

    if ncel > LIMITE_AVISO_CELULAS:
        print(f"      ATENÇÃO: acima de {LIMITE_AVISO_CELULAS:,} células o pré-processamento")
        print("      fica lento e o download do usuário fica pesado. Considere")
        print("      --reamostrar 30, ou recortar em regiões hidrográficas menores.")

    mascara = np.isfinite(mde)
    if nodata is not None:
        mascara &= mde != nodata
    if not mascara.any():
        raise SystemExit("O MDE está todo vazio. Confira o recorte e o valor de nodata.")

    # --- 2. Uso do solo -----------------------------------------------------
    print("[2/6] alinhando o uso do solo sobre a grade do MDE")
    uso_bruto = alinhar(args.uso, perfil)
    uso, nao_mapeados = remapear_uso(uso_bruto)
    uso[~mascara] = CLASSE_PADRAO

    if nao_mapeados:
        print("      códigos do MapBiomas sem tradução (viraram pastagem):")
        for cod, n in sorted(nao_mapeados.items(), key=lambda x: -x[1])[:8]:
            print(f"        código {cod:4d}  {n / ncel * 100:5.2f} % da área")
        print("      Ajuste MAPA_MAPBIOMAS no topo deste arquivo antes de usar em projeto.")

    # --- 3 e 4. Hidrologia --------------------------------------------------
    print("[3/6] preenchendo depressões")
    mde_cond = preencher(mde, mascara)

    print("[4/6] calculando a direção de fluxo D8")
    receptor = direcao_d8(mde_cond, mascara, celula)

    # --- 5. Codificação -----------------------------------------------------
    print("[5/6] codificando o pacote")
    d8 = np.zeros(ncel, dtype=np.uint8)
    idx = np.arange(ncel)
    tem = receptor >= 0
    di = receptor[tem] // col - idx[tem] // col
    dj = receptor[tem] % col - idx[tem] % col
    for codigo, (ddi, ddj) in enumerate(DESLOC, start=1):
        d8[idx[tem][(di == ddi) & (dj == ddj)]] = codigo

    mde_i16 = np.clip(np.nan_to_num(mde_cond, nan=0), -32000, 32000).round().astype("<i2")

    bruto = mde_i16.tobytes() + d8.tobytes() + uso.ravel().tobytes()
    arquivo_bin = saida / f"{args.id}.bin.gz"
    arquivo_bin.write_bytes(gzip.compress(bruto, 6))

    K, a, b, c = (float(x) for x in args.idf.split(","))
    x0, y0 = transform.c, transform.f

    # rasterio emite "+south=True" e "+no_defs=True"; a sintaxe correta é a
    # flag sozinha. Alguns leitores de proj4 recusam a forma com "=True".
    proj4 = crs.to_proj4().replace("+south=True", "+south").replace("+no_defs=True", "+no_defs")

    meta = {
        "id": args.id, "nome": args.nome,
        "lin": lin, "col": col, "celula": celula,
        "origem": [x0, y0],
        "proj4": proj4, "epsg": crs.to_epsg(),
        "idf": {"K": K, "a": a, "b": b, "c": c},
        "bin": arquivo_bin.name,
        "fonte": {"mde": Path(args.mde).name, "uso": Path(args.uso).name,
                  "gerado_em": time.strftime("%Y-%m-%d")},
    }
    (saida / f"{args.id}.json").write_text(json.dumps(meta, ensure_ascii=False, indent=1),
                                           encoding="utf-8")

    # --- 6. Manifesto -------------------------------------------------------
    print("[6/6] atualizando o manifesto")
    manifesto = saida / "manifesto.json"
    regioes = []
    if manifesto.exists():
        regioes = [r for r in json.loads(manifesto.read_text(encoding="utf-8"))["regioes"]
                   if r["id"] != args.id]
    regioes.append({"id": args.id, "nome": args.nome, "meta": f"{args.id}.json"})
    regioes.sort(key=lambda r: r["nome"])
    manifesto.write_text(json.dumps({"regioes": regioes}, ensure_ascii=False, indent=1),
                         encoding="utf-8")

    # --- Relatório ----------------------------------------------------------
    mb = arquivo_bin.stat().st_size / 1e6
    print(f"\nPronto em {time.time() - t0:.0f} s")
    print(f"  {arquivo_bin}   {mb:.2f} MB  (é o que o usuário baixa)")
    print(f"  {saida / (args.id + '.json')}")
    print(f"  {manifesto}")
    print("\nComposição do uso do solo:")
    for cod, nome in NOMES.items():
        pct = (uso == cod).sum() / ncel * 100
        if pct >= 0.05:
            print(f"  {nome:18s} {pct:5.1f} %")
    if mb > 8:
        print(f"\nATENÇÃO: {mb:.1f} MB é um download pesado para conexão de prefeitura.")
        print("Recorte em regiões menores ou use --reamostrar 30.")


if __name__ == "__main__":
    main()
