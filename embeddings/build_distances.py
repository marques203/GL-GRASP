"""Modulo de embeddings do GL-GRASP (ver FIGURA 3.1 do TG):

    Tecnica de GRL -> Reducao de dimensionalidade -> Distancias estruturais

Le uma instancia do C-IGDP, gera embeddings para todos os vertices de IG
usando uma (ou todas) das 8 tecnicas da tabela do TG, projeta para 2D via
PCA e calcula a distancia Euclidiana entre cada par de vertices adjacentes
(arco), exportando o resultado em um arquivo texto para ser consumido
depois pelo modulo de otimizacao (C++).

Uso:
    python build_distances.py <instancia.txt> [saida.txt] [--technique T] [--dim N] [--param k=v ...]

--technique aceita: spectral, hope, node2vec, sdne, line, grarep, netsmf,
prone, ou "all" (roda as 8 e gera um arquivo de saida por tecnica).

Se <saida.txt> nao for informado (ou --technique all), o resultado vai
para embeddings/distances/<instancia>.<tecnica>.dist.txt.

Formato do arquivo de saida:
    n_arcos
    layer_u id_u layer_v id_v distancia
    ...
"""
import argparse
import time
import traceback
from pathlib import Path

import numpy as np
from sklearn.decomposition import PCA

from hdag_io import read_instance
from embedding_techniques import TECHNIQUES, compute_embeddings

DEFAULT_OUTPUT_DIR = Path(__file__).parent / "distances"


def project_to_2d(embeddings: dict) -> dict:
    keys = list(embeddings.keys())
    matrix = np.stack([embeddings[k] for k in keys])

    coords = PCA(n_components=2).fit_transform(matrix)

    return {k: coords[i] for i, k in enumerate(keys)}


def compute_arc_distances(instance, coords_2d: dict) -> list:
    distances = []
    for (level_u, id_u), (level_v, id_v) in instance.arcs:
        p = coords_2d[(level_u, id_u)]
        q = coords_2d[(level_v, id_v)]
        d = float(np.linalg.norm(p - q))
        distances.append((level_u, id_u, level_v, id_v, d))
    return distances


def write_distances(path: Path, distances: list) -> None:
    with open(path, "w") as f:
        f.write(f"{len(distances)}\n")
        for level_u, id_u, level_v, id_v, d in distances:
            f.write(f"{level_u} {id_u} {level_v} {id_v} {d:.6f}\n")


def default_output_path(instance_path: str, technique: str) -> Path:
    DEFAULT_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    return DEFAULT_OUTPUT_DIR / f"{Path(instance_path).stem}.{technique}.dist.txt"


def run_one(instance, instance_path: str, technique: str, output_path: Path, dimension, params) -> None:
    t0 = time.time()
    embeddings = compute_embeddings(instance, technique, dimension=dimension, **params)
    dim = next(iter(embeddings.values())).shape[0]

    coords_2d = project_to_2d(embeddings)
    distances = compute_arc_distances(instance, coords_2d)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    write_distances(output_path, distances)

    print(
        f"[{technique}] dim={dim} arcos={len(distances)} tempo={time.time() - t0:.1f}s -> {output_path}"
    )


def parse_params(pairs: list) -> dict:
    params = {}
    for item in pairs or []:
        key, _, value = item.partition("=")
        try:
            value = float(value) if "." in value else int(value)
        except ValueError:
            pass
        params[key] = value
    return params


def main():
    parser = argparse.ArgumentParser(description="Gera a matriz de distancias estruturais (GL-GRASP - modulo Python)")
    parser.add_argument("instance", help="arquivo de instancia do C-IGDP")
    parser.add_argument(
        "output",
        nargs="?",
        default=None,
        help=f"arquivo de saida (default: {DEFAULT_OUTPUT_DIR}/<instancia>.<tecnica>.dist.txt; ignorado se --technique all)",
    )
    parser.add_argument(
        "--technique", "-t", default="hope", choices=TECHNIQUES + ["all"], help="tecnica de embedding (default: hope)"
    )
    parser.add_argument("--dim", type=int, default=None, help="dimensao do embedding (default: especifico por tecnica)")
    parser.add_argument(
        "--param", action="append", metavar="k=v", help="sobrescreve um parametro especifico da tecnica (repetivel)"
    )
    args = parser.parse_args()
    params = parse_params(args.param)

    instance = read_instance(args.instance)
    print(f"instancia: {instance.num_levels} niveis, {instance.num_nodes} vertices, {len(instance.arcs)} arcos")

    techniques = TECHNIQUES if args.technique == "all" else [args.technique]

    failures = []
    for technique in techniques:
        output_path = (
            Path(args.output)
            if (args.output and args.technique != "all")
            else default_output_path(args.instance, technique)
        )
        try:
            run_one(instance, args.instance, technique, output_path, args.dim, params)
        except Exception as exc:  # noqa: BLE001 - queremos continuar as outras tecnicas mesmo se uma falhar
            print(f"[{technique}] FALHOU: {exc}")
            traceback.print_exc()
            failures.append(technique)

    if len(techniques) > 1:
        ok = len(techniques) - len(failures)
        print(f"\nresumo: {ok}/{len(techniques)} tecnicas concluidas" + (f", falharam: {failures}" if failures else ""))


if __name__ == "__main__":
    main()
