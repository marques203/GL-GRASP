"""Modulo de embeddings do GL-GRASP (ver FIGURA 3.1 do TG):

    Tecnica de GRL -> Reducao de dimensionalidade -> Distancias estruturais

Le uma instancia do C-IGDP, gera embeddings HOPE (CogDL) para todos os
vertices de IG, projeta para 2D via PCA e calcula a distancia Euclidiana
entre cada par de vertices adjacentes (arco), exportando o resultado em um
arquivo texto para ser consumido depois pelo modulo de otimizacao (C++).

Uso:
    python build_distances.py <instancia.txt> [saida.txt] [--dim N] [--beta B]

Se <saida.txt> nao for informado, o resultado vai para
embeddings/distances/<nome_da_instancia>.dist.txt (a pasta e criada se
nao existir).

Formato do arquivo de saida:
    n_arcos
    layer_u id_u layer_v id_v distancia
    ...
"""
import argparse
from pathlib import Path

import numpy as np
from sklearn.decomposition import PCA

from hdag_io import read_instance
from hope_embeddings import compute_hope_embeddings

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


def write_distances(path: str, distances: list) -> None:
    with open(path, "w") as f:
        f.write(f"{len(distances)}\n")
        for level_u, id_u, level_v, id_v, d in distances:
            f.write(f"{level_u} {id_u} {level_v} {id_v} {d:.6f}\n")


def main():
    parser = argparse.ArgumentParser(description="Gera a matriz de distancias estruturais (GL-GRASP - modulo Python)")
    parser.add_argument("instance", help="arquivo de instancia do C-IGDP")
    parser.add_argument(
        "output",
        nargs="?",
        default=None,
        help=f"arquivo de saida com as distancias por arco (default: {DEFAULT_OUTPUT_DIR}/<instancia>.dist.txt)",
    )
    parser.add_argument("--dim", type=int, default=None, help="dimensao do embedding HOPE (default: 2n-1)")
    parser.add_argument("--beta", type=float, default=0.01, help="parametro beta do indice de Katz (default: 0.01)")
    args = parser.parse_args()

    if args.output is None:
        DEFAULT_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        output_path = DEFAULT_OUTPUT_DIR / f"{Path(args.instance).stem}.dist.txt"
    else:
        output_path = Path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)

    instance = read_instance(args.instance)
    print(f"instancia: {instance.num_levels} niveis, {instance.num_nodes} vertices, {len(instance.arcs)} arcos")

    embeddings = compute_hope_embeddings(instance, dimension=args.dim, beta=args.beta)
    print(f"embeddings HOPE gerados, dimensao = {next(iter(embeddings.values())).shape[0]}")

    coords_2d = project_to_2d(embeddings)

    distances = compute_arc_distances(instance, coords_2d)
    write_distances(output_path, distances)
    print(f"{len(distances)} distancias de arco escritas em {output_path}")


if __name__ == "__main__":
    main()
