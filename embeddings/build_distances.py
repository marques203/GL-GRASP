"""Modulo de embeddings do GL-GRASP (ver FIGURA 3.1 do TG):

    Tecnica de GRL -> Reducao de dimensionalidade -> Distancias estruturais

Le uma instancia do C-IGDP, gera embeddings para todos os vertices de IG
usando uma (ou todas) das 8 tecnicas da tabela do TG, projeta para n-D via
PCA e calcula a distancia entre cada par de vertices adjacentes (arco),
exportando o resultado em um arquivo texto para ser consumido depois pelo
modulo de otimizacao (C++).

Uso:
    python build_distances.py <instancia.txt> [saida.txt] [--technique T]
        [--dim N] [--pca-dim N] [--metric M] [--param k=v ...]

--technique aceita: spectral, hope, node2vec, sdne, line, grarep, netsmf,
prone, ou "all" (roda as 8 e gera um arquivo de saida por tecnica).

--dim e a dimensao do EMBEDDING; --pca-dim e a dimensao da projecao PCA
aplicada sobre ele. Sao independentes: o artigo usa embedding de k dimensoes
(2n-1 para HOPE) projetado sobre um plano, ou seja --pca-dim 2 (default).

Se <saida.txt> nao for informado (ou --technique all), o resultado vai
para embeddings/distances/<instancia>.<tecnica>.dist.txt. Quando --metric
ou --pca-dim diferem do default, o nome recebe sufixos correspondentes,
para que execucoes distintas nao sobrescrevam umas as outras.

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

DEFAULT_PCA_DIM = 2
DEFAULT_METRIC = "euclidean"

METRICS = ["euclidean", "manhattan", "chebyshev", "cosine"]


def distance_metric(p, q, metric: str = DEFAULT_METRIC) -> float:
    if metric == "euclidean":
        return float(np.linalg.norm(p - q))
    elif metric == "manhattan":
        return float(np.linalg.norm(p - q, ord=1))
    elif metric == "chebyshev":
        return float(np.linalg.norm(p - q, ord=np.inf))
    elif metric == "cosine":
        norm_p = np.linalg.norm(p)
        norm_q = np.linalg.norm(q)
        if norm_p == 0.0 or norm_q == 0.0:
            return 1.0
        return float(1.0 - np.dot(p, q) / (norm_p * norm_q))
    else:
        raise ValueError(f"Metrica desconhecida: {metric}")


def project_to_nd(embeddings: dict, n: int = DEFAULT_PCA_DIM) -> dict:
    keys = list(embeddings.keys())
    matrix = np.stack([embeddings[k] for k in keys])

    # O PCA nao consegue extrair mais componentes do que min(amostras, features);
    # pedir mais do que isso levanta excecao no scikit-learn.
    max_components = min(matrix.shape)
    if n > max_components:
        raise ValueError(
            f"--pca-dim {n} excede o maximo possivel para esta instancia/tecnica "
            f"({max_components} = min({matrix.shape[0]} vertices, {matrix.shape[1]} dims do embedding))"
        )

    coords = PCA(n_components=n).fit_transform(matrix)

    return {k: coords[i] for i, k in enumerate(keys)}


def compute_arc_distances(instance, coords_nd: dict, metric: str = DEFAULT_METRIC) -> list:
    distances = []
    for (level_u, id_u), (level_v, id_v) in instance.arcs:
        p = coords_nd[(level_u, id_u)]
        q = coords_nd[(level_v, id_v)]
        d = distance_metric(p, q, metric)
        distances.append((level_u, id_u, level_v, id_v, d))
    return distances


def write_distances(path: Path, distances: list) -> None:
    with open(path, "w") as f:
        f.write(f"{len(distances)}\n")
        for level_u, id_u, level_v, id_v, d in distances:
            f.write(f"{level_u} {id_u} {level_v} {id_v} {d:.6f}\n")


def default_output_path(instance_path: str, technique: str,
                        pca_dimension: int = DEFAULT_PCA_DIM,
                        metric: str = DEFAULT_METRIC) -> Path:
    DEFAULT_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Sufixos so aparecem quando o valor difere do default, para nao quebrar
    # os nomes ja usados nas execucoes anteriores.
    suffix = ""
    if pca_dimension != DEFAULT_PCA_DIM:
        suffix += f".pca{pca_dimension}"
    if metric != DEFAULT_METRIC:
        suffix += f".{metric}"

    return DEFAULT_OUTPUT_DIR / f"{Path(instance_path).stem}.{technique}{suffix}.dist.txt"


def run_one(instance, instance_path: str, technique: str, output_path: Path, dimension, params,
            pca_dimension: int = DEFAULT_PCA_DIM, metric: str = DEFAULT_METRIC) -> None:
    t0 = time.time()
    embeddings = compute_embeddings(instance, technique, dimension=dimension, **params)
    dim = next(iter(embeddings.values())).shape[0]

    coords_nd = project_to_nd(embeddings, n=pca_dimension)
    distances = compute_arc_distances(instance, coords_nd, metric)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    write_distances(output_path, distances)

    print(
        f"[{technique}] dim={dim} pca={pca_dimension} metrica={metric} "
        f"arcos={len(distances)} tempo={time.time() - t0:.1f}s -> {output_path}"
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
        "--pca-dim", type=int, default=DEFAULT_PCA_DIM,
        help=f"dimensao da projecao PCA aplicada sobre o embedding (default: {DEFAULT_PCA_DIM}, como no artigo)"
    )
    parser.add_argument(
        "--metric", default=DEFAULT_METRIC, choices=METRICS,
        help=f"metrica de distancia entre vertices adjacentes (default: {DEFAULT_METRIC})"
    )
    parser.add_argument(
        "--param", action="append", metavar="k=v", help="sobrescreve um parametro especifico da tecnica (repetivel)"
    )
    args = parser.parse_args()
    params = parse_params(args.param)

    if args.pca_dim < 1:
        parser.error("--pca-dim deve ser um inteiro positivo")

    instance = read_instance(args.instance)
    print(f"instancia: {instance.num_levels} niveis, {instance.num_nodes} vertices, {len(instance.arcs)} arcos")

    techniques = TECHNIQUES if args.technique == "all" else [args.technique]

    failures = []
    for technique in techniques:
        output_path = (
            Path(args.output)
            if (args.output and args.technique != "all")
            else default_output_path(args.instance, technique, args.pca_dim, args.metric)
        )
        try:
            run_one(instance, args.instance, technique, output_path, args.dim, params,
                    args.pca_dim, args.metric)
        except Exception as exc:  # noqa: BLE001 - queremos continuar as outras tecnicas mesmo se uma falhar
            print(f"[{technique}] FALHOU: {exc}")
            traceback.print_exc()
            failures.append(technique)

    if len(techniques) > 1:
        ok = len(techniques) - len(failures)
        print(f"\nresumo: {ok}/{len(techniques)} tecnicas concluidas" + (f", falharam: {failures}" if failures else ""))


if __name__ == "__main__":
    main()
