"""Geracao de node embeddings para o GL-GRASP, com as 8 tecnicas da tabela do TG
(4 "Base", do artigo original, + 4 "Nova", extensao proposta na pesquisa):

    Base: spectral, hope, node2vec, sdne
    Nova: line, grarep, netsmf, prone

Todas usam a implementacao de referencia do CogDL (cogdl.models.emb.*),
apenas adaptando a instancia do C-IGDP (niveis + ids locais) para o Graph
que o CogDL espera, igual ja era feito soh para o HOPE.
"""
import numpy as np

# cogdl 0.6 (node2vec/line/netsmf, via cogdl.utils.alias_setup) usa o alias
# np.int, removido a partir do numpy 1.24. Sem isso, essas tecnicas quebram
# com AttributeError ao serem chamadas.
if not hasattr(np, "int"):
    np.int = int

import torch
from cogdl.data import Graph
from cogdl.models.emb.spectral import Spectral
from cogdl.models.emb.hope import HOPE
from cogdl.models.emb.node2vec import Node2vec
from cogdl.models.emb.sdne import SDNE
from cogdl.models.emb.line import LINE
from cogdl.models.emb.grarep import GraRep
from cogdl.models.emb.netsmf import NetSMF
from cogdl.models.emb.prone import ProNE

from hdag_io import Instance

TECHNIQUES = ["spectral", "hope", "node2vec", "sdne", "line", "grarep", "netsmf", "prone"]


def default_dimension(n: int, technique: str) -> int:
    """Dimensao default por tecnica, seguindo o artigo (Secao 5.2) quando aplicavel."""
    if technique == "spectral":
        return n - 1  # SPEC: k = n - 1 (usa so os vetores singulares esquerdos)
    if technique == "hope":
        return 2 * n - 1  # HOPE: k = 2n - 1 (esquerdo + direito)
    return min(128, max(2, n - 1))  # default do CogDL para os demais, limitado pelo tamanho do grafo


def _build_model(technique: str, dimension: int, params: dict):
    """Monta o modelo do CogDL para a tecnica pedida.

    `params` permite sobrescrever qualquer parametro especifico da tecnica
    (ver --param na CLI de build_distances.py); os defaults abaixo foram
    reduzidos em relacao ao CogDL (pensado para grafos de milhoes de nos)
    para rodar em segundos/minutos nas instancias do C-IGDP (centenas de nos).
    """
    p = dict(params)

    if technique == "spectral":
        return Spectral(dimension)

    if technique == "hope":
        return HOPE(dimension, p.pop("beta", 0.01))

    if technique == "node2vec":
        return Node2vec(
            dimension,
            p.pop("walk_length", 40),
            p.pop("walk_num", 10),
            p.pop("window_size", 5),
            p.pop("worker", 4),
            p.pop("iteration", 5),
            p.pop("p", 1.0),
            p.pop("q", 1.0),
        )

    if technique == "sdne":
        hidden_size1 = p.pop("hidden_size1", min(256, max(dimension * 2, 8)))
        return SDNE(
            hidden_size1,
            dimension,
            p.pop("dropout", 0.2),
            p.pop("alpha", 1e-1),
            p.pop("beta", 5.0),
            p.pop("nu1", 1e-4),
            p.pop("nu2", 1e-3),
            p.pop("epochs", 20),
            p.pop("lr", 5e-3),
            True,  # cpu
        )

    if technique == "line":
        return LINE(
            dimension,
            p.pop("walk_length", 40),
            p.pop("walk_num", 10),
            p.pop("negative", 5),
            p.pop("batch_size", 128),
            p.pop("alpha", 0.025),
            p.pop("order", 3),
        )

    if technique == "grarep":
        return GraRep(dimension, p.pop("step", 5))

    if technique == "netsmf":
        return NetSMF(
            dimension,
            p.pop("window_size", 10),
            p.pop("negative", 1),
            p.pop("num_round", 10),
            p.pop("worker", 4),
        )

    if technique == "prone":
        return ProNE(dimension, p.pop("step", 5), p.pop("mu", 0.2), p.pop("theta", 0.5))

    raise ValueError(f"tecnica desconhecida: {technique!r} (opcoes: {TECHNIQUES})")


def _build_cogdl_graph(instance: Instance):
    """Mapeia (level, id_local) -> indice global 0..n-1 e monta o Graph do CogDL."""
    nodes = list(instance.all_nodes())
    index_of = {key: i for i, key in enumerate(nodes)}

    edges = [(index_of[u], index_of[v]) for u, v in instance.arcs]
    edge_index = torch.tensor(edges, dtype=torch.long).t().contiguous()

    graph = Graph(edge_index=edge_index)
    graph.num_nodes = len(nodes)

    return graph, nodes


def compute_embeddings(instance: Instance, technique: str, dimension: int = None, **params) -> dict:
    """Retorna {(level, id_local): np.ndarray} com o embedding de cada vertice de IG,
    usando a tecnica pedida (ver TECHNIQUES)."""
    graph, nodes = _build_cogdl_graph(instance)

    if dimension is None:
        dimension = default_dimension(instance.num_nodes, technique)

    model = _build_model(technique, dimension, params)
    features = model.forward(graph)

    return {key: features[i] for i, key in enumerate(nodes)}
