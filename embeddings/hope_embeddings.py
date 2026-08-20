"""Geracao de node embeddings via HOPE (CogDL), conforme Secao 4.1 do artigo
(Charytitsch e Nascimento, 2026): decoder baseado no indice de Katz,
embeddings obtidos por SVD, preservando transitividade assimetrica.

Usa a implementacao de referencia cogdl.models.emb.hope.HOPE diretamente,
apenas adaptando a instancia (niveis + ids locais) para o formato de grafo
que o CogDL espera (Graph com edge_index 0-indexado contiguo).
"""
import numpy as np
import torch
from cogdl.data import Graph
from cogdl.models.emb.hope import HOPE

from hdag_io import Instance


def _build_cogdl_graph(instance: Instance):
    """Mapeia (level, id_local) -> indice global 0..n-1 e monta o Graph do CogDL."""
    nodes = list(instance.all_nodes())
    index_of = {key: i for i, key in enumerate(nodes)}

    edges = [(index_of[u], index_of[v]) for u, v in instance.arcs]
    edge_index = torch.tensor(edges, dtype=torch.long).t().contiguous()

    graph = Graph(edge_index=edge_index)
    graph.num_nodes = len(nodes)

    return graph, nodes


def compute_hope_embeddings(instance: Instance, dimension: int = None, beta: float = 0.01) -> dict:
    """Retorna {(level, id_local): np.ndarray} com o embedding HOPE de cada vertice de IG.

    dimension: dimensao do embedding (default segue o artigo, k = 2n - 1, onde
        n e o numero total de vertices de IG).
    beta: parametro de decaimento do indice de Katz (default do CogDL, 0.01).
    """
    graph, nodes = _build_cogdl_graph(instance)

    n = instance.num_nodes
    if dimension is None:
        dimension = 2 * n - 1

    model = HOPE(dimension=dimension, beta=beta)
    features = model.forward(graph)

    return {key: features[i] for i, key in enumerate(nodes)}
