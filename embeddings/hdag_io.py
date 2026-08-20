"""Leitura do formato de instancia do C-IGDP (mesmo formato lido por HDAG::read_instance em HDAG.cpp).

Formato do arquivo:
  linha 1: L (numero de niveis)
  linha 2: L inteiros, tamanho de cada nivel
  em seguida, L blocos de linhas (um bloco por nivel, na ordem):
    cada linha: o u [v1 v2 ...]
      o: 1 se vertice original, 0 se incremental
      u: id do vertice, LOCAL ao nivel (0..tamanho_do_nivel-1), nao eh a
         mesma coisa que a posicao (ordem) da linha no bloco
      v1 v2 ...: (apenas para niveis 0..L-2) ids locais dos vertices do
         PROXIMO nivel aos quais u tem arco
"""
from dataclasses import dataclass, field


@dataclass(frozen=True)
class NodeKey:
    level: int
    node_id: int


@dataclass
class Instance:
    num_levels: int
    level_sizes: list
    # (level, id) -> flag (1 = original, 0 = incremental)
    node_flag: dict
    # lista de arcos ((level_u, id_u), (level_v, id_v)); sempre level_v == level_u + 1
    arcs: list = field(default_factory=list)

    def all_nodes(self):
        for level, size in enumerate(self.level_sizes):
            for node_id in range(size):
                yield (level, node_id)

    @property
    def num_nodes(self):
        return sum(self.level_sizes)


def read_instance(path: str) -> Instance:
    with open(path, "r") as f:
        lines = [line for line in f]

    num_levels = int(lines[0].split()[0])
    level_sizes = [int(x) for x in lines[1].split()[:num_levels]]

    node_flag = {}
    arcs = []

    row = 2
    for level in range(num_levels):
        for _ in range(level_sizes[level]):
            parts = lines[row].split()
            row += 1
            o = int(parts[0])
            u = int(parts[1])
            node_flag[(level, u)] = o
            if level < num_levels - 1:
                for tok in parts[2:]:
                    v = int(tok)
                    arcs.append(((level, u), (level + 1, v)))

    return Instance(num_levels=num_levels, level_sizes=level_sizes, node_flag=node_flag, arcs=arcs)
