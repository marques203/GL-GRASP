"""Gera instancias sinteticas para o teste T4 de tests/validate_changes.py.

Nenhuma das 18 instancias reais serve para T4: todas tem pelo menos um
vertice incremental sem vizinho original, logo a flag do item 2.1 esta sempre
ativa nelas e nao da para testar o ESCOPO da flag. As duas fixtures abaixo
cobrem os dois lados:

  synth_anchored    : todo incremental tem ao menos um vizinho ORIGINAL, logo
                      G(u) nunca fica indefinido e "inf" tem de ser
                      identico a "maxdis".
  synth_unanchored  : mesma estrutura, mas alguns incrementais tem vizinhos
                      somente incrementais, logo G(u) fica indefinido e as
                      duas semanticas tem de divergir. Serve de controle
                      positivo: mostra que T4 tem poder de deteccao.

Uso:  python tests/make_fixtures.py
"""
import math
import os
import random

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_INST = os.path.join(HERE, "fixtures")
OUT_DIST = os.path.join(HERE, "fixtures")

N_LEVELS = 5
N_ORIG = 6
N_INCR = 4
LEVEL_SIZE = N_ORIG + N_INCR


def build(unanchored_count, isolate=False):
    """Monta (linhas_da_instancia, arcos).

    Identidades: em cada nivel, os ids 0..N_ORIG-1 sao originais e
    N_ORIG..LEVEL_SIZE-1 sao incrementais. A ORDEM DAS LINHAS dentro de cada
    nivel e pi_0, e os originais precedem os incrementais, como nas instancias
    reais. Embaralhamos os ids dos originais para que pi_0 nao seja a
    identidade (e o caso geral).
    """
    rnd = random.Random(20260106)

    orig_ids = list(range(N_ORIG))
    incr_ids = list(range(N_ORIG, LEVEL_SIZE))

    # ordem de pi_0 por nivel
    order = []
    for l in range(N_LEVELS):
        o = orig_ids[:]
        rnd.shuffle(o)
        order.append(o + incr_ids[:])

    # os incrementais que ficarao SEM vizinho original
    deprived = set()
    for l in range(1, N_LEVELS - 1):
        for j in range(unanchored_count):
            if j < N_INCR:
                deprived.add((l, N_ORIG + j))

    # arcos como dict (l, id_u) -> lista de id_v no nivel l+1
    fwd = {(l, i): [] for l in range(N_LEVELS) for i in range(LEVEL_SIZE)}

    def add(l, u, v):
        if l + 1 < N_LEVELS and v not in fwd[(l, u)]:
            fwd[(l, u)].append(v)

    # 1) malha entre originais, suficiente para gerar cruzamentos
    for l in range(N_LEVELS - 1):
        for u in orig_ids:
            add(l, u, (u + 1) % N_ORIG)
            add(l, u, (u * 2 + 3) % N_ORIG)

    # 2) cada incremental recebe vizinhos; quem esta em "deprived" so se
    #    conecta a outros INCREMENTAIS, para que G(u) fique indefinido.
    for l in range(N_LEVELS):
        for u in incr_ids:

            if (l, u) in deprived:
                pool = incr_ids
            else:
                pool = orig_ids

            # arco para frente
            if l + 1 < N_LEVELS:
                add(l, u, pool[(u + l) % len(pool)])

            # arco para tras (vizinho no nivel l-1 aponta para u)
            if l > 0:
                w = pool[(u + l + 1) % len(pool)]
                add(l - 1, w, u)

    # 3) os "deprived" nao podem receber arco de nenhum original
    for (l, u) in deprived:
        if l > 0:
            for w in orig_ids:
                if u in fwd[(l - 1, w)]:
                    fwd[(l - 1, w)].remove(u)

    # 4) um incremental TOTALMENTE isolado (grau zero). Um vertice que tem
    #    algum vizinho acaba ancorado quando esse vizinho e posicionado; so o
    #    de grau zero mantem G(u) indefinido ate o fim, que e o caso em que
    #    position_without_anchor realmente decide a posicao.
    isolated = None
    if isolate:
        l_iso, u_iso = N_LEVELS // 2, LEVEL_SIZE - 1
        fwd[(l_iso, u_iso)] = []
        for w in range(LEVEL_SIZE):
            if l_iso > 0 and u_iso in fwd[(l_iso - 1, w)]:
                fwd[(l_iso - 1, w)].remove(u_iso)
        isolated = (l_iso, u_iso)

    lines = [str(N_LEVELS), " ".join(str(LEVEL_SIZE) for _ in range(N_LEVELS))]

    arcs = []
    for l in range(N_LEVELS):
        for u in order[l]:
            flag = 1 if u < N_ORIG else 0
            nbs = sorted(fwd[(l, u)])
            lines.append(" ".join([str(flag), str(u)] + [str(v) for v in nbs]))
            for v in nbs:
                arcs.append((l, u, l + 1, v))

    return lines, arcs, deprived, isolated


def write_dist(path, arcs, outlier=None):
    """Distancias deterministicas e bem espalhadas em (0, 1].

    Com "outlier", UM arco recebe uma distancia muito maior que as demais.
    Isso isola o mecanismo da secao 2.1: com a semantica "maxdis" esse valor
    passa a ser o max G da Eq. (12) sempre que houver candidato sem ancora,
    inflando xi e engordando a RCL; com "inf" ele e so mais um arco.
    """
    with open(path, "w") as f:
        f.write("%d\n" % len(arcs))
        for n, (lu, iu, lv, iv) in enumerate(arcs):
            if outlier is not None and n == 0:
                d = outlier
            else:
                # valores distintos e nao monotonos, para que min e max de G
                # nao coincidam trivialmente
                d = 0.05 + 0.9 * abs(math.sin(1.7 * n + 0.3 * lu + 0.11 * iu))
            f.write("%u %u %u %u %.6f\n" % (lu, iu, lv, iv, d))


def main():
    os.makedirs(OUT_INST, exist_ok=True)

    specs = (
        # nome,              deprived/nivel, isolado, outlier de distancia
        ("synth_anchored",   0,              False,   None),
        ("synth_unanchored", 2,              True,    50.0),
    )

    for name, n_dep, isolate, outlier in specs:

        lines, arcs, deprived, isolated = build(n_dep, isolate)

        inst_path = os.path.join(OUT_INST, name + ".txt")
        with open(inst_path, "w") as f:
            f.write("\n".join(lines) + "\n")

        write_dist(os.path.join(OUT_DIST, name + ".hope.dist.txt"), arcs, outlier)

        # confere a propriedade que a fixture deve ter
        orig = set(range(N_ORIG))
        anchored = set()
        for (lu, iu, lv, iv) in arcs:
            if iv in orig:
                anchored.add((lu, iu))
            if iu in orig:
                anchored.add((lv, iv))

        n_unanch = sum(
            1 for l in range(N_LEVELS) for i in range(N_ORIG, LEVEL_SIZE)
            if (l, i) not in anchored)

        print("%-20s %d niveis, %d vertices/nivel, %d arcos, "
              "%d incrementais sem ancora, isolado=%s, max_dis=%.2f"
              % (name, N_LEVELS, LEVEL_SIZE, len(arcs), n_unanch,
                 isolated, outlier if outlier else 0.95))

        if n_dep == 0 and n_unanch != 0:
            print("    AVISO: a fixture 'anchored' nao devia ter vertice sem ancora")
        if n_dep > 0 and n_unanch == 0:
            print("    AVISO: a fixture 'unanchored' devia ter vertice sem ancora")


if __name__ == "__main__":
    main()
