"""Analise A aplicada as 8 tecnicas de embedding.

Junta os resultados de experiment_wilcoxon.sh (GRASP2 e hope) com os de
experiment_all_techniques.sh (as outras 7 tecnicas), agrega as 30 repeticoes
por grupo (instancia x k) e aplica o teste de Wilcoxon pareado de cada
tecnica contra o GRASP2 -- a mesma metodologia de Charytitsch e Nascimento
(2026).

Gera as tabelas em LaTeX prontas para o TG.

Uso: python analyze_all_techniques.py > docs/tabelas_geradas.tex
"""
import collections
import csv
import re
import statistics as st
import sys

from scipy.stats import wilcoxon

TECHNIQUES = ["spectral", "hope", "node2vec", "sdne", "line", "grarep", "netsmf", "prone"]
BASELINE = "grasp2"


def load():
    """Devolve runs[(instancia, k, algo)] -> lista de (custo, tempo)."""
    runs = collections.defaultdict(list)
    for path in ("results_wilcoxon/runs.csv", "results_all/runs.csv"):
        try:
            f = open(path)
        except FileNotFoundError:
            print(f"% AVISO: {path} nao encontrado", file=sys.stderr)
            continue
        for r in csv.DictReader(f):
            if not r["best_cost"]:
                continue
            # experiment_wilcoxon.sh grava a tecnica hope sob o rotulo "glgrasp"
            algo = "hope" if r["algo"] == "glgrasp" else r["algo"]
            runs[(r["instance"], r["k"], algo)].append(
                (float(r["best_cost"]), float(r["total_time_s"]))
            )
    return runs


def means(runs):
    """Agrega as repeticoes: m[(instancia,k)][algo] = (custo medio, tempo medio, dp)."""
    m = collections.defaultdict(dict)
    for (inst, k, algo), vals in runs.items():
        costs = [v[0] for v in vals]
        times = [v[1] for v in vals]
        m[(inst, k)][algo] = (
            st.mean(costs),
            st.mean(times),
            st.stdev(costs) if len(costs) > 1 else 0.0,
            len(costs),
        )
    return m


def instance_key(inst):
    mt = re.match(r"incgraph_(\d+)_([\d.]+)_5_30_([\d.]+)_1", inst)
    return int(mt.group(1)), mt.group(2), mt.group(3)


def main():
    m = means(load())

    # grupos onde o baseline e TODAS as tecnicas tem resultado
    groups = sorted(
        (instance_key(i) + (int(k),), i, k)
        for (i, k) in m
        if BASELINE in m[(i, k)] and all(t in m[(i, k)] for t in TECHNIQUES)
    )
    print(f"% grupos comparaveis: {len(groups)}")
    reps = {m[(i, k)][a][3] for _, i, k in groups for a in [BASELINE] + TECHNIQUES}
    print(f"% repeticoes por configuracao: {sorted(reps)}")
    print()

    # ------------------------------------------------------------------
    # Tabela 1: resumo por tecnica (a tabela principal do experimento)
    # ------------------------------------------------------------------
    base_c = [m[(i, k)][BASELINE][0] for _, i, k in groups]
    base_t = [m[(i, k)][BASELINE][1] for _, i, k in groups]

    rows = []
    for t in TECHNIQUES:
        cs = [m[(i, k)][t][0] for _, i, k in groups]
        ts = [m[(i, k)][t][1] for _, i, k in groups]
        gaps = [(c - b) / b * 100 for c, b in zip(cs, base_c)]
        wins = sum(1 for c, b in zip(cs, base_c) if c < b)
        ties = sum(1 for c, b in zip(cs, base_c) if c == b)
        loss = len(cs) - wins - ties
        _, p_cost = wilcoxon(base_c, cs)
        _, p_time = wilcoxon(base_t, ts)
        rows.append((st.mean(gaps), t, st.median(gaps), wins, ties, loss,
                     p_cost, sum(base_t) / sum(ts), p_time))

    print(r"\begin{table}[htbp]")
    print(r"\centering\small")
    print(r"\caption{Desempenho das oito t\'ecnicas de \textit{embedding} no GL-GRASP, "
          r"comparadas ao GRASP2. Cada valor agrega 30 repeti\c{c}\~oes por grupo "
          r"(inst\^ancia $\times$ $k$); $p$ vem do teste de Wilcoxon pareado sobre os "
          f"{len(groups)} grupos. Gap positivo indica mais cruzamentos que o GRASP2 "
          r"(pior); \textit{speedup} maior que 1 indica execu\c{c}\~ao mais r\'apida.}")
    print(r"\label{tab:tecnicas}")
    print(r"\begin{tabular}{lrrrrrrr}")
    print(r"\toprule")
    print(r"T\'ecnica & Gap m\'edio & Gap mediano & V/E/D & $p$ (custo) & "
          r"\textit{Speedup} & $p$ (tempo) \\")
    print(r"\midrule")
    for mg, t, md, w, e, l, pc, sp, pt in sorted(rows):
        star_c = r"$^{\dagger}$" if pc < 0.05 else ""
        star_t = r"$^{\dagger}$" if pt < 0.05 else ""
        print(f"\\texttt{{{t}}} & {mg:+.4f}\\% & {md:+.4f}\\% & {w}/{e}/{l} & "
              f"{pc:.4f}{star_c} & {sp:.2f}$\\times$ & {pt:.4f}{star_t} \\\\")
    print(r"\bottomrule")
    print(r"\end{tabular}")
    print()
    print(r"\vspace{2mm}")
    print(r"\footnotesize V/E/D: vit\'orias/empates/derrotas contra o GRASP2. "
          r"$^{\dagger}$ diferen\c{c}a estatisticamente significativa ao n\'ivel de 5\%.")
    print(r"\end{table}")
    print()

    # ------------------------------------------------------------------
    # Tabela 2: recortes por caracteristica da instancia
    # ------------------------------------------------------------------
    cuts = [
        ("N\\'iveis", lambda key: key[0]),
        ("Densidade", lambda key: key[1]),
        ("$k$", lambda key: key[3]),
    ]

    print(r"\begin{table}[htbp]")
    print(r"\centering\small")
    print(r"\caption{Gap m\'edio (\%) de cada t\'ecnica em rela\c{c}\~ao ao GRASP2, "
          r"por caracter\'istica da inst\^ancia.}")
    print(r"\label{tab:recortes}")
    print(r"\begin{tabular}{ll" + "r" * len(TECHNIQUES) + "}")
    print(r"\toprule")
    print(r"Recorte & Valor & " + " & ".join(f"\\texttt{{{t}}}" for t in TECHNIQUES) + r" \\")
    print(r"\midrule")
    for label, fn in cuts:
        buckets = collections.defaultdict(list)
        for key, i, k in groups:
            buckets[fn(key)].append((i, k))
        for first, val in enumerate(sorted(buckets)):
            cells = []
            for t in TECHNIQUES:
                gaps = [
                    (m[(i, k)][t][0] - m[(i, k)][BASELINE][0]) / m[(i, k)][BASELINE][0] * 100
                    for (i, k) in buckets[val]
                ]
                cells.append(f"{st.mean(gaps):+.3f}")
            name = label if first == 0 else ""
            print(f"{name} & {val} & " + " & ".join(cells) + r" \\")
        print(r"\midrule")
    print(r"\bottomrule")
    print(r"\end{tabular}")
    print(r"\end{table}")


if __name__ == "__main__":
    main()
