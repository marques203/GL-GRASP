"""Retoma o experimento das 8 tecnicas de onde parou.

Le results_all/runs.csv, descobre quais combinacoes (tecnica, instancia, k,
semente) ainda nao foram executadas e roda apenas essas, anexando ao mesmo
arquivo. E idempotente: pode ser reexecutado quantas vezes for preciso, o que
importa se a maquina for suspensa ou reiniciada no meio (foi o que derrubou a
primeira tentativa).

Alem do resumo em runs.csv, guarda o historico completo de cada execucao em
results_all/traj/<tecnica>/<instancia>_k<k>_s<semente>.txt, que e o que
permite calcular o NPI (primal integral) e comparar os resultados pos fase
construtiva com os pos busca local. Formato, uma linha por iteracao do GRASP:

    phi  semente  custo_construcao  custo_busca_local  tempo_acumulado_s

e a ULTIMA linha e um resumo no mesmo formato (melhor custo de construcao,
melhor custo final e tempo total), nao uma iteracao.

Uso: python resume_experiment.py
"""
import csv
import os
import subprocess
import sys

BIN = "./C-IGDP.exe"
DIST_DIR = "embeddings/distances"
RESULTS = "results_all/runs.csv"
TMP_BEST = "results_all/_resume_best.txt"
TRAJ_DIR = "results_all/traj"

ALPHA = "0"
LS = "best"
MAXIT = "100"
TIME_LIMIT = "300"
ETA_MAX = "20"
REPS = 30

# Leitura do Step 5 e tratamento de G indefinido (ver README). Explicitos aqui
# para que o registro do experimento nao dependa dos defaults do binario.
POSITION_RULE = "figure"
UNANCHORED_RULE = "maxdis"
FEASIBLE_SEARCH = "reference"

# As 8 tecnicas. O hope ficava de fora porque vivia em
# results_wilcoxon/runs.csv, mas aqueles resultados sao do codigo anterior ao
# alinhamento com o codigo de referencia e nao tem historico por iteracao.
# Os resultados antigos destas 7 tecnicas estao em
# results_all/runs_pre_alinhamento.csv.
TECHNIQUES = ["hope", "spectral", "node2vec", "sdne", "line", "grarep", "netsmf", "prone"]

INSTANCES = [
    f"incgraph_{lev}_{dens}_5_30_{ratio}_1"
    for lev in (6, 13, 20)
    for dens in ("0.06", "0.17", "0.30")
    for ratio in ("1.20", "1.60")
]


def traj_path(technique, inst, k, seed):
    """Arquivo com o historico por iteracao de uma execucao."""
    return f"{TRAJ_DIR}/{technique}/{inst}_k{k}_s{seed}.txt"


def already_done():
    done = set()
    if not os.path.exists(RESULTS):
        return done
    with open(RESULTS) as f:
        for r in csv.DictReader(f):
            done.add((r["algo"], r["instance"], r["k"], int(r["seed"])))
    return done


def main():
    done = already_done()

    todo = [
        (t, inst, k, seed)
        for t in TECHNIQUES
        for inst in INSTANCES
        for k in ("1", "2", "3")
        for seed in range(1, REPS + 1)
        if (t, inst, k, seed) not in done
    ]

    print(f"ja executadas: {len(done)} | faltando: {len(todo)}")
    if not todo:
        print("nada a fazer")
        return

    new_file = not os.path.exists(RESULTS)
    out = open(RESULTS, "a", newline="")
    writer = csv.writer(out)
    if new_file:
        writer.writerow(["instance", "k", "algo", "seed", "best_cost", "total_time_s"])

    last_group = None
    for n, (t, inst, k, seed) in enumerate(todo, 1):

        if (t, inst, k) != last_group:
            last_group = (t, inst, k)
            print(f"  [{n}/{len(todo)}] {t} {inst} k={k}", flush=True)

        dist = f"{DIST_DIR}/{inst}.{t}.dist.txt"

        # O historico vai direto para o arquivo definitivo. O binario abre esse
        # arquivo em modo append, logo uma sobra de execucao interrompida tem
        # de ser apagada antes, senao as linhas se acumulariam.
        traj = traj_path(t, inst, k, seed)
        os.makedirs(os.path.dirname(traj), exist_ok=True)
        for p in (TMP_BEST, traj):
            if os.path.exists(p):
                os.remove(p)

        subprocess.run(
            [BIN, f"instance/{inst}.txt", "glgrasp", k, ALPHA, LS, MAXIT, inst,
             TMP_BEST, traj, TIME_LIMIT, dist, ETA_MAX, str(seed),
             POSITION_RULE, UNANCHORED_RULE, FEASIBLE_SEARCH],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )

        if os.path.exists(TMP_BEST):
            parts = open(TMP_BEST).read().split()
            # formato: k alpha custo tempo_total tempo_ate_melhor
            writer.writerow([inst, k, t, seed, parts[-3], parts[-2]])
        else:
            # main.cpp pula quando k excede o numero minimo de nos incrementais
            writer.writerow([inst, k, t, seed, "", "SKIPPED"])

        out.flush()   # grava a cada execucao: se cair, nao se perde nada

    out.close()
    print("concluido")


if __name__ == "__main__":
    main()
