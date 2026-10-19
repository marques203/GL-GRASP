"""Validacao de integracao das mudancas dos Tiers 1 e 2 no GL-GRASP.

Itens cobertos:
  1.1 sorteio uniforme na RCL          -> T5 (invariantes), unidade em C++
  1.2 regra de posicao do Step 5       -> T3 (a flag muda o resultado)
  1.3 semente da construcao por rep.   -> T2 (repeticoes agora independentes)
  1.4 limite inferior da Claim 2       -> T5 (nenhuma solucao inviavel)
  2.1 G indefinido na Eq. (12)         -> T3 e T4 (a flag muda o resultado, e
                                          SO em instancias com vertice sem
                                          ancora)

Uso:
    python tests/validate_changes.py
    python tests/validate_changes.py --quick     (subconjunto menor)
"""
import argparse
import csv
import itertools
import os
import re
import subprocess
import sys
import shutil
import tempfile
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(ROOT, "C-IGDP.exe")
INST_DIR = os.path.join(ROOT, "instance")
DIST_DIR = os.path.join(ROOT, "embeddings", "distances")

# Os arquivos temporarios NAO podem ficar dentro do projeto: a pasta esta no
# OneDrive e o sync trava o arquivo no meio da bateria (WinError 32).
TMP = tempfile.mkdtemp(prefix="glgrasp_validate_")

_run_counter = itertools.count()

POSITION_RULES = ["figure", "text"]
UNANCHORED_RULES = ["maxdis", "inf"]
SEARCH_RULES = ["reference", "spiral"]

failures = []
checks = 0


def check(ok, name, detail=""):
    global checks
    checks += 1
    print("  [%s] %s %s" % ("ok  " if ok else "FAIL", name, detail))
    if not ok:
        failures.append(name)


# ---------------------------------------------------------------------------
# Execucao do binario
# ---------------------------------------------------------------------------
def run(instance, k, seed, pos_rule="figure", unanch="maxdis",
        search="reference", alpha="0", maxit="100", ls="best", eta_max="20",
        timelimit="300", base_dir=None, construction=False):
    """Roda o GL-GRASP e devolve (custo_inicial, melhor_custo, stderr, rc).

    base_dir permite apontar para tests/fixtures em vez de instance/.

    Com construction=True, o segundo elemento passa a ser o custo da PRIMEIRA
    construcao (primeira linha do arquivo completo) em vez do melhor custo.
    A diferenca importa: o melhor custo e min(custo do desenho original da
    instancia, custos encontrados), e numa instancia em que a construcao
    raramente bate o desenho original ele fica travado no custo da instancia,
    escondendo qualquer variacao da construcao.
    """
    os.makedirs(TMP, exist_ok=True)
    # nome unico por execucao: o binario ABRE em modo append, logo reusar o
    # mesmo arquivo acumularia linhas de execucoes anteriores
    tag = next(_run_counter)
    best = os.path.join(TMP, "best_%d.txt" % tag)
    full = os.path.join(TMP, "full_%d.txt" % tag)

    inst_dir = base_dir or INST_DIR
    dist_dir = base_dir or DIST_DIR

    cmd = [BIN, os.path.join(inst_dir, instance + ".txt"), "glgrasp",
           str(k), alpha, ls, maxit, instance, best, full, timelimit,
           os.path.join(dist_dir, instance + ".hope.dist.txt"),
           eta_max, str(seed), pos_rule, unanch, search]

    pr = subprocess.run(cmd, capture_output=True, text=True)

    initial = None
    m = re.search(r"instance cost\s+(\d+)", pr.stdout)
    if m:
        initial = int(m.group(1))

    cost = None
    if os.path.exists(best):
        parts = open(best).read().split()
        if len(parts) >= 3:
            cost = float(parts[-3])
        os.remove(best)
    if os.path.exists(full):
        if construction:
            lines = [ln.split() for ln in open(full).read().splitlines() if ln.strip()]
            # formato: phi  semente  custo_construcao  custo_busca_local
            # tempo; a ultima linha e o resumo final, as anteriores sao uma
            # por iteracao. A coluna 2 e o custo da construcao.
            cost = float(lines[0][2]) if len(lines) >= 2 else None
        os.remove(full)

    return initial, cost, pr.stderr, pr.returncode


def min_incremental_per_level(instance, base_dir=None):
    """min_l (numero de incrementais no nivel l).

    main.cpp nao executa nada quando k excede esse valor (guard
    k <= MIN_INCREMENTAL_NUMBER), e nesse caso nao ha arquivo de saida. Os
    testes tem de PULAR essas combinacoes, nao contar como falha.
    """
    path = os.path.join(base_dir or INST_DIR, instance + ".txt")
    n_levels, sizes, originais, _ = read_instance(path)
    return min(sizes[l] - len(originais[l]) for l in range(n_levels))


# ---------------------------------------------------------------------------
# Quais instancias tem vertice incremental SEM NENHUM vizinho original?
# Esse e exactamente o caso em que G(u) fica indefinido, ou seja, em que a
# flag do item 2.1 pode mudar algo.
# ---------------------------------------------------------------------------
def read_instance(path):
    """Devolve (niveis, originais, arcos) com ids locais ao nivel.

    originais[l] = conjunto de ids originais do nivel l
    arcos = lista de (l_u, id_u, id_v) com v no nivel l_u + 1
    """
    tokens = open(path).read().split("\n")
    tokens = [t for t in tokens if t.strip() != ""]

    n_levels = int(tokens[0].split()[0])
    sizes = [int(x) for x in tokens[1].split()]

    originais = [set() for _ in range(n_levels)]
    arcos = []

    pos = 2
    for l in range(n_levels):
        for _ in range(sizes[l]):
            f = tokens[pos].split()
            pos += 1
            flag, node_id = int(f[0]), int(f[1])
            if flag == 1:
                originais[l].add(node_id)
            for v in f[2:]:
                arcos.append((l, node_id, int(v)))

    return n_levels, sizes, originais, arcos


def has_unanchored(instance):
    """True se existe incremental sem vizinho ORIGINAL em nivel adjacente."""
    n_levels, sizes, originais, arcos = read_instance(
        os.path.join(INST_DIR, instance + ".txt"))

    # vizinhos originais de cada (nivel, id)
    anchored = set()
    for (lu, iu, iv) in arcos:
        lv = lu + 1
        if iv in originais[lv]:
            anchored.add((lu, iu))
        if iu in originais[lu]:
            anchored.add((lv, iv))

    for l in range(n_levels):
        for i in range(sizes[l]):
            if i not in originais[l] and (l, i) not in anchored:
                return True
    return False


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()

    if not os.path.exists(BIN):
        print("binario nao encontrado: %s  (rode gmake)" % BIN)
        return 1

    available = sorted(
        os.path.basename(p)[:-len(".hope.dist.txt")]
        for p in os.listdir(DIST_DIR) if p.endswith(".hope.dist.txt"))
    if not available:
        print("nenhum arquivo .hope.dist.txt em %s" % DIST_DIR)
        return 1

    instances = available[:3] if args.quick else available
    ks = [1, 2] if args.quick else [1, 2, 3]
    seeds = list(range(1, 4)) if args.quick else list(range(1, 11))

    print("\n=== Validacao de integracao: GL-GRASP Tier 1 + Tier 2 ===")
    print("instancias: %d | k: %s | sementes: %s\n"
          % (len(instances), ks, seeds))

    # classifica as instancias
    unanch_map = {inst: has_unanchored(inst) for inst in instances}
    com = [i for i in instances if unanch_map[i]]
    sem = [i for i in instances if not unanch_map[i]]
    print("  instancias COM incremental sem ancora: %d" % len(com))
    print("  instancias SEM incremental sem ancora: %d" % len(sem))

    # combinacoes (instancia, k) que o binario de fato executa
    min_inc = {inst: min_incremental_per_level(inst) for inst in instances}
    pairs = [(inst, k) for inst in instances for k in ks if k <= min_inc[inst]]
    skipped = [(inst, k) for inst in instances for k in ks if k > min_inc[inst]]
    print("  pares (instancia, k) executaveis: %d | pulados pelo guard de k: %d\n"
          % (len(pairs), len(skipped)))

    all_stderr = []
    all_rc = []

    # -----------------------------------------------------------------
    print("T1. Determinismo (mesmos argumentos -> mesmo custo)")
    det_ok = True
    n = 0
    for inst, k in [p for p in pairs if p[0] in instances[:3]]:
        for pr, ur in itertools.product(POSITION_RULES, UNANCHORED_RULES):
            a = run(inst, k, 7, pr, ur)
            b = run(inst, k, 7, pr, ur)
            all_stderr += [a[2], b[2]]
            all_rc += [a[3], b[3]]
            n += 1
            if a[1] is None or a[1] != b[1]:
                det_ok = False
                print("        %s k=%d %s/%s: %s vs %s"
                      % (inst, k, pr, ur, a[1], b[1]))
    check(det_ok, "T1 reproducao exata com semente fixa",
          "(%d configuracoes)" % n)

    # -----------------------------------------------------------------
    print("\nT2. Item 1.3: repeticoes independentes")
    # alpha fixo (nao sorteado) + 1 iteracao + sem busca local isola a
    # semente da CONSTRUCAO. Com o codigo antigo (setseed(i % 100)) todas as
    # sementes davam o MESMO custo aqui.
    n_varied = 0
    for inst in instances[:5]:
        costs = set()
        for s in range(1, 11):
            _, c, se, rc = run(inst, 1, s, alpha="0.5", maxit="1", ls="no",
                               construction=True)
            all_stderr.append(se)
            all_rc.append(rc)
            costs.add(c)
        if len(costs) > 1:
            n_varied += 1
        else:
            print("        %s: custo constante em 10 sementes (%s)"
                  % (inst, costs))
    check(n_varied == len(instances[:5]),
          "T2 semente da construcao depende da repeticao",
          "(%d/%d instancias variaram)" % (n_varied, len(instances[:5])))

    # a mesma semente ainda tem de reproduzir
    r1 = run(instances[0], 1, 3, alpha="0.5", maxit="1", ls="no",
                               construction=True)
    r2 = run(instances[0], 1, 3, alpha="0.5", maxit="1", ls="no",
                               construction=True)
    check(r1[1] == r2[1] and r1[1] is not None,
          "T2 ainda reprodutivel com a mesma semente", "")

    # -----------------------------------------------------------------
    print("\nT3. As flags novas tem efeito observavel")
    diff_pos = 0
    total_pos = 0
    for inst, k in pairs:
        for s in seeds[:3]:
            a = run(inst, k, s, "figure", "inf")
            b = run(inst, k, s, "text", "inf")
            all_stderr += [a[2], b[2]]
            all_rc += [a[3], b[3]]
            total_pos += 1
            if a[1] != b[1]:
                diff_pos += 1
    check(diff_pos > 0, "T3 position-rule figure != text",
          "(difere em %d/%d execucoes)" % (diff_pos, total_pos))

    diff_un = 0
    total_un = 0
    for inst, k in [p for p in pairs if p[0] in com]:
        for s in seeds[:3]:
            a = run(inst, k, s, "figure", "inf")
            b = run(inst, k, s, "figure", "maxdis")
            all_stderr += [a[2], b[2]]
            all_rc += [a[3], b[3]]
            total_un += 1
            if a[1] != b[1]:
                diff_un += 1
    if com:
        check(diff_un > 0, "T3 unanchored-rule inf != maxdis",
              "(difere em %d/%d execucoes, so instancias com vertice sem ancora)"
              % (diff_un, total_un))
    else:
        print("  [skip] T3 unanchored-rule: nenhuma instancia tem vertice sem ancora")

    # -----------------------------------------------------------------
    print("\nT4. Escopo da flag 2.1, em fixtures sinteticas")
    # Nenhuma das 18 instancias reais serve aqui: TODAS tem vertice
    # incremental sem ancora, logo nelas a flag esta sempre ativa e nao da
    # para testar onde ela NAO deve agir. As fixtures de tests/fixtures
    # cobrem os dois lados. Usamos maxit=1 e ls=no para isolar a construcao,
    # que e onde a flag atua.
    fix_dir = os.path.join(ROOT, "tests", "fixtures")
    needed = ["synth_anchored", "synth_unanchored"]
    missing = [f for f in needed
               if not os.path.exists(os.path.join(fix_dir, f + ".txt"))]

    if missing:
        print("  [skip] T4: fixtures ausentes %s "
              "(rode python tests/make_fixtures.py)" % missing)
    else:
        fix_seeds = list(range(1, 21))

        n_diff = {}
        for fx in needed:
            d = 0
            for s in fix_seeds:
                a = run(fx, 1, s, "figure", "inf", maxit="1", ls="no",
                        base_dir=fix_dir)
                b = run(fx, 1, s, "figure", "maxdis", maxit="1", ls="no",
                        base_dir=fix_dir)
                all_stderr += [a[2], b[2]]
                all_rc += [a[3], b[3]]
                if a[1] != b[1]:
                    d += 1
            n_diff[fx] = d

        check(n_diff["synth_anchored"] == 0,
              "T4(a) inf == maxdis sem vertice sem ancora",
              "(%d/%d sementes divergiram)"
              % (n_diff["synth_anchored"], len(fix_seeds)))

        # controle positivo: sem isto, T4(a) poderia passar por vacuidade
        check(n_diff["synth_unanchored"] > 0,
              "T4(b) contraste: inf != maxdis com vertice sem ancora",
              "(%d/%d sementes divergiram)"
              % (n_diff["synth_unanchored"], len(fix_seeds)))

    # -----------------------------------------------------------------
    print("\nT4b. Busca de posicao viavel: reference (g_vq) vs spiral")
    # maxit=1 e ls=no isolam a construcao: com a busca local ligada, as duas
    # regras costumam convergir para o mesmo otimo local nas instancias mais
    # faceis, o que esconderia a diferenca.
    diff_s = 0
    total_s = 0
    for inst, k in pairs[:12]:
        for s in seeds[:5]:
            a = run(inst, k, s, search="reference", maxit="1", ls="no")
            b = run(inst, k, s, search="spiral", maxit="1", ls="no")
            all_stderr += [a[2], b[2]]
            all_rc += [a[3], b[3]]
            total_s += 1
            if a[1] != b[1]:
                diff_s += 1
    check(diff_s > 0, "T4b feasible-search reference != spiral",
          "(difere em %d/%d construcoes)" % (diff_s, total_s))

    # o default do binario tem de ser a configuracao de referencia
    inst0, k0 = pairs[0]
    d_implicit = run(inst0, k0, 7)
    d_explicit = run(inst0, k0, 7, "figure", "maxdis", "reference")
    check(d_implicit[1] is not None and d_implicit[1] == d_explicit[1],
          "T4b default == 'figure maxdis reference'",
          "(%s vs %s)" % (d_implicit[1], d_explicit[1]))

    print("\nT5. Invariantes em toda a bateria")
    bad_stderr = [s for s in all_stderr if "NOT FEASIBLE" in s]
    check(not bad_stderr, "T5 nenhuma solucao inviavel",
          "(%d execucoes verificadas)" % len(all_stderr))

    other_err = sorted(set(s.strip() for s in all_stderr if s.strip()))
    check(not other_err, "T5 stderr vazio",
          "" if not other_err else "(%s)" % other_err[:3])

    check(all(rc == 0 for rc in all_rc), "T5 codigo de retorno 0 em todas",
          "(%d execucoes)" % len(all_rc))

    # -----------------------------------------------------------------
    print("\nT6. Sanidade do custo")
    sane = True
    n = 0
    for inst, k in pairs:
        init, cost, se, rc = run(inst, k, 1)
        all_stderr.append(se)
        n += 1
        if init is None or cost is None or not (0 <= cost <= init):
            sane = False
            print("        %s k=%d: inicial=%s final=%s" % (inst, k, init, cost))
    check(sane, "T6 0 <= custo final <= custo da instancia", "(%d execucoes)" % n)

    # O guard de k tem de se comportar do jeito documentado: nada executado.
    guard_ok = True
    for inst, k in skipped:
        init, cost, se, rc = run(inst, k, 1)
        if cost is not None or rc != 0:
            guard_ok = False
            print("        %s k=%d deveria ser pulado: custo=%s rc=%d"
                  % (inst, k, cost, rc))
    check(guard_ok, "T6 guard k > min incrementais nao produz saida",
          "(%d pares pulados)" % len(skipped))

    # -----------------------------------------------------------------
    print("\nT7. Qualidade: configuracao nova vs resultados armazenados (HOPE)")
    old_path = os.path.join(ROOT, "results_wilcoxon", "runs.csv")
    if os.path.exists(old_path):
        old = defaultdict(list)
        for r in csv.DictReader(open(old_path)):
            if r["algo"] == "glgrasp" and r["best_cost"]:
                old[(r["instance"], int(r["k"]))].append(float(r["best_cost"]))

        rows = []
        for inst, k in pairs:
            key = (inst, k)
            if key not in old:
                continue
            new_costs = []
            for s in seeds:
                _, c, se, rc = run(inst, k, s)
                all_stderr.append(se)
                if c is not None:
                    new_costs.append(c)
            if not new_costs:
                continue
            om = sum(old[key]) / len(old[key])
            nm = sum(new_costs) / len(new_costs)
            rows.append((inst, k, om, nm))

        if rows:
            better = sum(1 for _, _, om, nm in rows if nm < om)
            worse = sum(1 for _, _, om, nm in rows if nm > om)
            equal = len(rows) - better - worse
            rel = sum((nm - om) / om for _, _, om, nm in rows) / len(rows)

            print("    antigo = codigo pre-alinhamento | "
                  "novo = default atual (figure/maxdis/reference)\n")
            print("    %-42s %10s %10s %8s" % ("instancia", "k", "antigo", "novo"))
            for inst, k, om, nm in rows:
                flag = "" if nm == om else ("  melhor" if nm < om else "  pior")
                print("    %-42s %10d %10.1f %8.1f%s" % (inst, k, om, nm, flag))
            print("\n    melhor em %d, pior em %d, igual em %d (de %d pares)"
                  % (better, worse, equal, len(rows)))
            print("    variacao media do custo: %+.3f%%" % (100 * rel))
            print("    (medias: %d sementes no antigo, %d no novo -- e sanidade,"
                  " nao o experimento final)"
                  % (len(next(iter(old.values()))), len(seeds)))
            check(True, "T7 comparacao de qualidade registrada", "")
        else:
            print("  [skip] T7: sem pares comparaveis")
    else:
        print("  [skip] T7: %s nao existe" % old_path)

    # -----------------------------------------------------------------
    print("\n%d verificacoes, %d falha(s)" % (checks, len(failures)))
    for f in failures:
        print("  falhou: %s" % f)
    print()

    shutil.rmtree(TMP, ignore_errors=True)

    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
