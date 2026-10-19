#!/bin/bash
# compare_grasp2_glgrasp.sh
#
# Compara GRASP2 (baseline, Napoletano et al. 2019) com o GL-GRASP
# (Charytitsch e Nascimento 2026) sobre o mesmo conjunto de instancias.
#
# As execucoes sao SEQUENCIAIS de proposito: rodar os algoritmos ao mesmo
# tempo faria os processos disputarem CPU e contaminaria justamente a
# medicao de tempo que queremos comparar.
#
# Sao tres configuracoes por instancia/k:
#   grasp2        - baseline
#   glgrasp-e20   - GL-GRASP com eta_max = 20 (default do artigo)
#   glgrasp-e100  - GL-GRASP com eta_max = 100 = maxit, ou seja, criterio de
#                   parada por estagnacao DESLIGADO. Serve para comparar o
#                   custo por iteracao com o grasp2 em igualdade de numero
#                   de iteracoes (o e20 para antes, o que por si so ja
#                   reduz o tempo).
#
# Uso: ./compare_grasp2_glgrasp.sh

set -e

BIN=./C-IGDP.exe
DIST_DIR=embeddings/distances
TECHNIQUE=hope           # distancias euclidianas, PCA 2D
ALPHA=0
LS=best
MAXIT=100
TIME_LIMIT=300
KVALUES="1 2 3"

RESULTS_DIR=results_compare
mkdir -p "$RESULTS_DIR"

SUMMARY="$RESULTS_DIR/summary.csv"
echo "instance,k,algo,best_cost,total_time_s,time_to_best_s" > "$SUMMARY"

INSTANCES=$(for lev in 6 13 20; do
  for dens in 0.06 0.17 0.30; do
    for ratio in 1.20 1.60; do
      echo "incgraph_${lev}_${dens}_5_30_${ratio}_1"
    done
  done
done)

# roda uma configuracao e acrescenta a linha correspondente ao CSV
# $1 = nome da instancia, $2 = k, $3 = rotulo do algoritmo, $4.. = argumentos
run_cfg() {
  local name=$1 k=$2 label=$3
  shift 3

  local out="$RESULTS_DIR/${name}_k${k}_${label}_best.txt"
  local outc="$RESULTS_DIR/${name}_k${k}_${label}_complete.txt"
  rm -f "$out" "$outc"

  set +e
  "$BIN" "instance/${name}.txt" "$@" "$MAXIT" "$name" "$out" "$outc" "$TIME_LIMIT" "${EXTRA[@]}" > /dev/null 2>&1
  local status=$?
  set -e

  if [ -f "$out" ]; then
    read -r r_k r_alpha r_cost r_ttot r_tbest <<< "$(tail -n1 "$out")"
    echo "${name},${k},${label},${r_cost},${r_ttot},${r_tbest}" >> "$SUMMARY"
    printf "    %-13s custo=%-10s tempo=%ss\n" "$label" "$r_cost" "$r_ttot"
  elif [ $status -ne 0 ]; then
    echo "${name},${k},${label},,ERRO(exit=$status)," >> "$SUMMARY"
    printf "    %-13s ERRO (exit=%s)\n" "$label" "$status"
  else
    echo "${name},${k},${label},,SKIPPED," >> "$SUMMARY"
    printf "    %-13s SKIPPED (k > min. de nos incrementais)\n" "$label"
  fi
}

for name in $INSTANCES; do

  dist="$DIST_DIR/${name}.${TECHNIQUE}.dist.txt"

  if [ ! -f "instance/${name}.txt" ]; then
    echo "AVISO: instancia nao encontrada: $name" >&2
    continue
  fi
  if [ ! -f "$dist" ]; then
    echo "AVISO: distancias nao encontradas: $dist" >&2
    continue
  fi

  for k in $KVALUES; do

    echo "$name  k=$k"

    EXTRA=()
    run_cfg "$name" "$k" "grasp2" grasp2 "$k" "$ALPHA" "$LS"

    EXTRA=("$dist" 20)
    run_cfg "$name" "$k" "glgrasp-e20" glgrasp "$k" "$ALPHA" "$LS"

    EXTRA=("$dist" 100)
    run_cfg "$name" "$k" "glgrasp-e100" glgrasp "$k" "$ALPHA" "$LS"

  done

done

echo ""
echo "Concluido. Resumo em $SUMMARY"
