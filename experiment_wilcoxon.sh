#!/bin/bash
# experiment_wilcoxon.sh
#
# Experimento pareado GRASP2 x GL-GRASP, nos moldes de Charytitsch e
# Nascimento (2026):
#
#  - 30 repeticoes por configuracao (o artigo usa 30);
#  - eta = 100 e eta_max = 20 para AMBOS os algoritmos, como o artigo
#    especifica (ele aplica esses valores tambem ao GRASP2/GRASP3);
#  - a repeticao r usa a semente r nos DOIS algoritmos (numeros aleatorios
#    comuns): na mesma repeticao os dois enfrentam os mesmos sorteios, o que
#    cancela parte do ruido na diferenca entre eles e torna o teste de
#    Wilcoxon mais sensivel.
#
# Uso: ./experiment_wilcoxon.sh

set -e

BIN=./C-IGDP.exe
DIST_DIR=embeddings/distances
TECHNIQUE=hope
ALPHA=0
LS=best
MAXIT=100
TIME_LIMIT=300
ETA_MAX=20
REPS=30
KVALUES="1 2 3"

RESULTS_DIR=results_wilcoxon
mkdir -p "$RESULTS_DIR"

SUMMARY="$RESULTS_DIR/runs.csv"
echo "instance,k,algo,seed,best_cost,total_time_s" > "$SUMMARY"

INSTANCES=$(for lev in 6 13 20; do
  for dens in 0.06 0.17 0.30; do
    for ratio in 1.20 1.60; do
      echo "incgraph_${lev}_${dens}_5_30_${ratio}_1"
    done
  done
done)

TMP="$RESULTS_DIR/_tmp"
mkdir -p "$TMP"

for name in $INSTANCES; do

  dist="$DIST_DIR/${name}.${TECHNIQUE}.dist.txt"
  [ -f "instance/${name}.txt" ] || { echo "AVISO: sem instancia $name" >&2; continue; }
  [ -f "$dist" ]                || { echo "AVISO: sem distancias $name" >&2; continue; }

  for k in $KVALUES; do

    echo "$name k=$k"

    for seed in $(seq 1 $REPS); do

      # ---- GRASP2 ----
      out="$TMP/b.txt"; rm -f "$out" "$TMP/c.txt"
      set +e
      "$BIN" "instance/${name}.txt" grasp2 "$k" "$ALPHA" "$LS" "$MAXIT" "$name" \
             "$out" "$TMP/c.txt" "$TIME_LIMIT" "$ETA_MAX" "$seed" > /dev/null 2>&1
      set -e
      if [ -f "$out" ]; then
        read -r a b cost tt rest <<< "$(tail -n1 "$out")"
        echo "${name},${k},grasp2,${seed},${cost},${tt}" >> "$SUMMARY"
      else
        echo "${name},${k},grasp2,${seed},,SKIPPED" >> "$SUMMARY"
      fi

      # ---- GL-GRASP ----
      out="$TMP/b.txt"; rm -f "$out" "$TMP/c.txt"
      set +e
      "$BIN" "instance/${name}.txt" glgrasp "$k" "$ALPHA" "$LS" "$MAXIT" "$name" \
             "$out" "$TMP/c.txt" "$TIME_LIMIT" "$dist" "$ETA_MAX" "$seed" > /dev/null 2>&1
      set -e
      if [ -f "$out" ]; then
        read -r a b cost tt rest <<< "$(tail -n1 "$out")"
        echo "${name},${k},glgrasp,${seed},${cost},${tt}" >> "$SUMMARY"
      else
        echo "${name},${k},glgrasp,${seed},,SKIPPED" >> "$SUMMARY"
      fi

    done

  done
done

rm -rf "$TMP"

echo ""
echo "Concluido. Execucoes em $SUMMARY"
