#!/bin/bash
# experiment_all_techniques.sh
#
# Estende o experimento pareado (Analise A) as 8 tecnicas de embedding.
#
# Mesma configuracao de experiment_wilcoxon.sh: eta=100, eta_max=20 para
# todos, 30 repeticoes com as sementes 1..30, e a semente r usada por TODOS
# os metodos na repeticao r (numeros aleatorios comuns).
#
# O GRASP2 e a tecnica hope ja foram executados por experiment_wilcoxon.sh e
# nao sao refeitos aqui; a analise junta os dois arquivos.
#
# Uso: ./experiment_all_techniques.sh

set -e

BIN=./C-IGDP.exe
DIST_DIR=embeddings/distances
ALPHA=0
LS=best
MAXIT=100
TIME_LIMIT=300
ETA_MAX=20
REPS=30
KVALUES="1 2 3"
TECHNIQUES="spectral node2vec sdne line grarep netsmf prone"   # hope ja feito

RESULTS_DIR=results_all
mkdir -p "$RESULTS_DIR"

SUMMARY="$RESULTS_DIR/runs.csv"
if [ ! -f "$SUMMARY" ]; then
  echo "instance,k,algo,seed,best_cost,total_time_s" > "$SUMMARY"
fi

INSTANCES=$(for lev in 6 13 20; do
  for dens in 0.06 0.17 0.30; do
    for ratio in 1.20 1.60; do
      echo "incgraph_${lev}_${dens}_5_30_${ratio}_1"
    done
  done
done)

TMP="$RESULTS_DIR/_tmp"
mkdir -p "$TMP"

for t in $TECHNIQUES; do

  echo "########## tecnica: $t ##########"

  for name in $INSTANCES; do

    dist="$DIST_DIR/${name}.${t}.dist.txt"
    [ -f "instance/${name}.txt" ] || { echo "AVISO: sem instancia $name" >&2; continue; }
    [ -f "$dist" ]                || { echo "AVISO: sem distancias $name/$t" >&2; continue; }

    for k in $KVALUES; do

      echo "  $name k=$k"

      for seed in $(seq 1 $REPS); do

        out="$TMP/b.txt"; rm -f "$out" "$TMP/c.txt"
        set +e
        "$BIN" "instance/${name}.txt" glgrasp "$k" "$ALPHA" "$LS" "$MAXIT" "$name" \
               "$out" "$TMP/c.txt" "$TIME_LIMIT" "$dist" "$ETA_MAX" "$seed" > /dev/null 2>&1
        set -e

        if [ -f "$out" ]; then
          read -r a b cost tt rest <<< "$(tail -n1 "$out")"
          echo "${name},${k},${t},${seed},${cost},${tt}" >> "$SUMMARY"
        else
          echo "${name},${k},${t},${seed},,SKIPPED" >> "$SUMMARY"
        fi

      done

    done

  done

done

rm -rf "$TMP"

echo ""
echo "Concluido. Execucoes em $SUMMARY"
