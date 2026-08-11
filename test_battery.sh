#!/bin/bash
# test_battery.sh
# Bateria de testes adaptada do TESTING.sh original.
# Roda o algoritmo grasp3 sobre um conjunto representativo de instâncias
# (uma semente por combinação de níveis/densidade/razão) para k = 1, 2, 3,
# e monta um resumo CSV em results/summary.csv
#
# Uso: ./test_battery.sh

set -e

BIN=./C-IGDP.exe
ALGO=grasp3
ALPHA=0
LS=best
MAXIT=100
TIME_LIMIT=300
KVALUES="1 2 3"

RESULTS_DIR=results
mkdir -p "$RESULTS_DIR"

SUMMARY="$RESULTS_DIR/summary.csv"
echo "instance,k,alpha,best_cost,total_time_s,time_to_best_s" > "$SUMMARY"

INSTANCES=$(for lev in 6 13 20; do
  for dens in 0.06 0.17 0.30; do
    for ratio in 1.20 1.60; do
      echo "instance/incgraph_${lev}_${dens}_5_30_${ratio}_1.txt"
    done
  done
done)

for f in $INSTANCES; do

  if [ ! -f "$f" ]; then
    echo "AVISO: instancia nao encontrada: $f" >&2
    continue
  fi

  name=$(basename "$f" .txt)

  for k in $KVALUES; do

    outbest="$RESULTS_DIR/${name}_k${k}_best.txt"
    outcomplete="$RESULTS_DIR/${name}_k${k}_complete.txt"
    rm -f "$outbest" "$outcomplete"

    echo "Rodando: $name  k=$k"

    set +e
    "$BIN" "$f" "$ALGO" "$k" "$ALPHA" "$LS" "$MAXIT" "$name" \
      "$outbest" "$outcomplete" "$TIME_LIMIT" > /dev/null
    status=$?
    set -e

    if [ -f "$outbest" ]; then
      line=$(tail -n1 "$outbest")
      read -r r_k r_alpha r_cost r_ttot r_tbest <<< "$line"
      echo "$name,$r_k,$r_alpha,$r_cost,$r_ttot,$r_tbest" >> "$SUMMARY"
    elif [ $status -ne 0 ]; then
      # binario terminou com erro (ex.: crash)
      echo "$name,$k,,,ERRO(exit=$status),ERRO" >> "$SUMMARY"
    else
      # main.cpp pula a execucao quando k > numero minimo de nos
      # incrementais em algum nivel da instancia (I.get_MIN_INCREMENTAL_NUMBER)
      echo "$name,$k,,,SKIPPED(k>min_incremental),SKIPPED" >> "$SUMMARY"
    fi

  done
done

echo ""
echo "Concluido. Resumo em $SUMMARY"
