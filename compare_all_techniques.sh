#!/bin/bash
# compare_all_techniques.sh
#
# Roda o GL-GRASP com as 8 tecnicas de embedding sobre o mesmo conjunto de
# instancias, para comparar a QUALIDADE entre elas e contra o GRASP2.
#
# Usa eta_max = maxit (criterio de parada por estagnacao desligado), de modo
# que todas as configuracoes executem o mesmo numero de iteracoes do GRASP2.
# Sem isso, a comparacao de qualidade ficaria confundida com o numero de
# iteracoes efetivamente executadas.
#
# Os resultados do GRASP2 vem da bateria anterior (results_compare/), nao
# sao refeitos aqui.
#
# Uso: ./compare_all_techniques.sh

set -e

BIN=./C-IGDP.exe
DIST_DIR=embeddings/distances
ALPHA=0
LS=best
MAXIT=100
TIME_LIMIT=300
ETA_MAX=100
KVALUES="1 2 3"
TECHNIQUES="spectral hope node2vec sdne line grarep netsmf prone"

RESULTS_DIR=results_techniques
mkdir -p "$RESULTS_DIR"

SUMMARY="$RESULTS_DIR/summary.csv"
echo "instance,k,technique,best_cost,total_time_s,time_to_best_s" > "$SUMMARY"

INSTANCES=$(for lev in 6 13 20; do
  for dens in 0.06 0.17 0.30; do
    for ratio in 1.20 1.60; do
      echo "incgraph_${lev}_${dens}_5_30_${ratio}_1"
    done
  done
done)

# ---------------------------------------------------------------------------
# 1) Garante que existem as distancias das 8 tecnicas para cada instancia
# ---------------------------------------------------------------------------
echo "### Etapa 1: geracao dos embeddings ###"
for name in $INSTANCES; do

  missing=0
  for t in $TECHNIQUES; do
    [ -f "$DIST_DIR/${name}.${t}.dist.txt" ] || missing=1
  done

  if [ $missing -eq 1 ]; then
    echo "  gerando as 8 tecnicas: $name"
    ( cd embeddings && python build_distances.py "../instance/${name}.txt" -t all > /dev/null 2>&1 )
  else
    echo "  ja completo: $name"
  fi

done

# ---------------------------------------------------------------------------
# 2) Roda o GL-GRASP com cada tecnica
# ---------------------------------------------------------------------------
echo ""
echo "### Etapa 2: execucao do GL-GRASP ###"
for name in $INSTANCES; do
  for k in $KVALUES; do

    echo "$name  k=$k"

    for t in $TECHNIQUES; do

      dist="$DIST_DIR/${name}.${t}.dist.txt"
      out="$RESULTS_DIR/${name}_k${k}_${t}_best.txt"
      outc="$RESULTS_DIR/${name}_k${k}_${t}_complete.txt"
      rm -f "$out" "$outc"

      if [ ! -f "$dist" ]; then
        echo "${name},${k},${t},,SEM_DISTANCIAS," >> "$SUMMARY"
        continue
      fi

      set +e
      "$BIN" "instance/${name}.txt" glgrasp "$k" "$ALPHA" "$LS" "$MAXIT" "$name" \
             "$out" "$outc" "$TIME_LIMIT" "$dist" "$ETA_MAX" > /dev/null 2>&1
      status=$?
      set -e

      if [ -f "$out" ]; then
        read -r r_k r_alpha r_cost r_ttot r_tbest <<< "$(tail -n1 "$out")"
        echo "${name},${k},${t},${r_cost},${r_ttot},${r_tbest}" >> "$SUMMARY"
        printf "    %-9s custo=%-11s tempo=%ss\n" "$t" "$r_cost" "$r_ttot"
      elif [ $status -ne 0 ]; then
        echo "${name},${k},${t},,ERRO(exit=$status)," >> "$SUMMARY"
        printf "    %-9s ERRO (exit=%s)\n" "$t" "$status"
      else
        echo "${name},${k},${t},,SKIPPED," >> "$SUMMARY"
        printf "    %-9s SKIPPED\n" "$t"
      fi

    done

  done
done

echo ""
echo "Concluido. Resumo em $SUMMARY"
