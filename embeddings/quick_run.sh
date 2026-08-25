#!/bin/bash
# quick_run.sh - teste rapido do modulo de embeddings.
# Roda so as tecnicas mais rapidas (hope, spectral, prone - nao tem random
# walk nem treino de rede, entao ficam quase instantaneas mesmo com os
# defaults "pesados" da lib) na menor instancia do repo.
#
# Uso: ./quick_run.sh

cd "$(dirname "$0")"

INSTANCE=../instance/incgraph_6_0.06_5_30_1.20_1.txt

for t in hope spectral prone; do
    python build_distances.py "$INSTANCE" --technique "$t"
done
