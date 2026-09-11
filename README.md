# GRASP (Projeto de implementação do GL-GRASP)
Futuro Código do GL-GRASP implementado em C++, baseado no código do Napoletano do GRASP

Resolve o problema **IMLCM** (Incremental Multi-Level Crossing Minimization): dado um
grafo hierárquico (`HDAG`) já desenhado em níveis, insere um conjunto de nós novos
("incrementais") minimizando o número de cruzamentos de arestas, respeitando um
parâmetro `k` que limita o quanto cada nó original pode se deslocar da sua posição
original.

## Compilação

Pré-requisitos: `g++` (C++14) e `mingw32-make` (ou `make`).

```bash
mingw32-make
```

Gera o binário `C-IGDP.exe`. Por padrão o build **não** inclui os modos `cplex` e
`lsolver`, que exigem os SDKs pagos do IBM CPLEX e do Localsolver instalados. Para
habilitá-los:

```bash
mingw32-make WITH_CPLEX=1 CPLEX_DIR=/caminho/do/cplex
mingw32-make WITH_LOCALSOLVER=1 LOCALSOLVER_DIR=/caminho/do/localsolver
```

## Rodando um teste único

```bash
./C-IGDP.exe <instancia> <algoritmo> <k> <alpha> <busca_local> <max_it> <nome> <saida_best.txt> <saida_completa.txt> <tempo_limite>
```

Exemplo:

```bash
./C-IGDP.exe instance/incgraph_13_0.06_5_30_1.20_1.txt grasp3 1 0 best 100 teste teste_best.txt teste_k1.txt 300
```

- `algoritmo`: `grasp1`, `grasp2`, `grasp3`, `tabu`, `glgrasp` (ou `cplex`/`lsolver`, se compilados)
- `k`: deslocamento máximo permitido de um nó original em relação à sua posição
  original (quanto maior, mais liberdade o algoritmo tem para reduzir cruzamentos —
  mas só é aceito se todo nível da instância tiver pelo menos `k` nós incrementais,
  senão a execução é pulada silenciosamente)
- `alpha`: fator de ganância/aleatoriedade do GRASP, em `[0,1]`. **`alpha = 0` é um
  valor-sentinela**: ativa o modo "reactive GRASP", em que um alpha aleatório
  diferente é sorteado a cada iteração, em vez de usar um valor fixo
- `busca_local`: `first`, `first-prepro`, `best` ou `no`
- `tempo_limite`: em segundos (10 a 3600)

> **Aviso sobre reprodutibilidade:** o gerador de números aleatórios (`MTRand`) é
> semeado a partir do relógio do sistema a cada execução, então rodar o mesmo
> comando duas vezes produz resultados ligeiramente diferentes.

## GL-GRASP (`glgrasp`)

Implementa a metodologia de Charytitsch e Nascimento (2026): a fase de construção
do GRASP passa a ser guiada pelas **distâncias estruturais** obtidas de *node
embeddings*, em vez do número de cruzamentos adicionais. A busca local
(*swap* + *insert*) é a mesma do `grasp3`.

Exige um argumento a mais — o arquivo de distâncias gerado pelo módulo Python
(ver [`embeddings/`](embeddings/)):

```bash
./C-IGDP.exe <instancia> glgrasp <k> <alpha> <busca_local> <max_it> <nome> <saida_best.txt> <saida_completa.txt> <tempo_limite> <arquivo_distancias> [eta_max]
```

Exemplo completo (gerar as distâncias e rodar):

```bash
cd embeddings && python build_distances.py ../instance/incgraph_6_0.06_5_30_1.20_1.txt -t hope && cd ..
./C-IGDP.exe instance/incgraph_6_0.06_5_30_1.20_1.txt glgrasp 1 0 best 100 teste gl_best.txt gl_full.txt 60 embeddings/distances/incgraph_6_0.06_5_30_1.20_1.hope.dist.txt
```

- `arquivo_distancias`: caminho do `.dist.txt`. Trocar a técnica de *embedding*
  é só apontar para outro arquivo — é assim que se compara `hope`, `spectral`,
  `node2vec`, `sdne`, `line`, `grarep`, `netsmf` e `prone`.
- `eta_max` (opcional, padrão `20`): número máximo de iterações **sem melhora**
  antes de parar ($\eta_{max}$ do artigo). O `grasp3` não tem esse critério —
  ele sempre roda até `max_it` ou até o tempo-limite.

O programa recusa executar se o arquivo de distâncias não existir ou não
corresponder à instância (níveis/ids fora da faixa), em vez de produzir
resultados silenciosamente errados.

### O que muda em relação ao `grasp3`

| Elemento | `grasp3` (Napoletano et al., 2019) | `glgrasp` |
|---|---|---|
| Critério guloso $\mathcal{G}(u)$ | mínimo de cruzamentos adicionais entre todas as posições viáveis | menor distância de *embedding* aos vizinhos já posicionados (Eq. 11) |
| Limiar da RCL | $\xi = \min + \varphi(\max-\min)$ | mesma fórmula (Eq. 12) |
| Posição de inserção | a que minimiza cruzamentos | sorteada entre: mais próxima do vizinho **mais próximo**, do **mais distante** e da **média** (Step 5) |
| Parada | `max_it` / tempo-limite | idem + $\eta_{max}$ iterações sem melhora |
| Busca local | *swap* + *insert* | idêntica |

Comparação medida (mesma máquina, `k=1`, `alpha=0`, `best`, `max_it=100`,
técnica `hope`):

| Instância | `grasp3` | `glgrasp` |
|---|---|---|
| `incgraph_6_0.06_..._1.20_1` | 1073 em 0,17 s | **1070** em 0,06 s |
| `incgraph_6_0.30_..._1.60_1` | **23845** em 8,94 s | 23844 em 1,41 s |
| `incgraph_13_0.17_..._1.60_1` | **26784** em 14,43 s | 26786 em 2,62 s |
| `incgraph_20_0.30_..._1.60_1` | **179902** em 60,38 s | 179944 em 26,73 s |

Qualidade equivalente com 2–6× menos tempo, o que é coerente com o relatado no
artigo: a construção deixa de enumerar todas as posições viáveis com o custo de
cruzamento de cada uma.

## Bateria de testes (`test_battery.sh`)

```bash
./test_battery.sh
```

Roda o algoritmo `grasp3` sobre um conjunto de instâncias representativas para
`k = 1, 2, 3` e grava um resumo consolidado em `results/summary.csv`.

**Instâncias usadas:** os arquivos em `instance/` seguem o padrão
`incgraph_<L>_<densidade>_5_30_<razão>_<seed>.txt`, onde:

| Campo | Valores no dataset | O que controla |
|---|---|---|
| `L` | `6`, `13`, `20` | número de níveis do grafo hierárquico |
| `densidade` | `0.06`, `0.17`, `0.30` | densidade de arcos (nº de nós é o mesmo entre densidades; só o nº de arcos muda) |
| `razão` | `1.20`, `1.60` | fator de crescimento = (nós totais) / (nós originais) — controla quantos nós incrementais são inseridos |
| `seed` | `1` a `10` | réplica/semente do gerador da instância (mesmos parâmetros, grafo sorteado diferente) |

A bateria varre todas as combinações de `L` × `densidade` × `razão` (18 no total),
sempre usando a `seed 1`, e testa cada uma nos três valores de `k`
(`1 2 3`, ver `KVALUES` no script) — totalizando até 54 execuções (algumas são
puladas quando `k` excede o número mínimo de nós incrementais de algum nível da
instância).

**Parâmetros fixos da bateria** (topo do `test_battery.sh`):

| Variável | Valor | Significado |
|---|---|---|
| `ALGO` | `grasp3` | algoritmo testado |
| `ALPHA` | `0` | modo reactive GRASP (alpha aleatório por iteração) |
| `LS` | `best` | tipo de busca local |
| `MAXIT` | `100` | máximo de iterações do GRASP |
| `TIME_LIMIT` | `300` | tempo-limite por execução, em segundos |
| `KVALUES` | `1 2 3` | valores de `k` testados em cada instância |

**Saída:** para cada execução, os arquivos `results/<instancia>_k<k>_best.txt` e
`results/<instancia>_k<k>_complete.txt` são gerados pelo próprio `C-IGDP.exe`, e a
última linha de cada `_best.txt` é agregada em `results/summary.csv` com as colunas
`instance,k,alpha,best_cost,total_time_s,time_to_best_s`. Execuções puladas (por
`k` inviável) aparecem marcadas como `SKIPPED`; falhas reais do binário aparecem
como `ERRO`.

Para rodar em outro conjunto de instâncias, outro algoritmo ou outro `time_limit`,
basta editar as variáveis/laços no topo do `test_battery.sh`.

Artigo Referência:
https://www.sciencedirect.com/science/article/abs/pii/S0377221718308701
