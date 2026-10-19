# Testes das mudanças de alinhamento do GL-GRASP

Validam as mudanças feitas para alinhar nossa implementação ao código de
referência de Charytitsch e Nascimento (2026), `Cb_*.py`.

## Como rodar

```bash
gmake test                              # testes de unidade (C++)
python tests/make_fixtures.py           # gera as instâncias sintéticas
python tests/validate_changes.py        # validação de integração completa
python tests/validate_changes.py --quick  # subconjunto, ~2 min
```

`gmake test` compila `tests/test_glgrasp_units.exe` ligando o header de
produção `GLGRASP.h`, não uma cópia: se a assinatura de `gl_func_ale` ou de
`gl_rcl_pick_uniform` mudar, o teste deixa de compilar.

## O que é coberto

| Mudança | Teste |
|---|---|
| 1.1 sorteio uniforme na RCL | unidade: dreno exato, uniformidade (±2% em 1,6M amostras), determinismo, e **contraste** mostrando que o esquema "nível primeiro" da classe base erra por 300% |
| 1.2 regra de posição do Step 5 | unidade: suporte `{min, max, ponto médio}`, uniformidade entre as três, bordas (`a == b`, simetria, overflow). Integração: T3 confirma que `figure` e `text` dão resultados diferentes |
| 1.3 semente da construção por repetição | T2: com $\varphi$ fixo, `max_it = 1` e sem busca local, sementes diferentes têm de dar custos diferentes — com o código antigo (`setseed(i % 100)`) davam todos o mesmo |
| 1.4 limite inferior da Claim 2 | T5: nenhuma solução inviável em toda a bateria |
| 2.1 $\mathcal{G}$ indefinido na Eq. (12) | T3 (a flag muda o resultado) e T4 (a flag age **só** onde deve) |

## Por que existem fixtures sintéticas

Todas as 18 instâncias reais têm pelo menos um vértice incremental sem vizinho
original, logo a flag do item 2.1 está **sempre** ativa nelas e não dá para
testar onde ela *não* deve agir. `tests/make_fixtures.py` gera:

- `synth_anchored` — todo incremental tem vizinho original, logo
  $\mathcal{G}(u)$ nunca fica indefinido e `inf` tem de ser **idêntico** a
  `maxdis`;
- `synth_unanchored` — tem incrementais sem âncora, um vértice de grau zero e
  um arco com distância *outlier* (50,0 contra ~0,95 dos demais), que é o que
  isola o mecanismo da Eq. (12): com `maxdis` esse valor passa a ser
  $\max\mathcal{G}$ e infla $\xi$.

O segundo é **controle positivo**: sem ele, o teste do primeiro passaria por
vacuidade. Na prática a separação é total — 0/20 sementes divergem em
`synth_anchored`, 20/20 em `synth_unanchored`.

## Determinismo das técnicas de *embedding*

O GL-GRASP usa **um embedding fixo durante toda a execução**. Para que o
experimento seja reproduzível ponta a ponta, a geração é semeada por padrão
(`DEFAULT_SEED = 1` em `build_distances.py`): o mesmo comando sempre produz o
mesmo `.dist.txt`.

Quais técnicas precisavam disso foi **medido**, não suposto — duas execuções no
mesmo grafo, comparando o arquivo byte a byte:

| Técnica | Sem semente | Com semente fixa reproduz? |
|---|---|---|
| `spectral`, `hope`, `grarep`, `netsmf`, `prone` | determinísticas | sim |
| `sdne`, `line` | **estocásticas** | sim |
| `node2vec` | **estocástica** | só com `worker=1` |

Duas correções do que eu havia previsto antes de medir: `netsmf` é
determinística, e `node2vec` **não** fica reproduzível só com a semente — o
CogDL passa `worker=10` ao `Word2Vec` do gensim, e o treino com várias threads
varia independentemente da semente. Por isso o default do `worker` foi mudado
para 1 em `embedding_techniques.py`. Custa ~2× o tempo de geração do node2vec
(≈38 s contra ≈18 s por instância) e não altera o algoritmo, só o paralelismo.

Para verificar que continua valendo:

```bash
cd embeddings
for t in spectral hope node2vec sdne line grarep netsmf prone; do
  python build_distances.py ../instance/incgraph_6_0.06_5_30_1.20_1.txt -t $t >/dev/null 2>&1
  md5sum distances/incgraph_6_0.06_5_30_1.20_1.$t.dist.txt
done
```

Rodar duas vezes tem de dar os mesmos *hashes*.

## Regressão das outras heurísticas

As mudanças são locais ao `GLGRASP` (mais `EmbeddingDistances`, o `case 6` do
`main.cpp` e o validador do `glgrasp`); `GRASP::extract_from_RCL` da classe
base ficou intacto de propósito, para que `grasp1`/`grasp2`/`grasp3` continuem
idênticos ao baseline de Napoletano et al. Verificado em 40 execuções
amostradas de `results_wilcoxon/runs.csv`: todas bit-idênticas ao armazenado.
