# Módulo de embeddings (GL-GRASP)

Primeira metade da arquitetura proposta (FIGURA 3.1 do TG): gera embeddings
dos vértices do grafo incremental (IG), reduz para 2D via PCA e calcula as
distâncias estruturais entre vértices adjacentes, exportando-as para um
arquivo que o módulo de otimização (C++) vai consumir depois.

```
Técnica de GRL  ->  Redução de dimensionalidade (PCA)  ->  Distâncias estruturais
```

Suporta as 8 técnicas da tabela do TG (todas via implementação de
referência do CogDL):

| Técnica | Proximidade preservada | Determinística | Origem |
|---|---|---|---|
| `spectral` | Estrutura global de particionamento | Sim | Base |
| `hope` | Alta ordem assimétrica (Katz) | Sim | Base |
| `node2vec` | Homofilia e/ou equivalência estrutural | Não | Base |
| `sdne` | 1ª e 2ª ordem, otimização conjunta | Não | Base |
| `line` | 1ª e 2ª ordem, concatenadas | Não | Nova |
| `grarep` | Alta ordem separada por passo k | Sim | Nova |
| `netsmf` | Matriz esparsificada | Sim (mas usa amostragem aleatória internamente) | Nova |
| `prone` | Local e global (filtro passa-faixa) | Sim | Nova |

## Instalação

```bash
pip install -r requirements.txt
```

> `cogdl 0.6` importa `optuna`, que usa uma função do `matplotlib`
> (`matplotlib.cm.get_cmap`) removida a partir da versão 3.9. Por isso o
> `requirements.txt` fixa `matplotlib<3.9` — sem isso, `import cogdl` quebra.
>
> Além disso, `node2vec`/`line`/`netsmf` (via `cogdl.utils.alias_setup`) usam
> os aliases `np.int`/`np.float`/`np.bool`, removidos do NumPy a partir da
> 1.24. `embedding_techniques.py` aplica um shim de compatibilidade
> (`np.int = int` etc.) no import — não precisa fazer nada extra, mas é bom
> saber por que isso existe se aparecer de novo em outro ambiente.

## Uso

```bash
python build_distances.py <instancia.txt> [saida.txt] [--technique T] [--dim N] [--param k=v ...]
```

- `instancia.txt`: arquivo no mesmo formato lido por `HDAG::read_instance`
  (ver `../instance/`).
- `saida.txt` (opcional): se omitido, o resultado vai para
  `embeddings/distances/<instancia>.<tecnica>.dist.txt` (a pasta é criada
  automaticamente). Esses arquivos são regeneráveis e não são versionados
  (ver `.gitignore`). Ignorado quando `--technique all`.
- `--technique` / `-t`: uma das 8 da tabela acima, ou `all` para rodar
  todas de uma vez (um arquivo de saída por técnica; se uma falhar, as
  outras continuam — resumo impresso no final). Default: `hope`.
- `--dim`: dimensão do embedding. Default por técnica (ver
  `default_dimension()` em `embedding_techniques.py`): `hope` usa `2n-1` e
  `spectral` usa `n-1` (conforme o artigo, Seção 5.2); as demais usam
  `min(128, n-1)` (default de literatura do CogDL, limitado pelo tamanho
  do grafo).
- `--param k=v` (repetível): sobrescreve um parâmetro específico da
  técnica escolhida (ex.: `--param walk_num=20 --param window_size=10`
  para `node2vec`). Os nomes dos parâmetros aceitos por técnica estão
  documentados nos comentários de `_build_model()` em
  `embedding_techniques.py`.

Exemplos:

```bash
# uma tecnica so
python build_distances.py ../instance/incgraph_6_0.06_5_30_1.20_1.txt --technique spectral

# as 8 de uma vez, pra comparar depois
python build_distances.py ../instance/incgraph_6_0.06_5_30_1.20_1.txt --technique all

# ajustando um parametro especifico do node2vec
python build_distances.py ../instance/incgraph_6_0.06_5_30_1.20_1.txt --technique node2vec --param walk_num=20
```

Saída (`saida.txt`):
```
n_arcos
layer_u id_u layer_v id_v distancia
...
```
`layer_u`/`id_u` e `layer_v`/`id_v` usam a mesma numeração (nível +
id local) do arquivo de instância original — sempre `layer_v == layer_u + 1`.

## Arquivos

- `hdag_io.py`: leitura do formato de instância do C-IGDP.
- `embedding_techniques.py`: registro das 8 técnicas, adaptação da
  instância para o `Graph` do CogDL, e `compute_embeddings()` — interface
  única independente da técnica escolhida.
- `build_distances.py`: script principal — orquestra embedding → PCA →
  distâncias → exportação, com suporte a rodar uma técnica ou todas.

## Validado em

`--technique all`, 8/8 concluídas nas duas instâncias, sem falhas:

| Técnica | `incgraph_6_..._1.20_1` (138 vértices, 160 arcos) | `incgraph_20_..._1.60_1` (676 vértices, 4120 arcos) |
|---|---|---|
| `spectral` | 0.1s | 1.3s |
| `hope` | 0.1s | 1.4s |
| `node2vec` | 0.7s | 3.0s |
| `sdne` | 0.2s | 0.4s |
| `line` | 1.9s | 10.6s |
| `grarep` | 0.2s | 1.3s |
| `netsmf` | 17.9s | 21.7s |
| `prone` | 0.0s | 0.2s |

`netsmf` é a mais lenta mesmo em grafos pequenos porque usa
`multiprocessing.Pool` — no Windows, o overhead de criar processos domina
o tempo total nas instâncias menores.
