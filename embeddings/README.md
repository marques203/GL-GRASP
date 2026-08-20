# Módulo de embeddings (GL-GRASP)

Primeira metade da arquitetura proposta (FIGURA 3.1 do TG): gera embeddings
dos vértices do grafo incremental (IG), reduz para 2D via PCA e calcula as
distâncias estruturais entre vértices adjacentes, exportando-as para um
arquivo que o módulo de otimização (C++) vai consumir depois.

```
Técnica de GRL (HOPE)  ->  Redução de dimensionalidade (PCA)  ->  Distâncias estruturais
```

## Instalação

```bash
pip install -r requirements.txt
```

> `cogdl 0.6` importa `optuna`, que usa uma função do `matplotlib`
> (`matplotlib.cm.get_cmap`) removida a partir da versão 3.9. Por isso o
> `requirements.txt` fixa `matplotlib<3.9` — sem isso, `import cogdl` quebra.

## Uso

```bash
python build_distances.py <instancia.txt> <saida.txt> [--dim N] [--beta B]
```

- `instancia.txt`: arquivo no mesmo formato lido por `HDAG::read_instance`
  (ver `../instance/`).
- `--dim`: dimensão do embedding HOPE. Default segue o artigo, `k = 2n - 1`
  (`n` = número total de vértices de IG). Nota: a implementação do CogDL
  trunca internamente para `int(dim/2)` valores singulares e depois
  concatena os dois lados (esquerdo/direito), então a dimensão final do
  embedding sai como `2*floor(dim/2)`, tipicamente `2n-2` — é um
  comportamento da implementação de referência do CogDL, não um bug daqui.
- `--beta`: parâmetro de decaimento do índice de Katz (default `0.01`,
  igual ao default do CogDL).

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
- `hope_embeddings.py`: geração dos embeddings via `cogdl.models.emb.hope.HOPE`
  (decoder de Katz + SVD, Seção 4.1 do artigo).
- `build_distances.py`: script principal — orquestra embedding → PCA →
  distâncias → exportação.

## Validado em

- `incgraph_6_0.06_5_30_1.20_1.txt` (138 vértices, 160 arcos) — ~poucos segundos.
- `incgraph_20_0.30_5_30_1.60_1.txt` (676 vértices, 4120 arcos) — ~18s.
