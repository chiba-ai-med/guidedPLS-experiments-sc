# guidedPLS-experiments-sc

guided-PLS による single-cell マルチオミクス (scRNA-Seq / scATAC-Seq) 統合実験。
複数データセット (Organogenesis, Testis) に対する guided-PLS の評価と比較手法 (Seurat, Harmony, Scanorama) との性能比較を行う Snakemake ワークフロー。

## Workflow

統合 Snakefile (`workflow/Snakefile`) が以下のサブワークフローを include する構成。

- **workflow/preprocess_organogenesis.smk**: Organogenesis (E8.5/E9.5/E10.5) scRNA-Seq / scATAC-Seq の前処理
- **workflow/preprocess_testis.smk**: Testis (E18/P2-3/P6-7) scRNA-Seq / scATAC-Seq の前処理
- **workflow/prepare_inputs.smk**: guided-PLS 入力データ (X1, X2, Y1, Y2) の準備
- **workflow/run_guidedpls.smk**: guided-PLS 実行
- **workflow/run_comparison.smk**: 比較手法 (Seurat, Harmony, Scanorama) の実行
- **workflow/evaluate.smk**: 評価指標 (accuracy, ARI, NMI, per-class F1) の計算と可視化

![](https://github.com/chiba-ai-med/guidedPLS-experiments-sc/blob/main/plot/Snakefile.png?raw=true)

## Requirements

- Snakemake: 8.10.0
- Graphviz (DAG 生成用)
- R / Python (各 `src/` スクリプトの依存パッケージ)

## How to reproduce this workflow

### Run the pipeline

```bash
snakemake -s workflow/Snakefile --cores 8 -p
```

データセット単位での実行:

```bash
snakemake -s workflow/Snakefile --cores 8 -p dataset --config dataset=Organogenesis
```

### Regenerate the DAG and report

```bash
bash workflow/dag.sh      # → plot/Snakefile.png
bash workflow/report.sh   # → report/Snakefile.html
```

## Outputs

- `output/{dataset}/evaluation/metrics.csv` — accuracy / ARI / NMI のまとめ
- `output/{dataset}/evaluation/per_class_f1.csv` — クラスごとの F1
- `output/{dataset}/figures/` — 棒グラフ・confusion matrix
- `plot/Figures/` — 論文用に手動でピックアップした図を格納予定
