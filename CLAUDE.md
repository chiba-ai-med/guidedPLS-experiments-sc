# CLAUDE.md

Claude Code がこのレポジトリで作業するときに最初に読むコンテキスト。

## このレポジトリの位置づけ

論文 1 本を **3 レポジトリの結果をマージ** して書くプロジェクトの一部。
- `chiba-ai-med/guidedPLS-experiments-sc` ← **このリポ** (single-cell マルチオミクス: scRNA + scATAC)
- `chiba-ai-med/guidedPLS-experiments-bulk` (bulk オミクス)
- `chiba-ai-med/ImageRegistration-experiments3` (画像レジストレーション)

各レポジトリで小さく結果をまとめ、最後にマージする方針。本リポでは guided-PLS の single-cell マルチオミクス統合性能を、既存手法 (Seurat, Harmony, Scanorama) と比較する。

## 主要ドキュメントへの索引

セッション開始時に必要なら以下を読む:

- `README.md` — workflow 全体像、DAG 画像、再現コマンド
- `DATA_CONTEXT.md` — Organogenesis / Testis データセットの取得元、ファイル構成、メタデータ仕様
- `guidedPLS-context.md` — guidedPLS パッケージ (R) 本体の API リファレンス、依存関係
- `plot/Snakefile.png` — workflow の rulegraph
- `report/Snakefile.html` — Snakemake 公式 HTML レポート

## データセット

| データセット | RNA | ATAC | ガイド条件 | 比較手法 |
|---|---|---|---|---|
| Organogenesis | GSE186068/GSE186069 (E8.5/E9.5/E10.5) | CNP0003941 fragments | none, stage, germlayer, stage_germlayer | seurat, harmony, scanorama |
| Testis | E18/P0/P3/P6 | ArchR project (`Save-ArchR-Project.rds`) | none, stage_binned | **なし** (fragment バーコード形式の問題で GAM が作れない) |

生データ (`data/`, 36GB) と中間出力 (`output/`, 14GB) は **`.gitignore` で除外**。push 不可。再生成は `snakemake -s workflow/Snakefile --cores 8 -p`。

## 実行環境

- **conda env: `snakemake`** を使うこと (`conda activate snakemake`)。Snakemake 8.10.0 + graphviz 同梱。
- `base` 環境では `conda` の dependency solver が壊滅的に遅い (1 パッケージ追加で 30 分以上固まる)。新規パッケージ追加もこの env で完結させる。
- R 側は別 conda env (`r-guidedpls`, `r_4.3` 等) が用意済み。スクリプト側で適切な env を選ぶ。

## 現在の進捗 (2026-05-25 時点)

- ✅ Organogenesis: 4 guided-PLS variants + 3 baseline 手法すべて評価完了
  - `output/Organogenesis/evaluation/metrics.csv`, `per_class_f1.csv`
  - confusion matrix 9 枚, barplot 2 枚
- ✅ Testis: 2 guided-PLS variants 評価完了, baseline は未実施 (上記の通り)
  - `output/Testis/evaluation/metrics.csv`, `per_class_f1.csv`
- ⚠️ **現状の数値**: Organogenesis で guided-PLS の最良が accuracy 16.4% (stage_germlayer)、Seurat が 47.6% で負けている。パラメータ調整・前処理見直し・guide ラベル設計の余地が大きい。

## 開発上の注意

- 論文用に厳選した図は `plot/Figures/` に手動でコピーする予定 (まだ空)。`output/**/figures/` の自動生成図と区別する。
- `workflow/dag.sh` と `workflow/report.sh` で DAG 画像と HTML レポートを再生成できる。
- Snakefile は単一エントリ。`.smk` は `include:` 経由なので個別実行は不可。
