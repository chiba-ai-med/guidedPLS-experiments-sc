# CLAUDE.md

Claude Code がこのレポジトリで作業するときに最初に読むコンテキスト。

## このレポジトリの位置づけ

論文 1 本を **4 レポジトリの結果をマージ** して書くプロジェクトの一部。
- `chiba-ai-med/guidedPLS-experiments-sc` ← **このリポ** (single-cell マルチオミクス: scRNA + scATAC, Fig.3)
- `chiba-ai-med/guidedPLS-experiments-sim` (合成データベンチマーク, Fig.2)
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

| データセット | 用途 | RNA | ATAC | ガイド条件 | 比較手法 |
|---|---|---|---|---|---|
| **PBMC** | **Main** | 10x 10k Multiome (h5) | 同上 | germlayer (broad_lineage を流用: T/B/NK/Myeloid/Other) | seurat, harmony, scanorama |
| Organogenesis | Supp | GSE186068/GSE186069 (E8.5/E9.5/E10.5) | CNP0003941 fragments | stage, germlayer, stage_germlayer | seurat, harmony, scanorama |

PBMC は同一細胞ペアなので preprocess で cell barcode に `_R` / `_A` の suffix を付けて分離保存している (Harmony merge 等が duplicate cell name で死ぬのを防ぐため)。

旧データセット (drop 済):
- **Testis**: ArchR scATAC project は barcode 形式が Signac の GAM 構築と非互換で、比較手法が走らない。説明コストが高すぎるので scope 外。生データ・出力は `data/Testis/`, `output/Testis/` にそのまま残っている (config から外しただけ)。

生データ (`data/`) と中間出力 (`output/`) は **`.gitignore` で除外**。push 不可。再生成は `snakemake -s workflow/Snakefile --cores 8 -p`。

## 実行環境

- **conda env: `snakemake`** を使うこと (`conda activate snakemake`)。Snakemake 8.10.0 + graphviz 同梱。
- `base` 環境では `conda` の dependency solver が壊滅的に遅い (1 パッケージ追加で 30 分以上固まる)。新規パッケージ追加もこの env で完結させる。
- R 側は **`r-guidedpls`** env を使う (Rscript + python3 + scanorama すべて同居)。snakemake 自体は `snakemake` env、R/Python スクリプトは `r-guidedpls` env なので、フル実行は次のように PATH を通してから:
  ```bash
  export PATH=/home/koki/anaconda3/envs/r-guidedpls/bin:$PATH
  conda activate snakemake
  snakemake -s workflow/Snakefile --cores 4 -p
  ```
  Snakefile の各ルールは `Rscript` / `python3` をそのまま呼ぶので、これを忘れると `command not found` で死ぬ。

## 現在の進捗 (2026-05-28 時点)

### 戦略 pivot
- Organogenesis (41 celltype) で gPLS は accuracy 16.4% (Seurat 47.6%) と大きく負け、UMAP でモダリティ分離が露骨。原因は (1) PLS の潜在次元が Y rank で頭打ち (9 次元 vs ベースライン 30-100 次元)、(2) RNA HVG と ATAC raw peak で feature 空間が乖離、(3) gPLS の目的関数にモダリティ結合項がない、の 3 つ。
- 改善案として interaction-Y や post-hoc Harmony を試したが、Harmony 補正は accuracy をむしろ悪化させた (-0.022)。
- **方針転換**: **細胞型数の小さい PBMC scMultiome (10x 10k, ~10 celltype) で再評価**。同じ pipeline + 同じ 3 baseline で勝負して、「celltype 数が少ない設定では gPLS が動く」ことを示す。Organogenesis は Supp に降格 (「大規模 dataset では限界」の例)。

### Organogenesis (Supp 行き)
- ✅ 4 guide condition (`stage`, `germlayer`, `stage_germlayer`; `none` は SVD 退化で削除済) + 3 baseline 評価完了
- ✅ 視覚化拡張済: gPLS UMAP、手法間 UMAP 比較、per-class F1 heatmap、alluvial flow
- 数値: gPLS_stage_germlayer 16.4% vs Seurat 47.6%

### PBMC (Main 候補)
- ✅ 10x 10k Multiome データダウンロード `data/PBMC/` (~2.7GB)
- ✅ `src/preprocess_pbmc.R` + `workflow/preprocess_pbmc.smk` 作成
- 🔄 preprocess 実行中 / 完了次第 gPLS + baseline 走らせる

### Testis (drop)
- 比較手法非互換のため scope 外。生データ/出力は `data/Testis/`, `output/Testis/` に残置。

## 開発上の注意

- 論文用に厳選した図は `plot/Figures/` に手動でコピーする予定 (まだ空)。`output/**/figures/` の自動生成図と区別する。
- `workflow/dag.sh` と `workflow/report.sh` で DAG 画像と HTML レポートを再生成できる。
- Snakefile は単一エントリ。`.smk` は `include:` 経由なので個別実行は不可。
