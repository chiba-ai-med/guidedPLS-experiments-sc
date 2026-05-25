# GuidedPLS Context File

GuidedPLS（Guided Partial Least Squares）を新しいオミックスデータ・計算環境で利用するための完全なリファレンス。

## 概要

GuidedPLSは、**異なるモダリティの空間オミックスデータ間でクロスモーダルアライメント**を行う手法。
解剖学的アノテーション（anatomy）をガイドとして、ソース空間の特徴量をターゲット空間に輸送（warp）する。

**典型的なユースケース**:
- 脂質MSI → 空間トランスクリプトミクス（遺伝子発現）
- MALDI → Xenium
- 同じ組織の異なるスライス・異なる測定モダリティ間のアライメント

## R パッケージ

```r
# guidedPLS のインストール (GitHub経由)
# BiocManager経由でインストールする場合:
BiocManager::install("rikenbit/guidedPLS")

# または devtools で:
devtools::install_github("rikenbit/guidedPLS")
```

### 依存パッケージ

```r
library("guidedPLS")   # メイン手法
library("RANN")         # kNNWarping用（nn2関数）
library("ggplot2")      # 可視化
library("dplyr")        # データ操作
library("viridis")      # カラーパレット
library("fields")       # 空間データのsmoothing
library("tagcloud")     # smoothPalette
library("RColorBrewer") # 追加カラースキーム
```

### conda 環境構築例

```yaml
name: r-guidedpls
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - r-base=4.3
  - r-essentials
  - r-devtools
  - r-tidyverse
  - r-dplyr
  - r-fields
  - r-rcolorbrewer
  - r-matrix
  - r-irlba
  - r-biocmanager
  - compilers
  - git
```

conda環境構築後、R内で `BiocManager::install("rikenbit/guidedPLS")` と `install.packages("RANN")` を実行。

## 入力データ形式

GuidedPLSに必要な4つの入力行列:

| 引数 | 内容 | 形式 |
|------|------|------|
| `X1` (source expression) | ソース発現行列 | CSV, header=TRUE, n_source_spots × n_source_features |
| `X2` (target expression) | ターゲット発現行列 | CSV, header=TRUE, n_target_spots × n_target_features |
| `Y1` (source anatomy) | ソース解剖学ラベル | CSV, header=TRUE, n_source_spots × n_regions (one-hot) |
| `Y2` (target anatomy) | ターゲット解剖学ラベル | CSV, header=TRUE, n_target_spots × n_regions (one-hot) |

### 重要な要件

1. **ソースとターゲットの特徴量は異なっていてよい**（クロスモーダル）
2. **解剖学ラベルはソースとターゲットで共通のカテゴリ**を使う（one-hotエンコーディングの列が同じ）
3. スポット数（行数）はソースとターゲットで異なっていてよい
4. 座標 (x, y) はGuidedPLS自体には不要（後の可視化で使用）

### データ前処理の流れ

```
元のCSV（X, Y, region, feature1, feature2, ...）
    ↓
1. region列がNAまたは空文字の行を除去
2. X, Y座標で行をソート
3. 座標(x, y)を分離して保存
4. region列を分離
5. regionをone-hotエンコーディング（ソースとターゲットで同じカラム集合）
6. 発現行列をCSV保存（header=TRUE）
7. one-hot anatomy行列をCSV保存（header=TRUE）
```

**one-hotエンコーディングのポイント**: ソースとターゲットで全ラベルの和集合をカラムとして使う:

```python
import pandas as pd

all_labels = sorted(pd.Index(s_source.unique()).union(s_target.unique()))
ohe_source = pd.get_dummies(s_source).reindex(columns=all_labels, fill_value=0).astype(int)
ohe_target = pd.get_dummies(s_target).reindex(columns=all_labels, fill_value=0).astype(int)
```

## GuidedPLS 実行コード

```r
source_all_exp <- as.matrix(read.csv("data/source/all_exp.csv", header=TRUE))
target_all_exp <- as.matrix(read.csv("data/target/all_exp.csv", header=TRUE))
source_anatomy <- as.matrix(read.csv("data/source/anatomy.csv", header=TRUE))
target_anatomy <- as.matrix(read.csv("data/target/anatomy.csv", header=TRUE))

# GuidedPLS実行
out <- guidedPLS(
    X1 = source_all_exp,
    X2 = target_all_exp,
    Y1 = source_anatomy,
    Y2 = target_anatomy,
    cortest = TRUE,
    verbose = TRUE
)

# kNNWarping: 潜在空間でのk近傍によるソース→ターゲット輸送
out2 <- kNNWarping(out, source_all_exp, source_anatomy,
    r = 18, k = 11)

# 結果の保存
write.table(out2$warped_exp, "output/guidedpls/warped.txt",
    col.names = TRUE, row.names = FALSE, sep = ",")
save(out, out2, file = "output/guidedpls/guidedpls.RData")
```

### kNNWarping関数の実装

```r
kNNWarping <- function(out, source_all_exp, source_anatomy, r = 5, k = 5) {
    scoreX1 <- apply(out$scoreX1, 2, scale)   # ソースのPLSスコアを標準化
    scoreX2 <- apply(out$scoreX2, 2, scale)   # ターゲットのPLSスコアを標準化
    r <- min(r, ncol(scoreX1), ncol(scoreX2)) # 使用次元数の調整

    # ターゲット各スポットに対し、ソース潜在空間でk近傍を探索
    nn <- RANN::nn2(
        data  = scoreX1[, seq(r)],
        query = scoreX2[, seq(r)],
        k = k
    )
    idx <- nn$nn.idx[, 1:k]

    # 近傍ソーススポットの発現量を平均して輸送後の発現量とする
    warped_exp <- t(apply(idx, 1, function(idxs) {
        colMeans(source_all_exp[idxs, , drop = FALSE])
    }))

    # 近傍ソーススポットの解剖学ラベルの多数決
    warped_anatomy <- apply(idx, 1, function(idxs) {
        tmp <- colMeans(source_anatomy[idxs, , drop = FALSE])
        names(tmp)[which.max(tmp)]
    })

    list(warped_exp = warped_exp, warped_anatomy = warped_anatomy)
}
```

### パラメータガイド

| パラメータ | 説明 | 推奨値 | 備考 |
|-----------|------|--------|------|
| `r` | 使用するPLS次元数 | 18 | データの複雑さに応じて調整。小さいデータなら5-10 |
| `k` | k近傍の数 | 11 | 大きくすると平滑化が強まる。5-20程度 |
| `cortest` | 相関検定 | TRUE | |
| `verbose` | 詳細出力 | TRUE | デバッグ時に有用 |

## 評価方法

輸送後のソース発現と、ターゲット側の既知マーカー遺伝子間のPearson相関係数(CC)で評価:

```r
# マーカーペアの定義（ドメイン知識に基づく）
source_markers <- c("Feature_A", "Feature_B")  # ソース側の対応する特徴量
target_markers <- c("Gene_X", "Gene_Y")        # ターゲット側のマーカー遺伝子

# 相関行列の計算
cor_combination <- function(warped_exp, target_all_exp,
    source_markers, target_markers) {
    cor_mat <- matrix(NA,
        nrow = length(source_markers),
        ncol = length(target_markers),
        dimnames = list(source_markers, target_markers))
    for (s in source_markers) {
        for (t in target_markers) {
            cor_mat[s, t] <- cor(
                warped_exp[, s],
                target_all_exp[, t],
                method = "pearson")
        }
    }
    unlist(cor_mat)
}

cc <- cor_combination(out2$warped_exp, target_all_exp,
    source_markers, target_markers)
```

**マーカー選定の考え方**: ソースとターゲットで生物学的に対応関係が既知の特徴量ペアを選ぶ。
例: 脂質HexCer ↔ ミエリン関連遺伝子Mog（どちらもミエリンに富む領域で高発現）

## 可視化

### 輸送後の空間発現プロット

```r
library("tagcloud")
library("viridis")

# カラーマッピング
.mycolor <- function(z) {
    smoothPalette(z,
        palfunc = colorRampPalette(viridis(100), alpha = TRUE))
}

# 空間プロット（y軸反転）
.plot_tissue_section <- function(x, y, z, cex = 1) {
    plot(x, -y, col = .mycolor(z), pch = 16, cex = cex,
        xaxt = "n", yaxt = "n", xlab = "", ylab = "", axes = FALSE)
}

# warped発現の可視化（ターゲット座標上にプロット）
warped_exp <- read.csv("output/guidedpls/warped.txt", header = TRUE)
target_x <- unlist(read.csv("data/target/x.csv", header = FALSE))
target_y <- unlist(read.csv("data/target/y.csv", header = FALSE))

for (i in seq_len(ncol(warped_exp))) {
    filename <- paste0("plot/guidedpls/", colnames(warped_exp)[i], ".png")
    png(filename, width = 1200, height = 1200, bg = "transparent")
    .plot_tissue_section(target_x, target_y, warped_exp[, i], cex = 3.5)
    dev.off()
}
```

### PLS潜在空間のペアプロット

```r
load("output/guidedpls/guidedpls.RData")  # out, out2 が読み込まれる

scoreX1 <- apply(out$scoreX1, 2, scale)
scoreX2 <- apply(out$scoreX2, 2, scale)
ndim <- min(10, ncol(scoreX1), ncol(scoreX2))

mat_all <- rbind(scoreX1[, 1:ndim], scoreX2[, 1:ndim])

# バッチ別カラーリング（赤=ソース、青=ターゲット）
col_batch <- c(rep("red", nrow(scoreX1)), rep("blue", nrow(scoreX2)))
png("plot/guidedpls/pairplot_batch.png", width = 2000, height = 2000)
pairs(mat_all, labels = paste("Score", 1:ndim),
    cex.labels = 3, pch = 16, cex = 2, col = col_batch)
dev.off()
```

## 出力ファイル

| ファイル | 内容 |
|---------|------|
| `warped.txt` | 輸送後の発現行列 (n_target_spots × n_source_features, CSV) |
| `guidedpls.RData` | guidedPLSオブジェクト（`out`, `out2`）。潜在空間スコア等を含む |
| `cc.csv` | マーカーペアごとの相関係数 |

### warped.txt の意味

- 行数 = ターゲットのスポット数
- 列数 = ソースの特徴量数
- 各行は「そのターゲットスポットに対応するソース発現量の推定値」
- ターゲット座標上にプロットすることで、ソースの特徴量をターゲット空間で可視化できる

## 新しいデータへの適用チェックリスト

1. **データ準備**
   - [ ] ソースとターゲットの発現行列CSV（header=TRUE, スポット×特徴量）
   - [ ] 共通の解剖学ラベルによるone-hot行列CSV（header=TRUE, スポット×リージョン）
   - [ ] 座標CSV（header=FALSE, スポット数の1列ベクトル）
   - [ ] 評価用マーカーペアの定義（ソース特徴量名とターゲット特徴量名の対応）

2. **環境構築**
   - [ ] R >= 4.3
   - [ ] `BiocManager::install("rikenbit/guidedPLS")`
   - [ ] `install.packages(c("RANN", "ggplot2", "dplyr", "viridis", "fields", "tagcloud", "RColorBrewer"))`

3. **パラメータ調整**
   - [ ] `r`（PLS次元数）: データ規模・複雑さに応じて 5-30
   - [ ] `k`（近傍数）: 5-20、大きいほど平滑化

4. **注意点**
   - one-hotのカラム集合はソースとターゲットで完全一致させる（片方にしかないラベルも0列として含める）
   - 発現値にNA/Inf/負値がないことを確認
   - ゼロ分散の列（全スポットで同一値）があるとPLSが不安定になる可能性あり
   - 大規模データ（数万スポット以上）でもメモリ1GB程度で動作可能

## 完全なパイプライン例（Snakemake不使用）

```bash
# 1. 環境構築
conda create -n guidedpls r-base=4.3 r-devtools r-dplyr r-fields r-rcolorbrewer r-biocmanager
conda activate guidedpls
R -e 'BiocManager::install("rikenbit/guidedPLS"); install.packages(c("RANN","ggplot2","viridis","tagcloud","fields"))'

# 2. 前処理（Pythonで実行）
python preprocess.py

# 3. GuidedPLS実行
Rscript guidedpls.R \
    data/source/all_exp.csv \
    data/target/all_exp.csv \
    data/source/anatomy.csv \
    data/target/anatomy.csv \
    output/guidedpls/warped.txt \
    output/guidedpls/guidedpls.RData

# 4. 評価
Rscript evaluate.R \
    output/guidedpls/warped.txt \
    data/target/all_exp.csv \
    output/guidedpls/cc.csv
```

## 過去の実験からの知見

- 251208データ（47K×40K spots）: GuidedPLSが唯一有効な手法だった（CC > 0）。OT手法（qGW, FRLC, LR-GW）はすべてCC ≈ 0
- kidneyデータ: IR手法（ANTsPy/SimpleITK）も一定の性能を示した
- kNNWarpingの`r=18, k=11`は251208データで経験的に決定されたパラメータ
- 解剖学ラベルの粒度（粗い vs 細かい）が結果に大きく影響する可能性あり
