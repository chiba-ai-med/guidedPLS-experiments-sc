# マウス胚発生 scRNA-Seq / scATAC-Seq データ統合プロジェクト

## 目的

マウス胚発生（E8.5〜E13.5）のscRNA-SeqおよびscATAC-Seqデータを取得し、
メタデータ（発生ステージ、胚葉 GermLayer、細胞型）を細胞に紐づけた行列データを構築する。

---

## データソース

### 1. scATAC-Seq（CNP0003941）

**論文**: "Single-cell chromatin accessibility profiling of cell-state-specific gene regulatory programs during mouse organogenesis"
（Frontiers in Neuroscience, 2023）

**FTPディレクトリ**:
```
https://ftp.cngb.org/pub/CNSA/data5/CNP0003941/Single_Cell/CSE0000217/
```

**利用可能なファイル**:

| ファイル | 内容 |
|---------|------|
| `mouse.embryo.E8.5.E9.5-E10.5.metadata.upload.txt` | メタデータ（101,599細胞） |
| `Mouse_embryo_E8.5_{1-5,7-10}.fragments.tsv.gz` | E8.5 fragments（9サンプル） |
| `Mouse_embryo_E9.5_{1-10}.fragments.tsv.gz` | E9.5 fragments（10サンプル） |
| `Mouse_embryo_E10.5_{1-10}.fragments.tsv.gz` | E10.5 fragments（10サンプル） |

**メタデータの構造**（TSV形式、101,599行 + ヘッダー）:

| カラム | 説明 | 例 |
|--------|------|-----|
| cell.id | 細胞バーコード | `Mouse_embryo_E8.5_1_BC2999_N02` |
| nCount_ATAC | ATACカウント数 | 1164 |
| nFeature_ATAC | 検出ピーク数 | 841 |
| TSSEnrichment | TSS濃縮スコア | 9.164 |
| nFrags | フラグメント数 | 7435 |
| DoubletScore | ダブレットスコア | 0 |
| DoubletEnrichment | ダブレット濃縮度 | 0.1 |
| Stage | 発生ステージ | E8.5 |
| library.id | ライブラリID | E8.5_1 |
| FRiP | Fraction of Reads in Peaks | 0.1565 |
| clusters | クラスタ番号 | 13 |
| Germ.layer | 胚葉 | Mesoderm |
| celltype | 細胞型 | Mesenchymal stromal cells |
| UMAP_1 | UMAP座標1 | -1.4405 |
| UMAP_2 | UMAP座標2 | 1.9529 |

**ステージ別細胞数**:
- E8.5: 29,152細胞
- E9.5: 31,230細胞
- E10.5: 41,217細胞

**胚葉別細胞数**:
- Mesoderm: 67,260
- Neuroectoderm: 28,347
- Surface ectoderm: 4,176
- ExE embryo: 1,309
- Endoderm: 507

**細胞型**: 41種類

**fragmentsファイルの形式**（10x Genomics標準、TSV 5列）:
```
chr1  3000139  3000191  Mouse_embryo_E8.5_1_BC2999_N02  1
染色体  開始位置   終了位置   細胞バーコード（=メタデータのcell.id）  カウント
```

**重要**: fragmentsは生データであり行列データではない。
行列データ（ピーク×細胞）を得るにはピークコーリング + FeatureMatrix が必要。

---

### 2. scRNA-Seq E8.5（GSE186069）

**論文**: "Systematic reconstruction of the cellular trajectories of mammalian embryogenesis"
（Nature Genetics, 2022）

**GEOページ**: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186069

**概要**: 12個の体節分解能E8.5マウス胚から約24万細胞のsci-RNA-seq3データ

**FTPディレクトリ**:
```
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE186nnn/GSE186069/suppl/
```

**利用可能なファイル**:

| ファイル | サイズ | 内容 |
|---------|--------|------|
| `GSE186069_gene_count.mtx.gz` | 2.1 GB | 遺伝子カウント行列（Matrix Market形式） |
| `GSE186069_cell_annotate.csv.gz` | 4.3 MB | 細胞アノテーション |
| `GSE186069_gene_annotate.csv.gz` | 944 KB | 遺伝子アノテーション |

**行列データの次元**: 49,585遺伝子 × 239,533細胞（715,237,685非ゼロ要素）

**cell_annotateの構造**（CSV形式、239,533行 + ヘッダー）:

| カラム | 説明 | 例 |
|--------|------|-----|
| sample | 細胞バーコード | `P2-01A.ATTCAAGCATGTTACGCAAG` |
| UMI_count | UMI数 | 7563 |
| gene_count | 検出遺伝子数 | 3405 |
| unmatched_rate | 非マッチ率 | 0.1095 |
| doublet_score | ダブレットスコア | 0.008 |
| somite_number | 体節数 | 11 |
| embryo_sex | 胚の性別 | F |
| development_stage | 発生ステージ | 8.5 |
| celltype | 細胞型 | Forebrain/midbrain |
| removed_by_low_quality_or_doublets | QC除外フラグ | No |
| RT_group | RTグループ | E8.5_test2-1 |

**細胞型**: 30種類（NAを含む85,220細胞あり）

---

### 3. scRNA-Seq E9.5-E13.5（GSE186068）

**GEOページ**: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186068

**概要**: 61個のE9.5-E13.5マウス胚から約245万細胞のsci-RNA-seq3深層シーケンスデータ

**FTPディレクトリ**:
```
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE186nnn/GSE186068/suppl/
```

**利用可能なファイル**:

| ファイル | サイズ | 内容 |
|---------|--------|------|
| `GSE186068_gene_count.loom.gz` | 8.4 GB | 遺伝子カウント行列（Loom形式） |
| `GSE186068_cell_annotate.csv.gz` | 47 MB | 細胞アノテーション |
| `GSE186068_gene_annotate.csv.gz` | 944 KB | 遺伝子アノテーション |

**cell_annotateの構造**（CSV形式、2,452,396行 + ヘッダー）:

| カラム | 説明 | 例 |
|--------|------|-----|
| sample | 細胞バーコード | `sci3-me-001.CGAGGCAATACGCCGTTCA` |
| UMI_num | UMI数 | 361 |
| gene_number | 検出遺伝子数 | 245 |
| unmatched_rate | 非マッチ率 | 0.2031 |
| embryo_id | 胚ID | 33 |
| embryo_sex | 胚の性別 | M |
| development_stage | 発生ステージ | 12.5 |
| doublet_score | ダブレットスコア | 0.0042 |
| celltype | 細胞型 | NA（QC除外細胞の場合） |
| removed_by_low_quality_or_doublets | QC除外フラグ | Yes |

**ステージ別細胞数**:
- E9.5: 179,944
- E10.5: 439,168
- E11.5: 722,488
- E12.5: 566,350
- E13.5: 526,885

**細胞型**: 68種類

---

## データ間の関係

### 細胞IDの対応

3つのデータセットの細胞IDフォーマットは**完全に異なる**。
同一細胞のマルチオミクスデータ（paired RNA+ATAC）ではない。

| データ | 細胞ID形式 |
|--------|-----------|
| scATAC-Seq (CNP0003941) | `Mouse_embryo_E8.5_1_BC2999_N02` |
| scRNA-Seq E8.5 (GSE186069) | `P2-01A.ATTCAAGCATGTTACGCAAG` |
| scRNA-Seq E9.5-E13.5 (GSE186068) | `sci3-me-001.CGAGGCAATACGCCGTTCA` |

### 共通細胞型（scATAC-Seq と scRNA-Seq E8.5 間）: 19種類

```
Amniochorionic mesoderm A
Amniochorionic mesoderm B
Endothelium
Extraembryonic mesoderm
Extraembryonic visceral endoderm
First heart field
Forebrain/midbrain
Gut
Hindbrain
Intermediate mesoderm
Neural crest
Neuromesodermal progenitors
Paraxial mesoderm A
Placodal area
Pre-epidermal keratinocytes
Primitive erythroid cells
Somatic mesoderm
Spinal cord
Splanchnic mesoderm
```

### 共通細胞型（scATAC-Seq と scRNA-Seq E9.5-E13.5 間）: 34種類

### ステージの重複

| ステージ | scATAC-Seq | scRNA-Seq (GSE186069) | scRNA-Seq (GSE186068) |
|---------|-----------|----------------------|----------------------|
| E8.5 | 29,152 | 239,533 | - |
| E9.5 | 31,230 | - | 179,944 |
| E10.5 | 41,217 | - | 439,168 |
| E11.5 | - | - | 722,488 |
| E12.5 | - | - | 566,350 |
| E13.5 | - | - | 526,885 |

---

## 行列データの取得方法

### scRNA-Seq: そのまま行列データとして利用可能

```bash
# E8.5: Matrix Market形式（scipy.io.mmreadやReadMtxで読み込み可能）
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE186nnn/GSE186069/suppl/GSE186069_gene_count.mtx.gz"
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE186nnn/GSE186069/suppl/GSE186069_cell_annotate.csv.gz"
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE186nnn/GSE186069/suppl/GSE186069_gene_annotate.csv.gz"

# E9.5-E13.5: Loom形式（loompy/loomRで読み込み可能）
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE186nnn/GSE186068/suppl/GSE186068_gene_count.loom.gz"
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE186nnn/GSE186068/suppl/GSE186068_cell_annotate.csv.gz"
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE186nnn/GSE186068/suppl/GSE186068_gene_annotate.csv.gz"
```

### scATAC-Seq: fragmentsから行列データへの変換が必要

```bash
# メタデータ
wget "https://ftp.cngb.org/pub/CNSA/data5/CNP0003941/Single_Cell/CSE0000217/mouse.embryo.E8.5.E9.5-E10.5.metadata.upload.txt"

# fragments（必要なサンプル分をダウンロード。全29ファイル、各数百MB）
wget "https://ftp.cngb.org/pub/CNSA/data5/CNP0003941/Single_Cell/CSE0000217/Mouse_embryo_E8.5_1.fragments.tsv.gz"
# ... 他のサンプルも同様
```

**fragmentsから行列データへの変換**（SignacとMACS2が必要）:

```R
library(Signac)
library(Seurat)
library(EnsDb.Mmusculus.v79)

# 1. ピークコーリング
peaks <- CallPeaks(object = "Mouse_embryo_E8.5_1.fragments.tsv.gz")

# 2. ピーク×細胞の行列作成
peak_matrix <- FeatureMatrix(
  fragments = "Mouse_embryo_E8.5_1.fragments.tsv.gz",
  features = peaks,
  cells = cell_barcodes
)

# 3. Gene Activity Matrix作成（RNA-Seqとの統合用）
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- "UCSC"

chrom_assay <- CreateChromatinAssay(
  counts = peak_matrix,
  fragments = "Mouse_embryo_E8.5_1.fragments.tsv.gz",
  genome = "mm10",
  annotation = annotations
)

atac <- CreateSeuratObject(counts = chrom_assay, assay = "ATAC")

gene_activities <- GeneActivity(atac, gene.annotation = annotations)
```

---

## 統合解析のアプローチ

### 前提

- 同一細胞のマルチオミクスデータではない
- 細胞型ラベルベースの統合が必要

### 方法1: ラベル転送（Label Transfer）

```R
# RNA-Seqを参照として、ATAC-Seqの細胞にラベルを転送
transfer.anchors <- FindTransferAnchors(
  reference = rna_seurat,
  query = atac_seurat,
  reduction = "cca"
)

predictions <- TransferData(
  anchorset = transfer.anchors,
  refdata = rna_seurat$celltype
)
```

### 方法2: 共埋め込み（Co-embedding）

```R
# WNN (Weighted Nearest Neighbor) 解析
combined <- FindMultiModalNeighbors(...)
combined <- RunUMAP(combined, nn.name = "weighted.nn")
```

### 方法3: 細胞型レベルの疑似バルク解析

```R
# 細胞型ごとにATAC/RNAプロファイルを集約
# 遺伝子制御ネットワークの推定
```

---

## 必要なソフトウェア

### Python
```bash
pip3 install MACS2        # ピークコーリング
pip3 install scanpy       # シングルセル解析（オプション）
pip3 install loompy       # Loom形式読み込み
pip3 install scipy        # Matrix Market読み込み
```

### R
```R
install.packages(c("Seurat", "Signac", "Matrix"))

BiocManager::install(c(
  "GenomicRanges",
  "GenomeInfoDb",
  "EnsDb.Mmusculus.v79",
  "BSgenome.Mmusculus.UCSC.mm10"
))
```

### システムツール
```bash
# tabix（fragmentsのインデックス作成、高速化用）
apt-get install tabix   # Ubuntu/Debian
brew install htslib      # macOS
```

---

## ディスク容量の目安

| データ | 圧縮サイズ | 展開後の目安 |
|--------|-----------|-------------|
| scATAC-Seq fragments（全29ファイル） | 約5-10 GB | 約50-100 GB |
| scATAC-Seq メタデータ | 5 MB | 17 MB |
| scRNA-Seq E8.5 行列 | 2.1 GB | 約10 GB |
| scRNA-Seq E8.5 アノテーション | 4.3 MB | 約20 MB |
| scRNA-Seq E9.5-E13.5 行列 | 8.4 GB | 約40 GB |
| scRNA-Seq E9.5-E13.5 アノテーション | 47 MB | 約200 MB |
| **合計** | **約16-21 GB** | **約100-150 GB** |

---

## 全データ一括ダウンロードスクリプト

```bash
#!/bin/bash
set -e

OUTDIR="./data"
mkdir -p "${OUTDIR}/scATAC" "${OUTDIR}/scRNA_E85" "${OUTDIR}/scRNA_E95_E135"

# === scATAC-Seq ===
BASE_ATAC="https://ftp.cngb.org/pub/CNSA/data5/CNP0003941/Single_Cell/CSE0000217"

# メタデータ
wget -P "${OUTDIR}/scATAC" \
  "${BASE_ATAC}/mouse.embryo.E8.5.E9.5-E10.5.metadata.upload.txt"

# fragments（全29サンプル）
for stage in E8.5 E9.5 E10.5; do
  for rep in 1 2 3 4 5 6 7 8 9 10; do
    FILE="Mouse_embryo_${stage}_${rep}.fragments.tsv.gz"
    URL="${BASE_ATAC}/${FILE}"
    # ファイルが存在するか確認してからダウンロード
    if wget --spider "${URL}" 2>/dev/null; then
      wget -P "${OUTDIR}/scATAC" "${URL}"
    fi
  done
done

# === scRNA-Seq E8.5 ===
BASE_RNA85="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE186nnn/GSE186069/suppl"

wget -P "${OUTDIR}/scRNA_E85" \
  "${BASE_RNA85}/GSE186069_gene_count.mtx.gz" \
  "${BASE_RNA85}/GSE186069_cell_annotate.csv.gz" \
  "${BASE_RNA85}/GSE186069_gene_annotate.csv.gz"

# === scRNA-Seq E9.5-E13.5 ===
BASE_RNA95="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE186nnn/GSE186068/suppl"

wget -P "${OUTDIR}/scRNA_E95_E135" \
  "${BASE_RNA95}/GSE186068_gene_count.loom.gz" \
  "${BASE_RNA95}/GSE186068_cell_annotate.csv.gz" \
  "${BASE_RNA95}/GSE186068_gene_annotate.csv.gz"

echo "ダウンロード完了"
ls -lhR "${OUTDIR}"
```

---

## 結論

1. **scRNA-Seq行列データ**: GSE186069（E8.5, mtx形式）とGSE186068（E9.5-E13.5, loom形式）から直接取得可能
2. **scATAC-Seq行列データ**: CNP0003941のfragmentsファイルからSignac/MACS2で変換して取得可能
3. **メタデータ**: 発生ステージ・胚葉・細胞型は全データセットに含まれる
4. **統合**: 同一細胞のマルチオミクスではないが、19〜34の共通細胞型でラベル転送ベースの統合が可能
