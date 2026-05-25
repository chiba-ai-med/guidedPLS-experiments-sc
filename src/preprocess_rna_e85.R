# =============================================================================
# scRNA-Seq E8.5 前処理 (GSE186069, Matrix Market形式)
# =============================================================================
# Usage: Rscript src/preprocess_rna_e85.R \
#   <mtx.gz> <cell_annotate.csv.gz> <gene_annotate.csv.gz> \
#   <out_rds> <min_genes> <seed>
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
mtx_file   <- args[1]
cells_file <- args[2]
genes_file <- args[3]
out_rds    <- args[4]
min_genes  <- as.integer(args[5])
seed       <- as.integer(args[6])

set.seed(seed)
cat("=== scRNA-Seq E8.5 preprocessing ===\n")

# --- 1. データ読み込み ---
cat("Reading Matrix Market file ...\n")
mat <- readMM(gzfile(mtx_file))

cell_ann <- read.csv(gzfile(cells_file), stringsAsFactors = FALSE)
gene_ann <- read.csv(gzfile(genes_file), stringsAsFactors = FALSE)

# 行列の向き確認: mtxは genes × cells
cat(sprintf("Matrix: %d x %d\n", nrow(mat), ncol(mat)))
cat(sprintf("Genes: %d, Cells: %d\n", nrow(gene_ann), nrow(cell_ann)))

# 行名・列名設定
if (nrow(mat) == nrow(gene_ann) && ncol(mat) == nrow(cell_ann)) {
  rownames(mat) <- make.unique(gene_ann$gene_short_name)
  colnames(mat) <- cell_ann$sample
} else if (nrow(mat) == nrow(cell_ann) && ncol(mat) == nrow(gene_ann)) {
  mat <- t(mat)
  rownames(mat) <- make.unique(gene_ann$gene_short_name)
  colnames(mat) <- cell_ann$sample
}

# --- 2. QCフィルタ ---
# 低品質・ダブレット除外
keep_qc <- cell_ann$removed_by_low_quality_or_doublets == "No"
keep_qc[is.na(keep_qc)] <- FALSE

# celltype がNAでない
keep_ct <- !is.na(cell_ann$celltype) & cell_ann$celltype != ""

keep <- keep_qc & keep_ct
mat <- mat[, keep]
cell_ann <- cell_ann[keep, ]
cat(sprintf("After QC filter: %d cells\n", ncol(mat)))

# 最小遺伝子数フィルタ
genes_per_cell <- colSums(mat > 0)
keep_genes <- genes_per_cell >= min_genes
mat <- mat[, keep_genes]
cell_ann <- cell_ann[keep_genes, ]
cat(sprintf("After min_genes filter: %d cells\n", ncol(mat)))

# ゼロ発現遺伝子除去
gene_keep <- rowSums(mat) > 0
mat <- mat[gene_keep, ]
cat(sprintf("Non-zero genes: %d\n", nrow(mat)))

# --- 3. Seuratオブジェクト作成 ---
rna <- CreateSeuratObject(counts = mat, project = "RNA_E85")

# メタデータ追加
rownames(cell_ann) <- cell_ann$sample
rna <- AddMetaData(rna, cell_ann)
rna$Stage <- "E8.5"

# 正規化
rna <- NormalizeData(rna)

cat(sprintf("Final RNA E8.5: %d genes x %d cells\n", nrow(rna), ncol(rna)))

# --- 4. 保存 ---
saveRDS(rna, out_rds)
cat("Done.\n")
