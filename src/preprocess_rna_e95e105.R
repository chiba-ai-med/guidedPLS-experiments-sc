# =============================================================================
# scRNA-Seq E9.5-E10.5 前処理 (GSE186068, Loom形式)
# =============================================================================
# Usage: Rscript src/preprocess_rna_e95e105.R \
#   <loom.gz> <cell_annotate.csv.gz> <gene_annotate.csv.gz> \
#   <out_rds> <min_genes> <seed>
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
loom_file  <- args[1]
cells_file <- args[2]
genes_file <- args[3]
out_rds    <- args[4]
min_genes  <- as.integer(args[5])
seed       <- as.integer(args[6])

set.seed(seed)
cat("=== scRNA-Seq E9.5-E10.5 preprocessing ===\n")

# --- 1. Loom解凍 ---
loom_uncompressed <- sub("\\.gz$", "", loom_file)
if (!file.exists(loom_uncompressed)) {
  cat("Decompressing loom file ...\n")
  system2("gunzip", args = c("-k", loom_file))
}

# --- 2. Loom読み込み ---
cat("Reading Loom file ...\n")
loom <- Connect(filename = loom_uncompressed, mode = "r")
mat <- loom[["matrix"]][, ]
row_attrs <- loom[["row_attrs"]]
col_attrs <- loom[["col_attrs"]]

# Loom形式: matrix は cells × genes の場合が多い
# 確認して必要なら転置
gene_ann <- read.csv(gzfile(genes_file), stringsAsFactors = FALSE)
cell_ann <- read.csv(gzfile(cells_file), stringsAsFactors = FALSE)

cat(sprintf("Matrix dim: %d x %d\n", nrow(mat), ncol(mat)))
cat(sprintf("Genes: %d, Cells: %d\n", nrow(gene_ann), nrow(cell_ann)))

# 行列を genes × cells に整形
if (nrow(mat) == nrow(cell_ann)) {
  mat <- t(mat)  # cells × genes → genes × cells
}
mat <- as(mat, "dgCMatrix")

rownames(mat) <- gene_ann$gene_short_name
colnames(mat) <- cell_ann$sample

loom$close_all()

# --- 3. ステージフィルタ: E9.5 と E10.5 のみ ---
cell_ann$Stage <- paste0("E", cell_ann$development_stage)
keep_stage <- cell_ann$Stage %in% c("E9.5", "E10.5")
mat <- mat[, keep_stage]
cell_ann <- cell_ann[keep_stage, ]
cat(sprintf("After stage filter (E9.5, E10.5): %d cells\n", ncol(mat)))

# --- 4. QCフィルタ ---
keep_qc <- cell_ann$removed_by_low_quality_or_doublets == "No"
keep_qc[is.na(keep_qc)] <- FALSE

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

# --- 5. Seuratオブジェクト作成 ---
rna <- CreateSeuratObject(counts = mat, project = "RNA_E95_E105")
rownames(cell_ann) <- cell_ann$sample
rna <- AddMetaData(rna, cell_ann)

rna <- NormalizeData(rna)

cat(sprintf("Final RNA E9.5-E10.5: %d genes x %d cells\n", nrow(rna), ncol(rna)))

# --- 6. 保存 ---
saveRDS(rna, out_rds)
cat("Done.\n")
