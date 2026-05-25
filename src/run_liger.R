# =============================================================================
# LIGER (rliger) + kNN Label Transfer
# =============================================================================
# Usage: Rscript src/run_liger.R \
#   <rna.rds> <atac.rds> <out_labels.csv> <seed> <k>
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(rliger)
  library(RANN)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
rna_rds    <- args[1]
atac_rds   <- args[2]
out_labels <- args[3]
seed       <- as.integer(args[4])
k_nn       <- as.integer(args[5])

set.seed(seed)
cat("=== LIGER + kNN Label Transfer ===\n")

# --- 1. データ読み込み ---
rna  <- readRDS(rna_rds)
atac <- readRDS(atac_rds)
DefaultAssay(atac) <- "GAM"

# --- 2. 共通遺伝子 ---
common_genes <- intersect(VariableFeatures(rna), rownames(atac))
if (length(common_genes) < 100) {
  common_genes <- intersect(rownames(rna), rownames(atac))
  common_genes <- head(sort(common_genes), 3000)
}
cat(sprintf("Common genes: %d\n", length(common_genes)))

# --- 3. LIGER入力準備 ---
# カウント行列を共通遺伝子に絞る (genes × cells)
rna_counts  <- GetAssayData(rna, layer = "counts")[common_genes, ]
atac_counts <- GetAssayData(atac, layer = "counts")[common_genes, ]

# rliger オブジェクト作成
liger_obj <- createLiger(list(RNA = rna_counts, ATAC = atac_counts))

# --- 4. LIGER パイプライン ---
cat("Running LIGER pipeline ...\n")
liger_obj <- normalize(liger_obj)
liger_obj <- selectGenes(liger_obj)
liger_obj <- scaleNotCenter(liger_obj)
liger_obj <- optimizeALS(liger_obj, k = 20, lambda = 5)
liger_obj <- quantile_norm(liger_obj)

# --- 5. kNN Label Transfer ---
cat("kNN label transfer ...\n")

# LIGER統合後の埋め込み
H_norm <- liger_obj@H.norm

# RNA/ATAC細胞のインデックス
rna_n  <- ncol(rna_counts)
atac_n <- ncol(atac_counts)
rna_idx  <- 1:rna_n
atac_idx <- (rna_n + 1):(rna_n + atac_n)

nn <- RANN::nn2(
  data  = H_norm[rna_idx, ],
  query = H_norm[atac_idx, ],
  k = k_nn
)

rna_celltype <- rna$celltype

predicted_celltype <- apply(nn$nn.idx, 1, function(idxs) {
  neighbor_types <- rna_celltype[idxs]
  tab <- table(neighbor_types)
  names(tab)[which.max(tab)]
})

prediction_score <- apply(nn$nn.idx, 1, function(idxs) {
  neighbor_types <- rna_celltype[idxs]
  tab <- table(neighbor_types)
  max(tab) / sum(tab)
})

# --- 6. 保存 ---
result <- data.frame(
  predicted_celltype = predicted_celltype,
  prediction_score   = prediction_score,
  stringsAsFactors   = FALSE
)
write.csv(result, out_labels, row.names = FALSE)

cat(sprintf("Predicted labels saved: %d cells\n", nrow(result)))
cat("Done.\n")
