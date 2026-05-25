# =============================================================================
# Seurat v5 Label Transfer (CCA-based)
# =============================================================================
# 参照上限として使用: ground truthと同フレームワークで生成されたラベルのため
# =============================================================================
# Usage: Rscript src/run_seurat_transfer.R \
#   <rna.rds> <atac.rds> <out_labels.csv> <seed>
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)
rna_rds    <- args[1]
atac_rds   <- args[2]
out_labels <- args[3]
seed       <- as.integer(args[4])

set.seed(seed)
cat("=== Seurat Label Transfer ===\n")

# --- 1. データ読み込み ---
rna  <- readRDS(rna_rds)
atac <- readRDS(atac_rds)

# Seurat v5: JoinLayers
if (inherits(rna[["RNA"]], "Assay5")) rna[["RNA"]] <- JoinLayers(rna[["RNA"]])

# ATAC: Gene Activity Matrixをデフォルトに
DefaultAssay(atac) <- "GAM"
if (inherits(atac[["GAM"]], "Assay5")) atac[["GAM"]] <- JoinLayers(atac[["GAM"]])

cat(sprintf("RNA: %d cells, ATAC: %d cells\n", ncol(rna), ncol(atac)))

# --- 2. 共通遺伝子でのアンカー探索 ---
# RNA側のHVGを使用
common_genes <- intersect(VariableFeatures(rna), rownames(atac))
if (length(common_genes) < 100) {
  common_genes <- intersect(rownames(rna), rownames(atac))
  common_genes <- head(sort(common_genes), 3000)
}
cat(sprintf("Common genes for anchors: %d\n", length(common_genes)))

# --- 3. Transfer Anchors ---
cat("Finding transfer anchors ...\n")
n_dims <- min(30, length(common_genes) - 1)
transfer.anchors <- FindTransferAnchors(
  reference  = rna,
  query      = atac,
  features   = common_genes,
  reduction  = "cca",
  dims       = 1:n_dims
)

# --- 4. Label Transfer ---
cat("Transferring labels ...\n")
# Seurat v5: weight.reduction を明示的に指定
predictions <- TransferData(
  anchorset        = transfer.anchors,
  refdata          = rna$celltype,
  dims             = 1:n_dims,
  weight.reduction = "cca"
)

# --- 5. 保存 ---
result <- data.frame(
  predicted_celltype = predictions$predicted.id,
  prediction_score   = predictions$prediction.score.max,
  stringsAsFactors   = FALSE
)
write.csv(result, out_labels, row.names = FALSE)

cat(sprintf("Predicted labels saved: %d cells\n", nrow(result)))
cat("Done.\n")
