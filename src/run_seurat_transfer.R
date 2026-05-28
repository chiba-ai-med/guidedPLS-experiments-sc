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
rna_rds        <- args[1]
atac_rds       <- args[2]
out_labels     <- args[3]
seed           <- as.integer(args[4])
out_embeddings <- if (length(args) >= 5) args[5] else NA_character_

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

# --- 6. Joint CCA embedding 保存 (UMAP用) ---
if (!is.na(out_embeddings)) {
  cat("Computing joint CCA embedding for UMAP ...\n")
  # RunCCA は両側が同じ assay name のとき安定するため、ATAC GAM を
  # RNA assay として再ラップしてから CCA を走らせる
  # (Harmony スクリプトと同じパターン)。
  rna_sub <- subset(rna, features = common_genes)
  rna_sub <- NormalizeData(rna_sub, verbose = FALSE)
  rna_sub <- ScaleData(rna_sub, features = common_genes, verbose = FALSE)
  rna_sub$modality <- "RNA"

  atac_counts <- LayerData(atac, layer = "counts")[common_genes, ]
  atac_as_rna <- CreateSeuratObject(counts = atac_counts)
  atac_as_rna <- AddMetaData(atac_as_rna, atac@meta.data)
  atac_as_rna$modality <- "ATAC"
  atac_as_rna <- NormalizeData(atac_as_rna, verbose = FALSE)
  atac_as_rna <- ScaleData(atac_as_rna, features = common_genes, verbose = FALSE)

  cca <- RunCCA(rna_sub, atac_as_rna, features = common_genes,
                num.cc = n_dims, verbose = FALSE)
  emb <- Embeddings(cca, "cca")
  mod <- c(rna_sub$modality, atac_as_rna$modality)[
    match(rownames(emb), c(colnames(rna_sub), colnames(atac_as_rna)))]
  ct  <- c(rna_sub$celltype, atac_as_rna$celltype)[
    match(rownames(emb), c(colnames(rna_sub), colnames(atac_as_rna)))]

  emb_df <- data.frame(
    cell_id  = rownames(emb),
    modality = mod,
    celltype = ct,
    emb,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  write.csv(emb_df, out_embeddings, row.names = FALSE)
  cat(sprintf("Embeddings saved: %s (%d x %d)\n",
              out_embeddings, nrow(emb_df), ncol(emb)))
}
cat("Done.\n")
