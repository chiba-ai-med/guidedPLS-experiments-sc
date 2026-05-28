# =============================================================================
# Harmony + kNN Label Transfer
# =============================================================================
# Usage: Rscript src/run_harmony_knn.R \
#   <rna.rds> <atac.rds> <out_labels.csv> <seed> <k>
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(RANN)
})

args <- commandArgs(trailingOnly = TRUE)
rna_rds        <- args[1]
atac_rds       <- args[2]
out_labels     <- args[3]
seed           <- as.integer(args[4])
k_nn           <- as.integer(args[5])
out_embeddings <- if (length(args) >= 6) args[6] else NA_character_

set.seed(seed)
cat("=== Harmony + kNN Label Transfer ===\n")

# --- 1. データ読み込み ---
rna  <- readRDS(rna_rds)
atac <- readRDS(atac_rds)
if (inherits(rna[["RNA"]], "Assay5")) rna[["RNA"]] <- JoinLayers(rna[["RNA"]])
DefaultAssay(atac) <- "GAM"
if (inherits(atac[["GAM"]], "Assay5")) atac[["GAM"]] <- JoinLayers(atac[["GAM"]])

# --- 2. データ統合 ---
# 共通遺伝子
common_genes <- intersect(VariableFeatures(rna), rownames(atac))
if (length(common_genes) < 100) {
  common_genes <- intersect(rownames(rna), rownames(atac))
  common_genes <- head(sort(common_genes), 3000)
}
cat(sprintf("Common genes: %d\n", length(common_genes)))

# マージ
rna$modality <- "RNA"
atac$modality <- "ATAC"

# ATAC側はGAMのRNA assayとして扱う
atac_rna_assay <- CreateAssayObject(
  counts = LayerData(atac, layer = "counts")[common_genes, ]
)
atac_for_merge <- CreateSeuratObject(counts = atac_rna_assay)
atac_for_merge <- AddMetaData(atac_for_merge, atac@meta.data)
atac_for_merge$modality <- "ATAC"
atac_for_merge <- NormalizeData(atac_for_merge)

# RNA側を共通遺伝子に絞る
rna_sub <- subset(rna, features = common_genes)
rna_sub$modality <- "RNA"

combined <- merge(rna_sub, atac_for_merge)
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined)
combined <- ScaleData(combined, features = common_genes)
combined <- RunPCA(combined, features = common_genes, npcs = 30)

# --- 3. Harmony バッチ補正 ---
cat("Running Harmony ...\n")
combined <- RunHarmony(combined, group.by.vars = "modality")

# --- 4. kNN Label Transfer ---
cat("kNN label transfer ...\n")
harmony_embeddings <- Embeddings(combined, reduction = "harmony")

rna_cells  <- which(combined$modality == "RNA")
atac_cells <- which(combined$modality == "ATAC")

nn <- RANN::nn2(
  data  = harmony_embeddings[rna_cells, ],
  query = harmony_embeddings[atac_cells, ],
  k = k_nn
)

rna_celltype <- combined$celltype[rna_cells]

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

# --- 5. 保存 ---
result <- data.frame(
  predicted_celltype = predicted_celltype,
  prediction_score   = prediction_score,
  stringsAsFactors   = FALSE
)
write.csv(result, out_labels, row.names = FALSE)

cat(sprintf("Predicted labels saved: %d cells\n", nrow(result)))

# --- 6. Embedding 保存 (UMAP用) ---
if (!is.na(out_embeddings)) {
  cat("Saving harmony embeddings ...\n")
  emb_df <- data.frame(
    cell_id  = rownames(harmony_embeddings),
    modality = combined$modality,
    celltype = combined$celltype,
    harmony_embeddings,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  write.csv(emb_df, out_embeddings, row.names = FALSE)
  cat(sprintf("Embeddings saved: %s (%d x %d)\n",
              out_embeddings, nrow(emb_df), ncol(harmony_embeddings)))
}
cat("Done.\n")
