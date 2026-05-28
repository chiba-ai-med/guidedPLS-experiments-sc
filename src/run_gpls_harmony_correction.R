# =============================================================================
# guidedPLS の latent score に対し、modality を batch として Harmony 補正を
# かけ、補正後空間で kNN label transfer を再実行する。
# =============================================================================
# Usage: Rscript src/run_gpls_harmony_correction.R \
#   <gpls_embeddings.csv> <rna_metadata.csv> <atac_metadata.csv> \
#   <out_embeddings.csv> <out_labels.csv> <k_nn>
# =============================================================================

suppressPackageStartupMessages({
  library(harmony)
  library(RANN)
})

args <- commandArgs(trailingOnly = TRUE)
emb_file        <- args[1]
rna_meta_file   <- args[2]
atac_meta_file  <- args[3]
out_emb_file    <- args[4]
out_labels_file <- args[5]
k_nn            <- as.integer(args[6])

dir.create(dirname(out_emb_file),    recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_labels_file), recursive = TRUE, showWarnings = FALSE)

cat("=== gPLS + Harmony post-hoc correction ===\n")
d <- read.csv(emb_file, stringsAsFactors = FALSE)
meta_cols <- intersect(c("cell_id", "modality", "celltype"), colnames(d))
emb <- as.matrix(d[, setdiff(colnames(d), meta_cols), drop = FALSE])
cat(sprintf("Input embedding: %d cells x %d dims\n", nrow(emb), ncol(emb)))
cat("Modality counts:\n"); print(table(d$modality))

set.seed(42)
harm <- HarmonyMatrix(
  data_mat  = emb,
  meta_data = data.frame(modality = d$modality),
  vars_use  = "modality",
  do_pca    = FALSE,
  verbose   = FALSE
)
cat(sprintf("Corrected embedding: %d cells x %d dims\n",
            nrow(harm), ncol(harm)))

# Save corrected embeddings
emb_df <- data.frame(
  cell_id  = d$cell_id,
  modality = d$modality,
  celltype = d$celltype,
  harm,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
colnames(emb_df)[4:ncol(emb_df)] <- paste0("gPLSh_", seq_len(ncol(harm)))
write.csv(emb_df, out_emb_file, row.names = FALSE)
cat(sprintf("Embeddings saved: %s\n", out_emb_file))

# --- kNN label transfer in the corrected space ---
rna_meta  <- read.csv(rna_meta_file,  stringsAsFactors = FALSE)
atac_meta <- read.csv(atac_meta_file, stringsAsFactors = FALSE)
rna_idx   <- which(d$modality == "RNA")
atac_idx  <- which(d$modality == "ATAC")
stopifnot(length(rna_idx)  == nrow(rna_meta))
stopifnot(length(atac_idx) == nrow(atac_meta))

cat(sprintf("kNN label transfer in corrected space (k = %d) ...\n", k_nn))
nn <- RANN::nn2(
  data  = harm[rna_idx, , drop = FALSE],
  query = harm[atac_idx, , drop = FALSE],
  k = k_nn
)

rna_celltype <- rna_meta$celltype
predicted <- apply(nn$nn.idx, 1, function(idxs) {
  v <- rna_celltype[idxs]
  t <- table(v)
  names(t)[which.max(t)]
})
score <- apply(nn$nn.idx, 1, function(idxs) {
  v <- rna_celltype[idxs]
  t <- table(v)
  max(t) / sum(t)
})

result <- data.frame(
  predicted_celltype = predicted,
  prediction_score   = score,
  stringsAsFactors   = FALSE
)
write.csv(result, out_labels_file, row.names = FALSE)
cat(sprintf("Predicted labels saved: %s (%d cells)\n",
            out_labels_file, nrow(result)))
cat("Done.\n")
