# =============================================================================
# guidedpls.RData の scoreX1 / scoreX2 を embeddings.csv にまとめる
# =============================================================================
# Usage: Rscript src/save_gpls_embeddings.R \
#   <guidedpls.RData> <rna_metadata.csv> <atac_metadata.csv> <out_embeddings.csv>
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
rdata_file     <- args[1]
rna_meta_file  <- args[2]
atac_meta_file <- args[3]
out_file       <- args[4]

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

load(rdata_file)
s1 <- out$scoreX1
s2 <- out$scoreX2

rna_meta  <- read.csv(rna_meta_file,  stringsAsFactors = FALSE)
atac_meta <- read.csv(atac_meta_file, stringsAsFactors = FALSE)

n1 <- nrow(s1); n2 <- nrow(s2)
stopifnot(n1 == nrow(rna_meta), n2 == nrow(atac_meta))

emb <- rbind(s1, s2)
colnames(emb) <- paste0("gPLS_", seq_len(ncol(emb)))

df <- data.frame(
  cell_id  = c(rna_meta$cell_id, atac_meta$cell_id),
  modality = c(rep("RNA", n1), rep("ATAC", n2)),
  celltype = c(rna_meta$celltype, atac_meta$celltype),
  emb,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
write.csv(df, out_file, row.names = FALSE)
cat(sprintf("Saved: %s (%d x %d)\n", out_file, nrow(df), ncol(emb)))
