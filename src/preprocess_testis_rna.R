# =============================================================================
# Testis scRNA-Seq 前処理
# =============================================================================
# Usage: Rscript src/preprocess_testis_rna.R \
#   <seurat.rds> <out_rds> <out_meta_csv> \
#   <subsample_n> <n_hvg> <seed>
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
input_rds    <- args[1]
out_rds      <- args[2]
out_meta_csv <- args[3]
subsample_n  <- as.integer(args[4])
n_hvg        <- as.integer(args[5])
seed         <- as.integer(args[6])

set.seed(seed)
cat("=== Testis scRNA-Seq preprocessing ===\n")

# --- 1. 読み込み + バージョン変換 ---
rna <- readRDS(input_rds)
rna <- UpdateSeuratObject(rna)
cat(sprintf("Loaded: %d genes x %d cells\n", nrow(rna), ncol(rna)))

# --- 2. ステージビニング ---
# sample: E18_1, P2_1, P7_1 → E18, P2-3, P6-7
stage_map <- c(
  "E18_1" = "E18",
  "P2_1"  = "P2-3",
  "P7_1"  = "P6-7"
)
rna@meta.data$Stage <- stage_map[rna@meta.data$sample]
cat("Stage distribution (RNA):\n")
print(table(rna$Stage))

# --- 3. サブサンプリング（ステージ層別）---
if (ncol(rna) > subsample_n) {
  cat(sprintf("Subsampling to %d cells ...\n", subsample_n))
  cells_keep <- c()
  stage_table <- table(rna$Stage)
  for (stg in names(stage_table)) {
    stg_cells <- colnames(rna)[rna$Stage == stg]
    n_take <- round(subsample_n * stage_table[stg] / sum(stage_table))
    n_take <- min(n_take, length(stg_cells))
    cells_keep <- c(cells_keep, sample(stg_cells, n_take))
  }
  rna <- subset(rna, cells = cells_keep)
}
cat(sprintf("After subsampling: %d cells\n", ncol(rna)))

# --- 4. 正規化・HVG ---
rna <- NormalizeData(rna)
rna <- FindVariableFeatures(rna, nfeatures = n_hvg)
rna <- ScaleData(rna)

# --- 5. 保存 ---
saveRDS(rna, out_rds)

meta_out <- rna@meta.data
meta_out$cell_id <- rownames(meta_out)
write.csv(meta_out, out_meta_csv, row.names = FALSE)

cat(sprintf("Final: %d genes x %d cells\n", nrow(rna), ncol(rna)))
cat("Done.\n")
