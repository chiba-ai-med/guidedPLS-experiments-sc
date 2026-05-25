# =============================================================================
# guided-PLS 実行 + kNNWarping によるラベル転送
# =============================================================================
# Usage: Rscript src/run_guidedpls.R \
#   <X1.csv> <X2.csv> <Y1.csv> <Y2.csv> \
#   <rna_metadata.csv> \
#   <out_labels.csv> <out_rdata> \
#   <r> <k>
# =============================================================================

suppressPackageStartupMessages({
  library(guidedPLS)
  library(RANN)
})

args <- commandArgs(trailingOnly = TRUE)
x1_file       <- args[1]
x2_file       <- args[2]
y1_file       <- args[3]
y2_file       <- args[4]
rna_meta_file <- args[5]
out_labels    <- args[6]
out_rdata     <- args[7]
r_dim         <- as.integer(args[8])
k_nn          <- as.integer(args[9])

cat("=== guided-PLS execution ===\n")

# --- 1. データ読み込み ---
cat("Loading data ...\n")
X1 <- as.matrix(read.csv(x1_file, header = TRUE))
X2 <- as.matrix(read.csv(x2_file, header = TRUE))
Y1 <- as.matrix(read.csv(y1_file, header = TRUE))
Y2 <- as.matrix(read.csv(y2_file, header = TRUE))

rna_meta <- read.csv(rna_meta_file, stringsAsFactors = FALSE)

cat(sprintf("X1: %d x %d (RNA)\n", nrow(X1), ncol(X1)))
cat(sprintf("X2: %d x %d (ATAC)\n", nrow(X2), ncol(X2)))
cat(sprintf("Y1: %d x %d\n", nrow(Y1), ncol(Y1)))
cat(sprintf("Y2: %d x %d\n", nrow(Y2), ncol(Y2)))

# --- 2. guided-PLS 実行 ---
cat("Running guidedPLS ...\n")
out <- guidedPLS(
  X1 = X1,
  X2 = X2,
  Y1 = Y1,
  Y2 = Y2,
  cortest = FALSE,
  verbose = TRUE
)

# --- 3. kNNWarping ---
cat("Running kNNWarping ...\n")

# スコアの標準化
scoreX1 <- apply(out$scoreX1, 2, scale)
scoreX2 <- apply(out$scoreX2, 2, scale)

# 利用可能な次元数の調整
r_use <- min(r_dim, ncol(scoreX1), ncol(scoreX2))
cat(sprintf("Using r = %d PLS dimensions\n", r_use))

# ターゲット(ATAC)各細胞に対し、ソース(RNA)潜在空間でk近傍を探索
nn <- RANN::nn2(
  data  = scoreX1[, seq(r_use)],
  query = scoreX2[, seq(r_use)],
  k = k_nn
)
idx <- nn$nn.idx

# --- 4. ラベル転送 (RNA celltype → ATAC) ---
cat("Transferring labels ...\n")

# RNA細胞のcelltype
rna_celltype <- rna_meta$celltype

# 各ATAC細胞について、k近傍RNA細胞のcelltypeの多数決
predicted_celltype <- apply(idx, 1, function(idxs) {
  neighbor_types <- rna_celltype[idxs]
  tab <- table(neighbor_types)
  names(tab)[which.max(tab)]
})

# 予測確信度 (多数決の割合)
prediction_score <- apply(idx, 1, function(idxs) {
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

save(out, nn, idx, file = out_rdata)

cat(sprintf("Predicted labels saved: %d cells\n", nrow(result)))
cat("Prediction distribution:\n")
print(table(predicted_celltype))
cat("Done.\n")
