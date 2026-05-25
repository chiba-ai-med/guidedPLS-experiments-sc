# =============================================================================
# 可視化: バープロット、Confusion Matrix
# =============================================================================
# Usage: Rscript src/plot_results.R \
#   <metrics.csv> <per_class_f1.csv> \
#   <guidedpls_dir> <comparison_dir> <atac_meta.csv> \
#   <conditions> <methods> <outdir>
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(RColorBrewer)
})

args <- commandArgs(trailingOnly = TRUE)
metrics_file    <- args[1]
per_class_file  <- args[2]
guidedpls_dir   <- args[3]
comparison_dir  <- args[4]
atac_meta_file  <- args[5]
conditions_str  <- args[6]
methods_str     <- args[7]
outdir          <- args[8]

conditions <- strsplit(conditions_str, ",")[[1]]
methods    <- strsplit(methods_str, ",")[[1]]

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# --- ヘルパー関数 (先に定義) ---
plot_confusion <- function(y_true, y_pred, method_name, outdir) {
  valid <- !is.na(y_true) & !is.na(y_pred)
  y_true <- y_true[valid]
  y_pred <- y_pred[valid]

  all_labels <- sort(union(unique(y_true), unique(y_pred)))
  y_true <- factor(y_true, levels = all_labels)
  y_pred <- factor(y_pred, levels = all_labels)

  conf <- as.data.frame(table(True = y_true, Predicted = y_pred))
  conf <- conf %>%
    group_by(True) %>%
    mutate(Proportion = Freq / sum(Freq)) %>%
    ungroup()

  p <- ggplot(conf, aes(x = Predicted, y = True, fill = Proportion)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "#2166AC") +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
      axis.text.y = element_text(size = 7)
    ) +
    labs(
      title = paste("Confusion Matrix:", method_name),
      x = "Predicted Cell Type",
      y = "True Cell Type"
    )

  fname <- file.path(outdir, paste0("confusion_", method_name, ".png"))
  ggsave(fname, p, width = 12, height = 10, dpi = 150)
}

# --- 1. メトリクス読み込み ---
metrics <- read.csv(metrics_file, stringsAsFactors = FALSE)
per_class <- read.csv(per_class_file, stringsAsFactors = FALSE)

# 手法の表示順序
# 比較手法がNONEまたは空の場合は除外
comp_methods <- methods[methods != "NONE" & methods != ""]
method_order <- c(
  paste0("guidedpls_", conditions),
  comp_methods
)
method_labels <- c(
  paste0("gPLS (", conditions, ")"),
  comp_methods
)
# 実在する手法のみ
keep <- method_order %in% metrics$method
method_order  <- method_order[keep]
method_labels <- method_labels[keep]

metrics$method <- factor(metrics$method,
  levels = method_order, labels = method_labels)

# --- 2. Accuracy / Balanced Accuracy バープロット ---
cat("Plotting accuracy barplot ...\n")

df_acc <- metrics %>%
  select(method, accuracy, balanced_accuracy) %>%
  pivot_longer(-method, names_to = "metric", values_to = "value")

# guided-PLS条件をハイライトするカラーマップ
n_gpls <- sum(grepl("^gPLS", method_labels))
n_comp <- length(method_labels) - n_gpls
colors <- c(
  colorRampPalette(c("#2166AC", "#92C5DE"))(n_gpls),
  brewer.pal(max(3, n_comp), "Set2")[seq_len(n_comp)]
)

p1 <- ggplot(df_acc, aes(x = method, y = value, fill = method)) +
  geom_col(width = 0.7) +
  facet_wrap(~metric, scales = "free_y",
    labeller = labeller(metric = c(
      accuracy = "Accuracy",
      balanced_accuracy = "Balanced Accuracy"
    ))) +
  scale_fill_manual(values = colors) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(x = NULL, y = "Score", title = "Label Transfer Performance")

ggsave(file.path(outdir, "barplot_accuracy.png"), p1,
  width = 10, height = 6, dpi = 150)

# --- 3. ARI / NMI バープロット ---
cat("Plotting ARI/NMI barplot ...\n")

df_ari <- metrics %>%
  select(method, ARI, NMI) %>%
  pivot_longer(-method, names_to = "metric", values_to = "value")

p2 <- ggplot(df_ari, aes(x = method, y = value, fill = method)) +
  geom_col(width = 0.7) +
  facet_wrap(~metric) +
  scale_fill_manual(values = colors) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(x = NULL, y = "Score", title = "Clustering Agreement Metrics")

ggsave(file.path(outdir, "barplot_ari_nmi.png"), p2,
  width = 10, height = 6, dpi = 150)

# --- 4. Confusion Matrix (各手法) ---
cat("Plotting confusion matrices ...\n")

atac_meta <- read.csv(atac_meta_file, stringsAsFactors = FALSE)
y_true <- atac_meta$celltype

for (cond in conditions) {
  method_id <- paste0("guidedpls_", cond)
  pred_file <- file.path(guidedpls_dir, cond, "predicted_labels.csv")
  if (!file.exists(pred_file)) next
  preds <- read.csv(pred_file, stringsAsFactors = FALSE)
  plot_confusion(y_true, preds$predicted_celltype, method_id, outdir)
}

for (m in methods) {
  pred_file <- file.path(comparison_dir, m, "predicted_labels.csv")
  if (!file.exists(pred_file)) next
  preds <- read.csv(pred_file, stringsAsFactors = FALSE)
  plot_confusion(y_true, preds$predicted_celltype, m, outdir)
}

cat("Done.\n")
