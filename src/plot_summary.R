# =============================================================================
# 総合図: 両データセットの結果を1枚にまとめる
# =============================================================================
# Usage: Rscript src/plot_summary.R <outdir>
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(RColorBrewer)
  library(gridExtra)
})

args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# --- 1. データ読み込み ---
org <- read.csv("output/Organogenesis/evaluation/metrics.csv", stringsAsFactors = FALSE)
tes <- read.csv("output/Testis/evaluation/metrics.csv", stringsAsFactors = FALSE)

org$dataset <- "Organogenesis"
tes$dataset <- "Testis"

# --- 2. 手法の表示名 ---
label_map <- c(
  "guidedpls_none"            = "gPLS\n(none)",
  "guidedpls_stage"           = "gPLS\n(stage)",
  "guidedpls_germlayer"       = "gPLS\n(germlayer)",
  "guidedpls_stage_germlayer" = "gPLS\n(stage+\ngermlayer)",
  "guidedpls_stage_binned"    = "gPLS\n(stage\nbinned)",
  "seurat"                    = "Seurat",
  "harmony"                   = "Harmony",
  "scanorama"                 = "Scanorama"
)

# gPLS条件かどうか
is_gpls <- function(m) grepl("^guidedpls_", m)

org$label <- label_map[org$method]
tes$label <- label_map[tes$method]
org$is_gpls <- is_gpls(org$method)
tes$is_gpls <- is_gpls(tes$method)

# 順序設定
org$label <- factor(org$label, levels = label_map[
  c("guidedpls_stage", "guidedpls_germlayer", "guidedpls_stage_germlayer",
    "seurat", "harmony", "scanorama")])
tes$label <- factor(tes$label, levels = label_map[
  c("guidedpls_none", "guidedpls_stage_binned")])

# --- 3. カラー ---
gpls_cols <- colorRampPalette(c("#2166AC", "#92C5DE"))(4)
comp_cols <- brewer.pal(3, "Set2")

org$fill <- ifelse(org$is_gpls,
  gpls_cols[match(org$method,
    c("guidedpls_stage", "guidedpls_germlayer", "guidedpls_stage_germlayer"))],
  comp_cols[match(org$method, c("seurat", "harmony", "scanorama"))]
)
tes$fill <- gpls_cols[match(tes$method,
  c("guidedpls_none", "guidedpls_stage_binned"))]

# --- 4. 4指標を long format に ---
metrics_long <- function(df) {
  df %>%
    select(label, dataset, balanced_accuracy, ARI, NMI, fill) %>%
    pivot_longer(cols = c(balanced_accuracy, ARI, NMI),
      names_to = "metric", values_to = "value") %>%
    mutate(metric = factor(metric,
      levels = c("balanced_accuracy", "ARI", "NMI"),
      labels = c("Balanced\nAccuracy", "ARI", "NMI")))
}

org_long <- metrics_long(org)
tes_long <- metrics_long(tes)

# --- 5. Organogenesis パネル ---
p_org <- ggplot(org_long, aes(x = label, y = value, fill = label)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  facet_wrap(~metric, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = setNames(org$fill, org$label)) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    strip.text = element_text(size = 11, face = "bold"),
    plot.title = element_text(size = 14, face = "bold")
  ) +
  labs(x = NULL, y = "Score", title = "Organogenesis (E8.5-E10.5)")

# --- 6. Testis パネル ---
p_tes <- ggplot(tes_long, aes(x = label, y = value, fill = label)) +
  geom_col(width = 0.5, show.legend = FALSE) +
  facet_wrap(~metric, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = setNames(tes$fill, tes$label)) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    strip.text = element_text(size = 11, face = "bold"),
    plot.title = element_text(size = 14, face = "bold")
  ) +
  labs(x = NULL, y = "Score", title = "Testis (E18-P7)")

# --- 7. 結合 ---
p_combined <- arrangeGrob(p_org, p_tes,
  nrow = 2, heights = c(1.2, 1))

ggsave(file.path(outdir, "summary_both_datasets.png"),
  p_combined, width = 12, height = 10, dpi = 200)

cat("Summary figure saved.\n")
