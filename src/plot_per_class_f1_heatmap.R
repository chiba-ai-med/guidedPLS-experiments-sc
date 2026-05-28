# =============================================================================
# Per-class F1 を method × celltype のヒートマップにする
# =============================================================================
# Usage: Rscript src/plot_per_class_f1_heatmap.R \
#   <per_class_f1.csv> <conditions> <methods> <outfile>
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
per_class_file <- args[1]
conditions_str <- args[2]
methods_str    <- args[3]
outfile        <- args[4]

conditions <- strsplit(conditions_str, ",")[[1]]
methods    <- strsplit(methods_str, ",")[[1]]
comp_methods <- methods[methods != "NONE" & methods != ""]

dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)

cat("=== Per-class F1 heatmap ===\n")
df <- read.csv(per_class_file, stringsAsFactors = FALSE)

# Set method order and labels
method_order <- c(paste0("guidedpls_", conditions), comp_methods)
method_labels <- c(paste0("gPLS (", conditions, ")"), comp_methods)
keep <- method_order %in% unique(df$method)
method_order  <- method_order[keep]
method_labels <- method_labels[keep]

df$method <- factor(df$method, levels = method_order, labels = method_labels)

# Order celltypes by mean F1 across methods (highest at top)
ct_order <- df %>%
  group_by(celltype) %>%
  summarise(m = mean(f1, na.rm = TRUE)) %>%
  arrange(m) %>%
  pull(celltype)
df$celltype <- factor(df$celltype, levels = ct_order)

p <- ggplot(df, aes(x = method, y = celltype, fill = f1)) +
  geom_tile(colour = "white", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.2f", f1)),
            colour = ifelse(df$f1 > 0.5, "white", "grey20"),
            size = 2.4) +
  scale_fill_gradient2(
    low = "#f7fbff", mid = "#6baed6", high = "#08306b",
    midpoint = 0.5, limits = c(0, 1),
    name = "F1"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 7),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold")
  ) +
  labs(
    x = NULL, y = NULL,
    title = "Per-class F1 score (method × cell type)",
    subtitle = "Cell types ordered by mean F1 across methods"
  )

# Auto height based on # celltypes
h <- max(6, 0.22 * length(ct_order) + 2)
ggsave(outfile, p, width = max(7, 0.7 * length(method_order) + 4), height = h, dpi = 150)
cat("Saved:", outfile, "\n")
cat("Done.\n")
