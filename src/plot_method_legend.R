# =============================================================================
# 手法名 → 色 (Dark2) の対応を、横方向に展開された凡例だけを描画した
# 透明背景 PNG として保存する
# =============================================================================
# Usage: Rscript src/plot_method_legend.R \
#   <conditions> <methods> <dataset_name> <outfile>
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(RColorBrewer)
  library(cowplot)
})

args <- commandArgs(trailingOnly = TRUE)
conditions_str <- args[1]
methods_str    <- args[2]
dataset_name   <- args[3]
outfile        <- args[4]

dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)

conditions   <- strsplit(conditions_str, ",")[[1]]
methods      <- strsplit(methods_str, ",")[[1]]
comp_methods <- methods[methods != "NONE" & methods != ""]

gpls_labels <- if (length(conditions) <= 1) "gPLS" else
  paste0("gPLS (", conditions, ")")
method_labels <- c(gpls_labels, comp_methods)

n_m <- length(method_labels)
colors <- if (n_m <= 8) brewer.pal(max(3, n_m), "Dark2")[seq_len(n_m)] else
  colorRampPalette(brewer.pal(8, "Dark2"))(n_m)
names(colors) <- method_labels

df <- data.frame(
  method = factor(method_labels, levels = method_labels),
  x = seq_along(method_labels), y = 1,
  stringsAsFactors = FALSE
)

# 2 行で並べる (2x2 layout when n_m == 4)
p <- ggplot(df, aes(x, y, fill = method)) +
  geom_point(size = 6, shape = 21, colour = "white", stroke = 0.5) +
  scale_fill_manual(values = colors, name = NULL,
                    guide = guide_legend(nrow = 2, byrow = TRUE,
                                         override.aes = list(size = 7))) +
  theme_void(base_size = 16) +
  theme(
    legend.position = "bottom",
    legend.key = element_rect(fill = "transparent", colour = NA),
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.text = element_text(size = 16),
    plot.background  = element_rect(fill = "transparent", colour = NA)
  )

# Extract legend only via cowplot
legend_grob <- cowplot::get_legend(p)
legend_plot <- cowplot::ggdraw(legend_grob) +
  theme(plot.background = element_rect(fill = "transparent", colour = NA))

# Aspect: width based on ceil(n/2) columns, height for 2 rows
n_col <- ceiling(n_m / 2)
w <- max(4, 2.5 * n_col)
ggsave(outfile, legend_plot,
       width = w, height = 1.8, dpi = 150, bg = "transparent")
cat("Saved:", outfile, "\n")
