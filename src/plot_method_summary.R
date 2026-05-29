# =============================================================================
# 手法ごとの kNN accuracy / 計算時間 / ピークメモリを 3 つの独立 PNG として
# 保存する (Dark2、legend なし、透明背景)。
# =============================================================================
# Usage: Rscript src/plot_method_summary.R \
#   <metrics.csv> <benchmark_dir> <conditions> <methods> <dataset_name> <outdir>
# 出力: <outdir>/method_accuracy.png, method_time.png, method_memory.png
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(RColorBrewer)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
metrics_file   <- args[1]
bench_dir      <- args[2]
conditions_str <- args[3]
methods_str    <- args[4]
dataset_name   <- args[5]
outdir         <- args[6]

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

conditions   <- strsplit(conditions_str, ",")[[1]]
methods      <- strsplit(methods_str, ",")[[1]]
comp_methods <- methods[methods != "NONE" & methods != ""]

# Method order: gPLS variants first, then baselines
gpls_keys   <- paste0("guidedpls_", conditions)
gpls_labels <- if (length(conditions) <= 1) "gPLS" else
  paste0("gPLS (", conditions, ")")
method_keys   <- c(gpls_keys, comp_methods)
method_labels <- c(gpls_labels, comp_methods)

# kNN accuracy from metrics.csv
metrics <- read.csv(metrics_file, stringsAsFactors = FALSE)
acc_map <- setNames(metrics$accuracy, metrics$method)
acc_vec <- acc_map[method_keys]

# Benchmark TSVs
bench_files <- c(
  setNames(paste0(bench_dir, "/guidedpls_", conditions, ".tsv"), gpls_keys),
  setNames(paste0(bench_dir, "/", comp_methods, ".tsv"),         comp_methods)
)
time_vec <- rep(NA_real_, length(method_keys))
mem_vec  <- rep(NA_real_, length(method_keys))
names(time_vec) <- names(mem_vec) <- method_keys
for (k in method_keys) {
  f <- bench_files[[k]]
  if (!file.exists(f)) next
  d <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  time_vec[k] <- d[["s"]][1]
  mem_vec[k]  <- d[["max_rss"]][1]
}

# Drop methods that have any missing value
keep <- !is.na(acc_vec) & !is.na(time_vec) & !is.na(mem_vec)
df <- data.frame(
  method   = factor(method_labels[keep], levels = method_labels[keep]),
  accuracy = acc_vec[keep],
  time_sec = time_vec[keep],
  mem_gb   = mem_vec[keep] / 1024,
  stringsAsFactors = FALSE
)

# Dark2 colors for the methods
n_m <- nlevels(df$method)
colors <- if (n_m <= 8) brewer.pal(max(3, n_m), "Dark2")[seq_len(n_m)] else
  colorRampPalette(brewer.pal(8, "Dark2"))(n_m)
names(colors) <- levels(df$method)

base_theme <- theme_minimal(base_size = 36) +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y  = element_text(size = 36),
    axis.title.y = element_text(size = 42),
    plot.title   = element_blank(),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    plot.background  = element_rect(fill = "transparent", colour = NA),
    panel.background = element_rect(fill = "transparent", colour = NA)
  )

p_acc <- ggplot(df, aes(method, accuracy, fill = method)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.3f", accuracy)),
            vjust = -0.4, size = 11) +
  scale_fill_manual(values = colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  base_theme +
  labs(x = NULL, y = "kNN accuracy")

p_time <- ggplot(df, aes(method, time_sec, fill = method)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.0f s", time_sec)),
            vjust = -0.4, size = 11) +
  scale_fill_manual(values = colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  base_theme +
  labs(x = NULL, y = "Wall time (s)")

p_mem <- ggplot(df, aes(method, mem_gb, fill = method)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.1f GB", mem_gb)),
            vjust = -0.4, size = 11) +
  scale_fill_manual(values = colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  base_theme +
  labs(x = NULL, y = "Peak memory (GB)")

ggsave(file.path(outdir, "method_accuracy.png"), p_acc,
       width = 8, height = 6, dpi = 150, bg = "transparent")
ggsave(file.path(outdir, "method_time.png"), p_time,
       width = 8, height = 6, dpi = 150, bg = "transparent")
ggsave(file.path(outdir, "method_memory.png"), p_mem,
       width = 8, height = 6, dpi = 150, bg = "transparent")
cat("Saved 3 files to:", outdir, "\n")
cat("\n=== summary ===\n")
print(df, row.names = FALSE)
