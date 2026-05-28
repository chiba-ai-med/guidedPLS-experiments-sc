# =============================================================================
# 手法ごとの計算時間 + メモリ使用量を Snakemake benchmark TSV から集約して
# 横並びの 2 panel bar plot にする
# =============================================================================
# Usage: Rscript src/plot_method_resources.R \
#   <benchmark_dir> <conditions> <methods> <dataset_name> <outfile>
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(RColorBrewer)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
bench_dir      <- args[1]
conditions_str <- args[2]
methods_str    <- args[3]
dataset_name   <- args[4]
outfile        <- args[5]

dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)

conditions   <- strsplit(conditions_str, ",")[[1]]
methods      <- strsplit(methods_str, ",")[[1]]
comp_methods <- methods[methods != "NONE" & methods != ""]

# Map method id → benchmark file path
all_runs <- c(
  setNames(paste0(bench_dir, "/guidedpls_", conditions, ".tsv"),
           paste0("gPLS (", conditions, ")")),
  setNames(paste0(bench_dir, "/", comp_methods, ".tsv"),
           comp_methods)
)

rows <- list()
for (label in names(all_runs)) {
  f <- all_runs[[label]]
  if (!file.exists(f)) { cat("MISSING:", f, "\n"); next }
  d <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  # snakemake benchmark schema: s, h:m:s, max_rss (MB), max_vms (MB), max_uss, max_pss, ...
  rows[[label]] <- data.frame(
    method = label,
    time_sec = d[["s"]][1],
    mem_mb   = d[["max_rss"]][1],
    stringsAsFactors = FALSE
  )
}
df <- do.call(rbind, rows)
df$method <- factor(df$method, levels = df$method)

# Color palette: gPLS shades + baselines Set2
n_gpls <- sum(grepl("^gPLS", levels(df$method)))
n_comp <- nlevels(df$method) - n_gpls
colors <- c(
  colorRampPalette(c("#2166AC","#92C5DE"))(max(1,n_gpls)),
  brewer.pal(max(3,n_comp), "Set2")[seq_len(n_comp)]
)
names(colors) <- levels(df$method)

base_theme <- theme_minimal(base_size = 18) +
  theme(
    axis.text.x  = element_text(angle = 30, hjust = 1, size = 14),
    axis.text.y  = element_text(size = 13),
    axis.title.y = element_text(size = 16),
    plot.title   = element_text(size = 17, face = "bold"),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

p_time <- ggplot(df, aes(method, time_sec, fill = method)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.0f s", time_sec)),
            vjust = -0.4, size = 5) +
  scale_fill_manual(values = colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  base_theme +
  labs(x = NULL, y = "Wall time (s)",
       title = sprintf("Computation time — %s", dataset_name))

p_mem <- ggplot(df, aes(method, mem_mb / 1024, fill = method)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.1f GB", mem_mb / 1024)),
            vjust = -0.4, size = 5) +
  scale_fill_manual(values = colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  base_theme +
  labs(x = NULL, y = "Peak memory (GB, max_rss)",
       title = sprintf("Peak memory — %s", dataset_name))

combined <- p_time + p_mem
ggsave(outfile, combined, width = 16, height = 6, dpi = 150)
cat("Saved:", outfile, "\n")
cat("\n=== resource summary ===\n")
print(df, row.names = FALSE)
