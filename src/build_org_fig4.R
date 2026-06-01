# =============================================================================
# Build all Fig 4-equivalent panels for Organogenesis under the
#   "gPLS_stage as guide -> germ_layer as eval label" framing.
# Outputs directly into plot/Figures/supplementary/Organogenesis/ with the
# SuppFigS6_organogenesis_* prefix. Does not touch the Snakemake pipeline.
# =============================================================================

suppressPackageStartupMessages({
  library(uwot); library(ggplot2); library(dplyr); library(tidyr)
  library(RColorBrewer); library(patchwork); library(RANN); library(cowplot)
  library(Matrix)
})

set.seed(42)

OUTDIR <- "plot/Figures/supplementary/Organogenesis"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# --- metadata ---
rna  <- read.csv("output/Organogenesis/preprocess/rna_metadata.csv", stringsAsFactors = FALSE)
atac <- read.csv("output/Organogenesis/preprocess/atac_metadata.csv", stringsAsFactors = FALSE)
ct2gl <- with(rna, tapply(germ_layer, celltype,
  function(v) names(sort(table(v), decreasing = TRUE))[1]))
true_gl <- atac$Germ.layer

# --- methods + paths (only clean methods: gPLS_stage + 3 baselines) ---
methods <- c("gPLS", "Seurat", "Harmony", "Scanorama")
paths_pred <- list(
  "gPLS"      = "output/Organogenesis/guidedpls/stage/predicted_labels.csv",
  "Seurat"    = "output/Organogenesis/comparison/seurat/predicted_labels.csv",
  "Harmony"   = "output/Organogenesis/comparison/harmony/predicted_labels.csv",
  "Scanorama" = "output/Organogenesis/comparison/scanorama/predicted_labels.csv"
)
paths_emb <- list(
  "gPLS"      = "output/Organogenesis/guidedpls/stage/embeddings.csv",
  "Seurat"    = "output/Organogenesis/comparison/seurat/embeddings.csv",
  "Harmony"   = "output/Organogenesis/comparison/harmony/embeddings.csv",
  "Scanorama" = "output/Organogenesis/comparison/scanorama/embeddings.csv"
)
paths_bench <- list(
  "gPLS"      = "output/Organogenesis/benchmark/guidedpls_stage.tsv",
  "Seurat"    = "output/Organogenesis/benchmark/seurat.tsv",
  "Harmony"   = "output/Organogenesis/benchmark/harmony.tsv",
  "Scanorama" = "output/Organogenesis/benchmark/scanorama.tsv"
)

# Dark2 colors
colors <- brewer.pal(max(3, length(methods)), "Dark2")[seq_along(methods)]
names(colors) <- methods

# --- metric calculation (germ_layer eval) ---
acc_v <- setNames(rep(NA_real_, length(methods)), methods)
per_class <- list()
for (m in methods) {
  p <- read.csv(paths_pred[[m]], stringsAsFactors = FALSE)
  pred_gl <- ct2gl[p$predicted_celltype]
  acc_v[m] <- mean(pred_gl == true_gl, na.rm = TRUE)
  for (gl in sort(unique(na.omit(true_gl)))) {
    tp <- sum(pred_gl == gl & true_gl == gl, na.rm = TRUE)
    fp <- sum(pred_gl == gl & true_gl != gl, na.rm = TRUE)
    fn <- sum(pred_gl != gl & true_gl == gl, na.rm = TRUE)
    pr <- ifelse((tp + fp) == 0, 0, tp / (tp + fp))
    re <- ifelse((tp + fn) == 0, 0, tp / (tp + fn))
    f1 <- ifelse((pr + re) == 0, 0, 2 * pr * re / (pr + re))
    per_class[[length(per_class) + 1]] <- data.frame(
      method = m, germ_layer = gl, f1 = f1, stringsAsFactors = FALSE)
  }
}
df_pc <- do.call(rbind, per_class)
df_pc$method <- factor(df_pc$method, levels = methods)

# Benchmark
time_v <- setNames(rep(NA_real_, length(methods)), methods)
mem_v  <- setNames(rep(NA_real_, length(methods)), methods)
for (m in methods) {
  f <- paths_bench[[m]]
  if (file.exists(f)) {
    d <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    time_v[m] <- d[["s"]][1]
    mem_v[m]  <- d[["max_rss"]][1]
  }
}
df_bar <- data.frame(method = factor(methods, levels = methods),
                     accuracy = acc_v, time_sec = time_v, mem_gb = mem_v / 1024)

# --- bar themes / panel builder ---
base_theme <- theme_minimal(base_size = 36) +
  theme(
    axis.text.x  = element_blank(), axis.ticks.x = element_blank(),
    axis.text.y  = element_text(size = 36),
    axis.title.y = element_text(size = 42),
    plot.title   = element_blank(),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    plot.background  = element_rect(fill = "transparent", colour = NA),
    panel.background = element_rect(fill = "transparent", colour = NA)
  )
mk_bar <- function(y_col, y_lab, fmt) {
  ggplot(df_bar, aes(method, .data[[y_col]], fill = method)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = sprintf(fmt, .data[[y_col]])),
              vjust = -0.4, size = 11) +
    scale_fill_manual(values = colors) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    base_theme + labs(x = NULL, y = y_lab)
}
ggsave(file.path(OUTDIR, "SuppFigS6_organogenesis_method_accuracy.png"),
       mk_bar("accuracy", "Germ-layer accuracy", "%.3f"),
       width = 8, height = 6, dpi = 150, bg = "transparent")
ggsave(file.path(OUTDIR, "SuppFigS6_organogenesis_method_time.png"),
       mk_bar("time_sec", "Wall time (s)", "%.0f s"),
       width = 8, height = 6, dpi = 150, bg = "transparent")
ggsave(file.path(OUTDIR, "SuppFigS6_organogenesis_method_memory.png"),
       mk_bar("mem_gb", "Peak memory (GB)", "%.1f GB"),
       width = 8, height = 6, dpi = 150, bg = "transparent")

# --- method legend (2x2) ---
{
  ldf <- data.frame(label = factor(methods, levels = methods),
                    x = seq_along(methods), y = 1)
  pleg <- ggplot(ldf, aes(x, y, fill = label)) +
    geom_point(size = 6, shape = 21, colour = "white", stroke = 0.5) +
    scale_fill_manual(values = colors, name = NULL,
                      guide = guide_legend(nrow = 2, byrow = TRUE,
                                           override.aes = list(size = 7))) +
    theme_void(base_size = 16) +
    theme(legend.position = "bottom",
          legend.key = element_rect(fill = "transparent", colour = NA),
          legend.background = element_rect(fill = "transparent", colour = NA),
          legend.text = element_text(size = 16),
          plot.background = element_rect(fill = "transparent", colour = NA))
  gleg <- cowplot::ggdraw(cowplot::get_legend(pleg)) +
    theme(plot.background = element_rect(fill = "transparent", colour = NA))
  ggsave(file.path(OUTDIR, "SuppFigS6_organogenesis_method_legend.png"),
         gleg, width = 5, height = 1.8, dpi = 150, bg = "transparent")
}

# --- per-class F1 heatmap (6 germ_layer x 4 methods) ---
{
  df_pc$method <- factor(df_pc$method, levels = methods)
  gl_order <- df_pc %>% group_by(germ_layer) %>% summarise(m = mean(f1)) %>%
    arrange(m) %>% pull(germ_layer)
  df_pc$germ_layer <- factor(df_pc$germ_layer, levels = gl_order)
  ph <- ggplot(df_pc, aes(method, germ_layer, fill = f1)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.2f", f1)),
              colour = ifelse(df_pc$f1 > 0.5, "white", "grey20"), size = 7) +
    scale_fill_gradient2(low = "#f7fbff", mid = "#6baed6", high = "#08306b",
                         midpoint = 0.5, limits = c(0, 1), name = "F1") +
    theme_minimal(base_size = 22) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          panel.grid = element_blank(),
          plot.title = element_blank(),
          plot.background = element_rect(fill = "transparent", colour = NA),
          panel.background = element_rect(fill = "transparent", colour = NA)) +
    labs(x = NULL, y = NULL)
  ggsave(file.path(OUTDIR, "SuppFigS6_organogenesis_per_class_f1_heatmap.png"),
         ph, width = 9, height = 5, dpi = 150, bg = "transparent")
}

# --- UMAP combined (3 rows × 4 cols) + 3 legends (modality / stage / germ) ---
umap_theme <- theme_void() +
  theme(legend.position = "none", strip.text = element_blank(),
        plot.background  = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.title = element_blank())

# Build per-method data frame with UMAP coords
umaps <- list()
for (m in methods) {
  d <- read.csv(paths_emb[[m]], stringsAsFactors = FALSE)
  meta_cols <- intersect(c("cell_id","modality","celltype"), colnames(d))
  emb <- as.matrix(d[, setdiff(colnames(d), meta_cols), drop = FALSE])
  set.seed(42)
  um <- umap(emb, n_neighbors = 30, min_dist = 0.3, metric = "cosine",
             n_threads = 4)
  # Join stage + germ_layer from atac/rna metadata
  rid <- match(d$cell_id, rna$cell_id)
  aid <- match(d$cell_id, atac$cell_id)
  stage <- ifelse(!is.na(rid), rna$Stage[rid],
           ifelse(!is.na(aid), atac$Stage[aid], NA))
  gl    <- ifelse(!is.na(rid), rna$germ_layer[rid],
           ifelse(!is.na(aid), atac$Germ.layer[aid], NA))
  umaps[[m]] <- data.frame(
    UMAP1 = um[, 1], UMAP2 = um[, 2], modality = d$modality,
    stage = stage, germ_layer = gl, method = m,
    stringsAsFactors = FALSE)
  cat(sprintf("UMAP %s: n=%d, stage_join=%d, gl_join=%d\n",
              m, nrow(d), sum(!is.na(stage)), sum(!is.na(gl))))
}
all_df <- bind_rows(umaps)
all_df$method <- factor(all_df$method, levels = methods)
w <- 6 * length(methods)

# Modality palette
mod_pal <- c(RNA = "#1f78b4", ATAC = "#e31a1c")
p_mod <- ggplot(all_df, aes(UMAP1, UMAP2, colour = modality)) +
  geom_point(size = 0.3, alpha = 0.55) +
  facet_wrap(~method, nrow = 1, scales = "free") +
  scale_colour_manual(values = mod_pal) + umap_theme

# Stage palette
st_levels <- sort(unique(all_df$stage[!is.na(all_df$stage)]))
st_pal <- brewer.pal(max(3, length(st_levels)), "YlOrRd")[seq_along(st_levels)]
names(st_pal) <- st_levels
p_st <- ggplot(all_df[!is.na(all_df$stage), ],
               aes(UMAP1, UMAP2, colour = stage)) +
  geom_point(size = 0.3, alpha = 0.55) +
  facet_wrap(~method, nrow = 1, scales = "free") +
  scale_colour_manual(values = st_pal) + umap_theme

# Germ layer palette
gl_levels <- sort(unique(all_df$germ_layer[!is.na(all_df$germ_layer)]))
gl_pal <- brewer.pal(max(3, length(gl_levels)), "Dark2")[seq_along(gl_levels)]
names(gl_pal) <- gl_levels
p_gl <- ggplot(all_df[!is.na(all_df$germ_layer), ],
               aes(UMAP1, UMAP2, colour = germ_layer)) +
  geom_point(size = 0.3, alpha = 0.55) +
  facet_wrap(~method, nrow = 1, scales = "free") +
  scale_colour_manual(values = gl_pal) + umap_theme

combined <- (p_mod / p_st / p_gl) &
  theme(plot.background = element_rect(fill = "transparent", colour = NA))
ggsave(file.path(OUTDIR, "SuppFigS6_organogenesis_method_umap.png"),
       combined, width = w, height = 16, dpi = 150, bg = "transparent")

# Standalone vertical legends
save_v_legend <- function(values, name, file) {
  df <- data.frame(label = factor(names(values), levels = names(values)),
                   x = seq_along(values), y = 1)
  pl <- ggplot(df, aes(x, y, fill = label)) +
    geom_point(size = 6, shape = 21, colour = "white", stroke = 0.5) +
    scale_fill_manual(values = values, name = name,
                      guide = guide_legend(ncol = 1, byrow = TRUE,
                                           override.aes = list(size = 7))) +
    theme_void(base_size = 16) +
    theme(legend.position = "right",
          legend.key = element_rect(fill = "transparent", colour = NA),
          legend.background = element_rect(fill = "transparent", colour = NA),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 15, face = "bold"),
          plot.background = element_rect(fill = "transparent", colour = NA))
  g <- cowplot::ggdraw(cowplot::get_legend(pl)) +
    theme(plot.background = element_rect(fill = "transparent", colour = NA))
  ggsave(file, g, width = 4,
         height = max(1.5, 0.45 * length(values) + 0.6),
         dpi = 150, bg = "transparent")
}
save_v_legend(mod_pal, "Modality",
              file.path(OUTDIR, "SuppFigS6_organogenesis_umap_legend_modality.png"))
save_v_legend(st_pal,  "Stage",
              file.path(OUTDIR, "SuppFigS6_organogenesis_umap_legend_stage.png"))
save_v_legend(gl_pal,  "Germ layer",
              file.path(OUTDIR, "SuppFigS6_organogenesis_umap_legend_germlayer.png"))

cat("\n=== Summary ===\n")
cat(sprintf("%-12s %10s %10s %10s\n", "method", "acc", "time(s)", "mem(GB)"))
for (i in seq_along(methods)) {
  cat(sprintf("%-12s %10.3f %10.1f %10.2f\n",
              methods[i], df_bar$accuracy[i], df_bar$time_sec[i],
              df_bar$mem_gb[i]))
}
cat("\nDone.\n")
