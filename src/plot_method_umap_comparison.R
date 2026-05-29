# =============================================================================
# 各手法の embeddings.csv から UMAP を作り、横並びで比較する
# =============================================================================
# Usage: Rscript src/plot_method_umap_comparison.R \
#   [--meta=<rna_meta.csv,atac_meta.csv>] \
#   <method1=path1.csv> <method2=path2.csv> ... <dataset_name> <outdir>
#
# --meta が与えられたら、その metadata の broad_lineage / germ_layer /
# Germ.layer 列を cell_id で join して、combined panel に行を追加する。
# =============================================================================

suppressPackageStartupMessages({
  library(uwot)
  library(ggplot2)
  library(dplyr)
  library(RColorBrewer)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)

# Optional --meta=<rna_meta.csv,atac_meta.csv>
meta_paths <- character(0)
meta_idx <- grep("^--meta=", args)
if (length(meta_idx)) {
  meta_paths <- strsplit(sub("^--meta=", "", args[meta_idx]), ",", fixed = TRUE)[[1]]
  args <- args[-meta_idx]
}

n <- length(args)
dataset_name <- args[n - 1]
outdir       <- args[n]
pairs        <- args[seq_len(n - 2)]

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("=== Method comparison UMAP — %s ===\n", dataset_name))

methods <- c(); paths <- c()
for (p in pairs) {
  kv <- strsplit(p, "=", fixed = TRUE)[[1]]
  stopifnot(length(kv) == 2)
  methods <- c(methods, kv[1]); paths <- c(paths, kv[2])
}

# Load broad_lineage lookup if --meta given
broad_lookup <- NULL
if (length(meta_paths)) {
  parts <- lapply(meta_paths, function(mf) {
    if (!file.exists(mf)) return(NULL)
    md <- read.csv(mf, stringsAsFactors = FALSE)
    bl_col <- intersect(c("broad_lineage", "germ_layer", "Germ.layer"),
                        colnames(md))
    if (length(bl_col) == 0) return(NULL)
    id_col <- intersect(c("cell_id", "cell.id"), colnames(md))
    if (length(id_col) == 0) id_col <- colnames(md)[1]
    data.frame(cell_id = md[[id_col[1]]],
               broad_lineage = md[[bl_col[1]]],
               stringsAsFactors = FALSE)
  })
  broad_lookup <- do.call(rbind, parts)
  cat(sprintf("broad_lineage lookup: %d cells\n", nrow(broad_lookup)))
}

# Compute UMAP per method
umaps <- list()
common_ct <- NULL
for (i in seq_along(methods)) {
  m <- methods[i]; pth <- paths[i]
  cat(sprintf("Loading %s from %s ... ", m, pth))
  if (!file.exists(pth)) { cat("MISSING; skip\n"); next }
  d <- read.csv(pth, stringsAsFactors = FALSE)
  meta_cols <- intersect(c("cell_id", "modality", "celltype"), colnames(d))
  emb <- as.matrix(d[, setdiff(colnames(d), meta_cols), drop = FALSE])
  cat(sprintf("%d x %d\n", nrow(emb), ncol(emb)))

  set.seed(42)
  um <- umap(emb, n_neighbors = 30, min_dist = 0.3,
             metric = "cosine", n_threads = 4)
  df_m <- data.frame(
    UMAP1 = um[, 1], UMAP2 = um[, 2],
    modality = d$modality, celltype = d$celltype, method = m,
    stringsAsFactors = FALSE
  )
  if (!is.null(broad_lookup) && "cell_id" %in% colnames(d)) {
    df_m$broad_lineage <-
      broad_lookup$broad_lineage[match(d$cell_id, broad_lookup$cell_id)]
  } else {
    df_m$broad_lineage <- NA_character_
  }
  umaps[[m]] <- df_m
  common_ct <- if (is.null(common_ct)) unique(d$celltype) else
    union(common_ct, unique(d$celltype))
}

if (length(umaps) == 0) stop("No embeddings found.")

all_df <- bind_rows(umaps)
all_df$method <- factor(all_df$method, levels = methods)
all_df <- all_df[!is.na(all_df$celltype), ]
have_broad <- !all(is.na(all_df$broad_lineage))

ct_levels <- sort(unique(all_df$celltype))
n_ct <- length(ct_levels)
ct_pal <- if (n_ct <= 8) brewer.pal(max(3, n_ct), "Set2")[seq_len(n_ct)] else
  colorRampPalette(brewer.pal(8, "Set2"))(n_ct)
names(ct_pal) <- ct_levels
mod_pal <- c(RNA = "#1f78b4", ATAC = "#e31a1c")

base_theme <- theme_void() +
  theme(
    legend.position = "none",
    strip.text = element_blank(),
    plot.background  = element_rect(fill = "transparent", colour = NA),
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.title       = element_blank()
  )

# Panel 1: by modality
p_mod <- ggplot(all_df, aes(UMAP1, UMAP2, colour = modality)) +
  geom_point(size = 0.3, alpha = 0.55) +
  facet_wrap(~method, nrow = 1, scales = "free") +
  scale_colour_manual(values = mod_pal, name = "Modality") +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  base_theme +
  labs(title = sprintf("%s — UMAP by modality across methods", dataset_name))

# Panel 2: by celltype
p_ct <- ggplot(all_df, aes(UMAP1, UMAP2, colour = celltype)) +
  geom_point(size = 0.3, alpha = 0.55) +
  facet_wrap(~method, nrow = 1, scales = "free") +
  scale_colour_manual(values = ct_pal, name = "Cell type") +
  guides(colour = guide_legend(ncol = 1, override.aes = list(size = 2, alpha = 1))) +
  base_theme +
  labs(title = sprintf("%s — UMAP by cell type across methods", dataset_name))

# --- vertical legend helper (transparent bg, same colour as the panel) ---
save_legend <- function(values, name, outfile) {
  df <- data.frame(label = factor(names(values), levels = names(values)),
                   x = seq_along(values), y = 1)
  ph <- ggplot(df, aes(x, y, fill = label)) +
    geom_point(size = 6, shape = 21, colour = "white", stroke = 0.5) +
    scale_fill_manual(values = values, name = name,
                      guide = guide_legend(ncol = 1, byrow = TRUE,
                                           override.aes = list(size = 7))) +
    theme_void(base_size = 16) +
    theme(
      legend.position = "right",
      legend.key = element_rect(fill = "transparent", colour = NA),
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 15, face = "bold"),
      plot.background = element_rect(fill = "transparent", colour = NA)
    )
  legend_grob <- cowplot::get_legend(ph)
  legend_plot <- cowplot::ggdraw(legend_grob) +
    theme(plot.background = element_rect(fill = "transparent", colour = NA))
  # Vertical layout: width fixed, height scales with item count
  h_out <- max(1.5, 0.45 * length(values) + 0.6)
  ggsave(outfile, legend_plot,
         width = 4, height = h_out, dpi = 150, bg = "transparent")
}
suppressPackageStartupMessages(library(cowplot))
save_legend(mod_pal, "Modality",
            file.path(outdir, "umap_legend_modality.png"))
save_legend(ct_pal,  "Cell type",
            file.path(outdir, "umap_legend_celltype.png"))

w <- max(12, 6 * length(methods))   # 横方向を広めに (4 → 6 per method)
ggsave(file.path(outdir, "umap_methods_by_modality.png"), p_mod,
       width = w, height = 5, dpi = 150, bg = "transparent")
ggsave(file.path(outdir, "umap_methods_by_celltype.png"), p_ct,
       width = w, height = 6, dpi = 150, bg = "transparent")

if (have_broad) {
  bl_df <- all_df[!is.na(all_df$broad_lineage), ]
  bl_levels <- sort(unique(bl_df$broad_lineage))
  n_bl <- length(bl_levels)
  bl_pal <- if (n_bl <= 8) brewer.pal(max(3, n_bl), "Dark2")[seq_len(n_bl)] else
    colorRampPalette(brewer.pal(8, "Dark2"))(n_bl)
  names(bl_pal) <- bl_levels
  save_legend(bl_pal, "Germlayer",
              file.path(outdir, "umap_legend_germlayer.png"))
  p_bl <- ggplot(bl_df, aes(UMAP1, UMAP2, colour = broad_lineage)) +
    geom_point(size = 0.3, alpha = 0.55) +
    facet_wrap(~method, nrow = 1, scales = "free") +
    scale_colour_manual(values = bl_pal, name = "Germlayer / broad lineage") +
    guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    base_theme +
    labs(title = sprintf("%s — UMAP by germlayer across methods", dataset_name))
  ggsave(file.path(outdir, "umap_methods_by_germlayer.png"), p_bl,
         width = w, height = 5, dpi = 150, bg = "transparent")
  combined_plot <- (p_mod / p_bl / p_ct) &
    theme(plot.background = element_rect(fill = "transparent", colour = NA))
  ggsave(file.path(outdir, "umap_methods_combined.png"), combined_plot,
         width = w, height = 16, dpi = 150, bg = "transparent")
} else {
  combined_plot <- (p_mod / p_ct) &
    theme(plot.background = element_rect(fill = "transparent", colour = NA))
  ggsave(file.path(outdir, "umap_methods_combined.png"), combined_plot,
         width = w, height = 11, dpi = 150, bg = "transparent")
}

cat("Saved UMAP comparison panels to:", outdir, "\n")
cat("Done.\n")
