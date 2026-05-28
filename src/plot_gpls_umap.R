#  =============================================================================
# guided-PLS の latent score を UMAP で可視化
# =============================================================================
# Usage: Rscript src/plot_gpls_umap.R \
#   <guidedpls.RData> <rna_metadata.csv> <atac_metadata.csv> \
#   <predicted_labels.csv> <dataset> <condition> <outdir>
# =============================================================================

suppressPackageStartupMessages({
  library(uwot)
  library(ggplot2)
  library(dplyr)
  library(RColorBrewer)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
rdata_file    <- args[1]
rna_meta_file <- args[2]
atac_meta_file <- args[3]
pred_file     <- args[4]
dataset_name  <- args[5]
condition     <- args[6]
outdir        <- args[7]

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("=== gPLS UMAP: %s / %s ===\n", dataset_name, condition))

load(rdata_file)
score1 <- out$scoreX1
score2 <- out$scoreX2

rna_meta  <- read.csv(rna_meta_file,  stringsAsFactors = FALSE)
atac_meta <- read.csv(atac_meta_file, stringsAsFactors = FALSE)
preds     <- read.csv(pred_file,      stringsAsFactors = FALSE)

stopifnot(nrow(score1) == nrow(rna_meta))
stopifnot(nrow(score2) == nrow(atac_meta))
stopifnot(nrow(score2) == nrow(preds))

n1 <- nrow(score1); n2 <- nrow(score2)
cat(sprintf("RNA: %d cells, ATAC: %d cells, latent dim: %d\n", n1, n2, ncol(score1)))

# Combine and run UMAP
combined <- rbind(score1, score2)
set.seed(42)
um <- umap(combined, n_neighbors = 30, min_dist = 0.3, metric = "cosine", n_threads = 4)

# Build a single tidy data frame
stage_col <- function(meta) {
  for (cand in c("Stage", "stage", "development_stage")) {
    if (cand %in% colnames(meta)) return(as.character(meta[[cand]]))
  }
  rep(NA_character_, nrow(meta))
}

df <- data.frame(
  UMAP1 = um[, 1],
  UMAP2 = um[, 2],
  modality = c(rep("RNA", n1), rep("ATAC", n2)),
  celltype_true = c(rna_meta$celltype, atac_meta$celltype),
  stage = c(stage_col(rna_meta), stage_col(atac_meta)),
  predicted = c(rep(NA_character_, n1), preds$predicted_celltype),
  stringsAsFactors = FALSE
)

# Drop cells with no celltype to keep palette clean
df_ct <- df[!is.na(df$celltype_true), ]

# A consistent ordered palette for celltypes
ct_levels <- sort(unique(df_ct$celltype_true))
n_ct <- length(ct_levels)
ct_pal <- if (n_ct <= 8) {
  brewer.pal(max(3, n_ct), "Set2")[seq_len(n_ct)]
} else {
  colorRampPalette(brewer.pal(8, "Set2"))(n_ct)
}
names(ct_pal) <- ct_levels

base_theme <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2, colour = "grey90"),
    legend.position = "right",
    legend.key.size = grid::unit(0.4, "cm"),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    plot.title = element_text(size = 11, face = "bold")
  )

# --- Panel 1: by modality ---
mod_pal <- c(RNA = "#1f78b4", ATAC = "#e31a1c")
p_mod <- ggplot(df, aes(UMAP1, UMAP2, colour = modality)) +
  geom_point(size = 0.4, alpha = 0.6) +
  scale_colour_manual(values = mod_pal, name = "Modality") +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  base_theme +
  labs(title = sprintf("%s / Z = %s — modality", dataset_name, condition))

# --- Panel 2: by true celltype, faceted by modality ---
p_ct <- ggplot(df_ct, aes(UMAP1, UMAP2, colour = celltype_true)) +
  geom_point(size = 0.4, alpha = 0.6) +
  facet_wrap(~modality) +
  scale_colour_manual(values = ct_pal, name = "Cell type (eval)") +
  guides(colour = guide_legend(ncol = 1, override.aes = list(size = 2, alpha = 1))) +
  base_theme +
  labs(title = sprintf("%s / Z = %s — true celltype (eval label)", dataset_name, condition))

# --- Panel 3: by stage ---
has_stage <- !all(is.na(df$stage))
if (has_stage) {
  stage_levels <- sort(unique(df$stage[!is.na(df$stage)]))
  stage_pal <- colorRampPalette(brewer.pal(min(9, max(3, length(stage_levels))), "YlOrRd"))(length(stage_levels))
  names(stage_pal) <- stage_levels
  p_stage <- ggplot(df[!is.na(df$stage), ], aes(UMAP1, UMAP2, colour = stage)) +
    geom_point(size = 0.4, alpha = 0.6) +
    facet_wrap(~modality) +
    scale_colour_manual(values = stage_pal, name = "Stage") +
    guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    base_theme +
    labs(title = sprintf("%s / Z = %s — developmental stage", dataset_name, condition))
}

# --- Panel 4: ATAC colored by predicted celltype ---
df_atac_pred <- df[df$modality == "ATAC" & !is.na(df$predicted), ]
df_atac_pred$correct <- df_atac_pred$predicted == df_atac_pred$celltype_true
p_pred <- ggplot(df_atac_pred, aes(UMAP1, UMAP2, colour = predicted)) +
  geom_point(size = 0.4, alpha = 0.6) +
  scale_colour_manual(values = ct_pal, name = "Predicted celltype") +
  guides(colour = guide_legend(ncol = 1, override.aes = list(size = 2, alpha = 1))) +
  base_theme +
  labs(title = sprintf("%s / Z = %s — ATAC predicted celltype", dataset_name, condition))

# --- Save individual panels ---
ggsave(file.path(outdir, "umap_by_modality.png"),  p_mod,   width = 7,  height = 5, dpi = 150)
ggsave(file.path(outdir, "umap_by_celltype.png"),  p_ct,    width = 12, height = 5, dpi = 150)
if (has_stage) {
  ggsave(file.path(outdir, "umap_by_stage.png"),   p_stage, width = 12, height = 5, dpi = 150)
}
ggsave(file.path(outdir, "umap_atac_predicted.png"), p_pred, width = 7, height = 5, dpi = 150)

# --- Combined panel: 4 列横並びで横幅広く (旧: 2x2) ---
if (has_stage) {
  combined_plot <- (p_mod | p_stage | p_ct | p_pred) +
    plot_annotation(
      title = sprintf("guided-PLS latent space — %s (Z = %s)", dataset_name, condition),
      theme = theme(plot.title = element_text(size = 13, face = "bold"))
    )
  ggsave(file.path(outdir, "umap_combined.png"), combined_plot,
    width = 26, height = 6, dpi = 150)
} else {
  combined_plot <- (p_mod | p_ct | p_pred) +
    plot_annotation(
      title = sprintf("guided-PLS latent space — %s (Z = %s)", dataset_name, condition),
      theme = theme(plot.title = element_text(size = 13, face = "bold"))
    )
  ggsave(file.path(outdir, "umap_combined.png"), combined_plot,
    width = 22, height = 6, dpi = 150)
}

cat("UMAP plots saved to:", outdir, "\n")
cat("Done.\n")
