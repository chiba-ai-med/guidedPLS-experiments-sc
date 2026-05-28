# =============================================================================
# True → Predicted ラベル流れの alluvial プロット
# =============================================================================
# Usage: Rscript src/plot_label_flow.R \
#   <predicted_labels.csv> <atac_metadata.csv> <method_label> <outfile>
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggalluvial)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
pred_file   <- args[1]
atac_meta_file <- args[2]
method_label <- args[3]
outfile     <- args[4]

dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)

cat(sprintf("=== Label flow: %s ===\n", method_label))

preds <- read.csv(pred_file, stringsAsFactors = FALSE)
atac  <- read.csv(atac_meta_file, stringsAsFactors = FALSE)

stopifnot(nrow(preds) == nrow(atac))

df <- data.frame(
  true = atac$celltype,
  pred = preds$predicted_celltype,
  stringsAsFactors = FALSE
) |>
  filter(!is.na(true), !is.na(pred)) |>
  count(true, pred, name = "Freq")

# To keep the plot readable when there are many celltypes, drop links with
# very small counts (< 0.5% of the dataset). This keeps the dominant flows.
total <- sum(df$Freq)
df <- df %>% filter(Freq >= max(1, 0.005 * total))

# Color by whether the link is on-diagonal (correct) vs off-diagonal (mistake)
df$correct <- df$true == df$pred

# Order celltype levels by total true count, then apply to both axes
true_order <- df %>%
  group_by(true) %>% summarise(s = sum(Freq)) %>% arrange(desc(s)) %>% pull(true)
pred_order <- df %>%
  group_by(pred) %>% summarise(s = sum(Freq)) %>% arrange(desc(s)) %>% pull(pred)
all_levels <- union(true_order, pred_order)
df$true <- factor(df$true, levels = all_levels)
df$pred <- factor(df$pred, levels = all_levels)

p <- ggplot(df,
       aes(y = Freq, axis1 = true, axis2 = pred)) +
  geom_alluvium(aes(fill = correct), width = 1/8, alpha = 0.7) +
  geom_stratum(width = 1/8, fill = "grey90", colour = "grey40") +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum)),
            size = 2.4) +
  scale_x_discrete(limits = c("True", "Predicted"), expand = c(.05, .05)) +
  scale_fill_manual(values = c(`TRUE` = "#2c7fb8", `FALSE` = "#e6550d"),
                    labels = c(`TRUE` = "Correct", `FALSE` = "Misclassified"),
                    name = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    legend.position = "top",
    plot.title = element_text(face = "bold")
  ) +
  labs(title = sprintf("Label flow — %s", method_label),
       subtitle = "Links <0.5% of all cells removed for readability")

# Auto height by number of strata
h <- max(6, 0.18 * length(all_levels) + 3)
ggsave(outfile, p, width = 8, height = h, dpi = 150)
cat("Saved:", outfile, "\n")
cat("Done.\n")
