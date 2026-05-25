# =============================================================================
# scRNA-Seq データ統合 (可変個数のRDSファイルをマージ)
# =============================================================================
# Usage: Rscript src/merge_rna.R \
#   <rna1.rds> [<rna2.rds> ...] \
#   <out_rds> <out_meta_csv> \
#   <subsample_n> <n_hvg> <seed>
#
# 末尾3引数は固定パラメータ。それ以外は入力RDSファイル。
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)

# 末尾3引数がパラメータ、残りがRDSファイル
n_args <- length(args)
seed        <- as.integer(args[n_args])
n_hvg       <- as.integer(args[n_args - 1])
subsample_n <- as.integer(args[n_args - 2])
out_meta_csv <- args[n_args - 3]
out_rds      <- args[n_args - 4]
rds_files    <- args[1:(n_args - 5)]

set.seed(seed)
cat("=== Merging RNA datasets ===\n")
cat(sprintf("Input RDS files: %s\n", paste(rds_files, collapse = ", ")))

# --- 1. 全RDSを読み込み ---
obj_list <- lapply(rds_files, function(f) {
  cat(sprintf("Reading %s ...\n", f))
  readRDS(f)
})

for (i in seq_along(obj_list)) {
  cat(sprintf("  %s: %d cells\n", basename(rds_files[i]), ncol(obj_list[[i]])))
}

# --- 2. マージ ---
if (length(obj_list) == 1) {
  rna <- obj_list[[1]]
} else {
  rna <- merge(obj_list[[1]], obj_list[-1])
}
cat(sprintf("Merged: %d cells\n", ncol(rna)))

# --- 3. サブサンプリング（ステージ層別） ---
if (ncol(rna) > subsample_n) {
  cat(sprintf("Subsampling to %d cells (stratified by stage) ...\n", subsample_n))
  cells_keep <- c()
  stage_table <- table(rna$Stage)
  for (stg in names(stage_table)) {
    stg_cells <- colnames(rna)[rna$Stage == stg]
    n_take <- round(as.numeric(subsample_n) * as.numeric(stage_table[stg]) / as.numeric(sum(stage_table)))
    n_take <- min(n_take, length(stg_cells))
    cells_keep <- c(cells_keep, sample(stg_cells, n_take))
  }
  rna <- subset(rna, cells = cells_keep)
}
cat(sprintf("After subsampling: %d cells\n", ncol(rna)))

# --- 4. HVG選択・スケーリング ---
rna <- FindVariableFeatures(rna, nfeatures = n_hvg)
rna <- ScaleData(rna)

# --- 5. メタデータにgerm_layerを追加（celltypeから導出） ---
celltype_to_germlayer <- c(
  # Mesoderm
  "First heart field"             = "Mesoderm",
  "Second heart field"            = "Mesoderm",
  "Intermediate mesoderm"         = "Mesoderm",
  "Paraxial mesoderm A"           = "Mesoderm",
  "Paraxial mesoderm B"           = "Mesoderm",
  "Somatic mesoderm"              = "Mesoderm",
  "Splanchnic mesoderm"           = "Mesoderm",
  "Lateral plate mesoderm"        = "Mesoderm",
  "Amniochorionic mesoderm A"     = "Mesoderm",
  "Amniochorionic mesoderm B"     = "Mesoderm",
  "Extraembryonic mesoderm"       = "Mesoderm",
  "Mesenchymal stromal cells"     = "Mesoderm",
  "Allantois"                     = "Mesoderm",
  "Cardiomyocytes"                = "Mesoderm",
  "Endothelium"                   = "Mesoderm",
  "Blood progenitors"             = "Mesoderm",
  "Primitive erythroid cells"     = "Mesoderm",
  "Erythroid cells"               = "Mesoderm",
  "Megakaryocytes"                = "Mesoderm",
  "Hematopoietic progenitors"     = "Mesoderm",
  "Neuromesodermal progenitors"   = "Mesoderm",
  "Limb mesoderm"                 = "Mesoderm",
  "Jaw and tooth progenitors"     = "Mesoderm",
  "Myocytes"                      = "Mesoderm",
  "Osteoblast progenitors"        = "Mesoderm",
  "Chondrocyte progenitors"       = "Mesoderm",
  # Neuroectoderm
  "Forebrain/midbrain"            = "Neuroectoderm",
  "Hindbrain"                     = "Neuroectoderm",
  "Spinal cord"                   = "Neuroectoderm",
  "Neural crest"                  = "Neuroectoderm",
  "Diencephalon"                  = "Neuroectoderm",
  "Telencephalon"                 = "Neuroectoderm",
  "Midbrain"                      = "Neuroectoderm",
  "Cerebellum"                    = "Neuroectoderm",
  "Neural tube"                   = "Neuroectoderm",
  "Dorsal root ganglia"           = "Neuroectoderm",
  "Sympathetic neurons"           = "Neuroectoderm",
  "Cranial neural crest"          = "Neuroectoderm",
  "Schwann cell precursors"       = "Neuroectoderm",
  "Retina"                        = "Neuroectoderm",
  # Surface ectoderm
  "Pre-epidermal keratinocytes"   = "Surface ectoderm",
  "Placodal area"                 = "Surface ectoderm",
  "Epidermis"                     = "Surface ectoderm",
  "Limb ectoderm"                 = "Surface ectoderm",
  "Otic placode"                  = "Surface ectoderm",
  "Lens"                          = "Surface ectoderm",
  # Endoderm
  "Gut"                           = "Endoderm",
  "Liver"                         = "Endoderm",
  "Lung"                          = "Endoderm",
  "Pancreas"                      = "Endoderm",
  "Thyroid"                       = "Endoderm",
  "Extraembryonic visceral endoderm" = "Endoderm",
  "Definitive endoderm"           = "Endoderm",
  "Foregut"                       = "Endoderm",
  "Midgut/Hindgut"                = "Endoderm",
  "Pharyngeal endoderm"           = "Endoderm",
  # ExE embryo
  "ExE ectoderm"                  = "ExE embryo",
  "Trophoblast"                   = "ExE embryo",
  "ExE endoderm"                  = "ExE embryo",
  "Parietal endoderm"             = "ExE embryo"
)

rna@meta.data$germ_layer <- celltype_to_germlayer[rna@meta.data$celltype]
rna@meta.data$germ_layer[is.na(rna@meta.data$germ_layer)] <- "Unknown"

cat("Germ layer distribution (RNA):\n")
print(table(rna$germ_layer))

# --- 6. 保存 ---
saveRDS(rna, out_rds)

meta_out <- rna@meta.data
meta_out$cell_id <- rownames(meta_out)
write.csv(meta_out, out_meta_csv, row.names = FALSE)

cat("Done.\n")
