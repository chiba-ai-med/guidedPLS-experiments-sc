# =============================================================================
# scATAC-Seq 前処理: fragments → Gene Activity Matrix
# =============================================================================
# Usage: Rscript src/preprocess_atac.R \
#   <metadata> <fragments_dir> <out_rds> <out_meta_csv> \
#   <min_peaks> <subsample_n> <seed>
# =============================================================================

suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
  library(GenomicRanges)
  library(EnsDb.Mmusculus.v79)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
metadata_file  <- args[1]
fragments_dir  <- args[2]
out_rds        <- args[3]
out_meta_csv   <- args[4]
min_peaks      <- as.integer(args[5])
subsample_n    <- as.integer(args[6])
seed           <- as.integer(args[7])

set.seed(seed)
cat("=== scATAC-Seq preprocessing ===\n")

# --- 1. メタデータ読み込み ---
meta <- read.delim(metadata_file, stringsAsFactors = FALSE)
cat(sprintf("Metadata: %d cells\n", nrow(meta)))

meta <- meta[meta$Stage %in% c("E8.5", "E9.5", "E10.5"), ]
cat(sprintf("After stage filter: %d cells\n", nrow(meta)))

meta <- meta[meta$DoubletScore == 0, ]
cat(sprintf("After doublet filter: %d cells\n", nrow(meta)))
rownames(meta) <- meta$cell.id

# --- 2. サブサンプリング（先に行い処理量を削減）---
if (nrow(meta) > subsample_n) {
  cat(sprintf("Subsampling metadata to %d cells ...\n", subsample_n))
  cells_keep <- c()
  stage_table <- table(meta$Stage)
  for (stg in names(stage_table)) {
    stg_cells <- meta$cell.id[meta$Stage == stg]
    n_take <- round(subsample_n * stage_table[stg] / sum(stage_table))
    n_take <- min(n_take, length(stg_cells))
    cells_keep <- c(cells_keep, sample(stg_cells, n_take))
  }
  meta <- meta[cells_keep, ]
}
cat(sprintf("Working with %d cells\n", nrow(meta)))

# --- 3. Fragmentsファイル ---
frag_files <- list.files(fragments_dir,
  pattern = "Mouse_embryo_E.+\\.fragments\\.tsv\\.gz$",
  full.names = TRUE)
cat(sprintf("Found %d fragment files\n", length(frag_files)))

# tabixインデックス確認
for (f in frag_files) {
  tbi <- paste0(f, ".tbi")
  if (!file.exists(tbi)) {
    cat(sprintf("Indexing %s ...\n", basename(f)))
    system2("tabix", args = c("-p", "bed", f))
  }
}

# --- 4. 結合fragments → 1つのFragmentオブジェクト ---
# 全fragmentsを結合した一時ファイルを作成
cat("Merging fragment files ...\n")
merged_frag <- file.path(fragments_dir, "merged_fragments.tsv.gz")
if (!file.exists(merged_frag)) {
  # サブサンプル細胞のバーコードセット
  cell_set <- meta$cell.id

  # awkで該当細胞のみフィルタしながら結合
  cell_file <- tempfile(fileext = ".txt")
  writeLines(cell_set, cell_file)

  cmd <- sprintf(
    "zcat %s | awk -F'\\t' 'NR==FNR{a[$1];next} $4 in a' %s - | sort -k1,1 -k2,2n | bgzip > %s",
    paste(frag_files, collapse = " "),
    cell_file,
    merged_frag
  )
  cat("Running fragment merge + filter + sort ...\n")
  system(cmd)
  system2("tabix", args = c("-p", "bed", merged_frag))
  unlink(cell_file)
}
cat("Merged fragments ready.\n")

# --- 5. ピークコーリング (MACS3) ---
cat("Calling peaks with MACS3 ...\n")
peaks_dir <- file.path(fragments_dir, "macs3_output")
peaks_file <- file.path(peaks_dir, "atac_peaks_peaks.narrowPeak")
if (!file.exists(peaks_file)) {
  dir.create(peaks_dir, showWarnings = FALSE, recursive = TRUE)
  macs3_cmd <- sprintf(
    "macs3 callpeak -t %s -f BED --nomodel --shift -100 --extsize 200 --keep-dup all -g mm -n atac_peaks --outdir %s",
    merged_frag, peaks_dir
  )
  cat(sprintf("Running: %s\n", macs3_cmd))
  system(macs3_cmd)
}

# narrowPeak → GRanges
np <- read.table(peaks_file, stringsAsFactors = FALSE)
peaks <- GRanges(
  seqnames = np$V1,
  ranges   = IRanges(start = np$V2, end = np$V3)
)
cat(sprintf("Peaks called: %d\n", length(peaks)))

frags <- CreateFragmentObject(
  path = merged_frag,
  cells = meta$cell.id
)

# --- 6. ピーク×細胞行列 ---
cat("Creating peak-cell matrix ...\n")
peak_matrix <- FeatureMatrix(
  fragments = frags,
  features  = peaks,
  cells     = meta$cell.id
)
cat(sprintf("Peak matrix: %d x %d\n", nrow(peak_matrix), ncol(peak_matrix)))

# QCフィルタ
n_peaks_per_cell <- colSums(peak_matrix > 0)
keep <- n_peaks_per_cell >= min_peaks
peak_matrix <- peak_matrix[, keep]
meta <- meta[colnames(peak_matrix), ]
cat(sprintf("After min_peaks filter: %d cells\n", ncol(peak_matrix)))

# --- 7. Seuratオブジェクト ---
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- "UCSC"

chrom_assay <- CreateChromatinAssay(
  counts     = peak_matrix,
  fragments  = frags,
  genome     = "mm10",
  annotation = annotations
)

atac <- CreateSeuratObject(counts = chrom_assay, assay = "ATAC")
atac <- AddMetaData(atac, meta)

# --- 8. Gene Activity Matrix ---
cat("Computing Gene Activity Matrix ...\n")
gene_activities <- GeneActivity(atac)
atac[["GAM"]] <- CreateAssayObject(counts = gene_activities)

DefaultAssay(atac) <- "GAM"
atac <- NormalizeData(atac)
atac <- FindVariableFeatures(atac)
atac <- ScaleData(atac)

cat(sprintf("Final ATAC: %d cells, %d GAM genes\n",
  ncol(atac), nrow(atac[["GAM"]])))

# --- 9. 保存 ---
saveRDS(atac, out_rds)

meta_out <- atac@meta.data
meta_out$cell_id <- rownames(meta_out)
write.csv(meta_out, out_meta_csv, row.names = FALSE)

cat("Done.\n")
