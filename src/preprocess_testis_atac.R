# =============================================================================
# Testis scATAC-Seq 前処理: fragments + ArchR peakSet → ピーク行列 + GAM
# =============================================================================
# Usage: Rscript src/preprocess_testis_atac.R \
#   <archr.rds> <rna_combined.rds> \
#   <fragments_dir> \
#   <out_rds> <out_meta_csv> \
#   <min_peaks> <subsample_n> <seed>
# =============================================================================

suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
  library(GenomicRanges)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
archr_rds      <- args[1]
rna_rds        <- args[2]
fragments_dir  <- args[3]
out_rds        <- args[4]
out_meta_csv   <- args[5]
min_peaks      <- as.integer(args[6])
subsample_n    <- as.integer(args[7])
seed           <- as.integer(args[8])

set.seed(seed)
cat("=== Testis scATAC-Seq preprocessing ===\n")

# --- 1. ArchR メタデータ + peakSet 抽出 ---
archr <- readRDS(archr_rds)
meta <- as.data.frame(archr@cellColData)
meta$cell_id <- rownames(meta)

# ステージビニング
stage_map <- c(
  "E18.2" = "E18",
  "P0.2"  = "P2-3",
  "P3.2"  = "P2-3",
  "P6.2"  = "P6-7"
)
# Sample がNAの場合、cell_idからサンプル名を復元
# cell_id = "P3.2#alignments.possorted.tagged_BC00028_N01.1"
na_sample <- is.na(meta$Sample)
if (any(na_sample)) {
  meta$Sample[na_sample] <- sub("#.*", "", meta$cell_id[na_sample])
  cat(sprintf("Recovered %d NA samples from cell_id\n", sum(na_sample)))
}

meta$Stage <- stage_map[meta$Sample]
meta$celltype <- meta$predictedGroup_Co

cat(sprintf("ArchR cells: %d\n", nrow(meta)))
cat("Stage distribution:\n")
print(table(meta$Stage))

# peakSet (222,919 peaks)
peaks <- archr@peakSet
cat(sprintf("peakSet: %d peaks\n", length(peaks)))

# --- 2. サブサンプリング（先にメタデータレベルで） ---
if (nrow(meta) > subsample_n) {
  cat(sprintf("Subsampling to %d cells ...\n", subsample_n))
  cells_keep <- c()
  stage_table <- table(meta$Stage)
  for (stg in names(stage_table)) {
    stg_cells <- meta$cell_id[meta$Stage == stg]
    n_take <- round(as.numeric(subsample_n) * as.numeric(stage_table[stg]) / as.numeric(sum(stage_table)))
    n_take <- min(n_take, length(stg_cells))
    cells_keep <- c(cells_keep, sample(stg_cells, n_take))
  }
  meta <- meta[meta$cell_id %in% cells_keep, ]
}
cat(sprintf("Working with %d cells\n", nrow(meta)))

# --- 3. Fragment files ---
frag_files <- list.files(fragments_dir,
  pattern = "\\.fragments\\.tsv\\.gz$",
  full.names = TRUE, recursive = TRUE)
cat(sprintf("Found %d fragment files\n", length(frag_files)))

# --- 4. バーコード形式の確認・マッピング ---
# ArchRのバーコード: "P3.2#alignments.possorted.tagged_BC04940_N03"
# Fragmentのバーコード: "alignments.possorted.tagged_BC04940_N03"
# → サンプルプレフィックスを除去してマッチさせる必要がある

# ArchRバーコードからサンプル名とfragmentバーコードを分離
barcode_parts <- strsplit(meta$cell_id, "#")
meta$sample_prefix <- sapply(barcode_parts, `[`, 1)
meta$frag_barcode <- sapply(barcode_parts, `[`, 2)

cat("Sample prefix distribution:\n")
print(table(meta$sample_prefix))
cat("Fragment barcode examples:", head(meta$frag_barcode, 3), "\n")

# サンプル→ファイルのマッピング
sample_to_file <- c()
for (f in frag_files) {
  # ファイル名からサンプル識別
  dirname_f <- basename(dirname(f))  # E18, P0, P3, P6
  sample_to_file[dirname_f] <- f
}
cat("Sample to file mapping:\n")
print(sample_to_file)

# --- 5. Fragments結合 + ソート ---
merged_frag <- file.path(fragments_dir, "merged_sorted_fragments.tsv.gz")
if (!file.exists(merged_frag)) {
  cat("Merging and sorting fragment files ...\n")

  # サブサンプル細胞のバーコードセット（fragmentファイル内の形式）
  cell_file <- tempfile(fileext = ".txt")
  writeLines(meta$frag_barcode, cell_file)

  # 全fragmentsを結合、該当細胞をフィルタ、ソート
  sort_tmp <- file.path(fragments_dir, "sort_tmp")
  dir.create(sort_tmp, showWarnings = FALSE)
  cmd <- sprintf(
    "zcat %s | awk -F'\\t' 'NR==FNR{a[$1];next} $4 in a' %s - | sort -T %s -k1,1 -k2,2n | bgzip > %s",
    paste(frag_files, collapse = " "),
    cell_file,
    sort_tmp,
    merged_frag
  )
  system(cmd)
  system2("tabix", args = c("-p", "bed", merged_frag))
  unlink(cell_file)
}
cat("Merged fragments ready.\n")

# fragment内のバーコードを確認
frag_head <- read.delim(
  pipe(sprintf("zcat %s | head -5", merged_frag)),
  header = FALSE, stringsAsFactors = FALSE)
cat("Fragment columns:", ncol(frag_head), "\n")
cat("Fragment barcode example:", frag_head$V4[1], "\n")

# --- 6. FeatureMatrix (peakSet × cells) ---
cat("Creating Fragment object ...\n")
frags <- CreateFragmentObject(
  path = merged_frag,
  cells = meta$frag_barcode
)

cat("Creating peak-cell matrix ...\n")
peak_matrix <- FeatureMatrix(
  fragments = frags,
  features  = peaks,
  cells     = meta$frag_barcode
)
cat(sprintf("Peak matrix: %d x %d\n", nrow(peak_matrix), ncol(peak_matrix)))

# QCフィルタ
n_peaks_per_cell <- colSums(peak_matrix > 0)
keep <- n_peaks_per_cell >= min_peaks
peak_matrix <- peak_matrix[, keep]
cat(sprintf("After min_peaks filter: %d cells\n", ncol(peak_matrix)))

# コロン名をArchRスタイルに戻す（メタデータと対応させるため）
# fragment barcode → ArchR cell_id のマッピング
frag_to_archr <- setNames(meta$cell_id, meta$frag_barcode)
new_colnames <- frag_to_archr[colnames(peak_matrix)]
# マッピングできなかったものは除外
keep_mapped <- !is.na(new_colnames)
peak_matrix <- peak_matrix[, keep_mapped]
colnames(peak_matrix) <- new_colnames[keep_mapped]
meta <- meta[meta$cell_id %in% colnames(peak_matrix), ]
rownames(meta) <- meta$cell_id

cat(sprintf("After barcode mapping: %d cells\n", ncol(peak_matrix)))

# --- 7. Seuratオブジェクト ---
chrom_assay <- CreateChromatinAssay(
  counts     = peak_matrix,
  fragments  = frags,
  genome     = "mm10"
)

atac <- CreateSeuratObject(counts = chrom_assay, assay = "ATAC")
atac <- AddMetaData(atac, meta)

# --- 8. Gene Activity Matrix (比較手法用) ---
cat("Computing Gene Activity Matrix ...\n")
tryCatch({
  library(EnsDb.Mmusculus.v79)
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
  seqlevelsStyle(annotations) <- "UCSC"
  Annotation(atac) <- annotations

  gene_activities <- GeneActivity(atac)
  atac[["GAM"]] <- CreateAssayObject(counts = gene_activities)
  DefaultAssay(atac) <- "GAM"
  atac <- NormalizeData(atac)
  atac <- FindVariableFeatures(atac)
  atac <- ScaleData(atac)

  cat(sprintf("Final ATAC: %d cells, %d peaks, %d GAM genes\n",
    ncol(atac), nrow(atac[["ATAC"]]), nrow(atac[["GAM"]])))
}, error = function(e) {
  cat(sprintf("GAM computation failed: %s\n", e$message))
  cat("Saving without GAM (ATAC peaks only).\n")
  cat(sprintf("Final ATAC: %d cells, %d peaks\n",
    ncol(atac), nrow(atac[["ATAC"]])))
})

# --- 9. 保存 ---
saveRDS(atac, out_rds)

meta_out <- atac@meta.data
meta_out$cell_id <- rownames(meta_out)
write.csv(meta_out, out_meta_csv, row.names = FALSE)

cat("Done.\n")
