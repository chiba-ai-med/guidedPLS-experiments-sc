# =============================================================================
# PBMC scMultiome (10x 10k) 前処理
# =============================================================================
# 同一細胞由来の RNA + ATAC を Seurat + Signac でロードし、
# RNA 側でクラスタリング → marker-gene scoring で celltype と broad_lineage を
# 付与、ATAC は TF-IDF + GAM を作成して保存する。
# =============================================================================
# Usage: Rscript src/preprocess_pbmc.R <h5> <fragments_tsv_gz> \
#   <out_rna.rds> <out_atac.rds> <out_rna_meta.csv> <out_atac_meta.csv> \
#   <subsample_n>
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(EnsDb.Hsapiens.v86)
  library(GenomeInfoDb)
})

args <- commandArgs(trailingOnly = TRUE)
h5_file       <- args[1]
frag_file     <- args[2]
out_rna_rds   <- args[3]
out_atac_rds  <- args[4]
out_rna_meta  <- args[5]
out_atac_meta <- args[6]
subsample_n   <- as.integer(args[7])

set.seed(42)
cat("=== PBMC scMultiome preprocessing ===\n")

# --- 1. h5 読み込み ---
counts <- Read10X_h5(h5_file)
rna_counts  <- counts[["Gene Expression"]]
atac_counts <- counts[["Peaks"]]
common_bc <- intersect(colnames(rna_counts), colnames(atac_counts))
cat(sprintf("Cells with paired RNA+ATAC: %d\n", length(common_bc)))

# --- 2. Seurat object (RNA) ---
pbmc <- CreateSeuratObject(counts = rna_counts[, common_bc], assay = "RNA")
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
cat(sprintf("Before QC: %d cells\n", ncol(pbmc)))

pbmc <- subset(pbmc,
  subset = nCount_RNA > 1000 & nCount_RNA < 25000 &
           nFeature_RNA > 500  & percent.mt < 20)
cat(sprintf("After QC: %d cells\n", ncol(pbmc)))

# --- 3. ATAC assay ---
annotation <- tryCatch({
  a <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  seqlevelsStyle(a) <- "UCSC"
  a
}, error = function(e) NULL)

atac_keep <- atac_counts[, colnames(pbmc)]
# Use peaks on canonical chromosomes only
chrs <- grepl("^chr[0-9XY]+", rownames(atac_keep))
atac_keep <- atac_keep[chrs, ]
cat(sprintf("Peaks (autosome+XY): %d\n", nrow(atac_keep)))

chrom_assay <- CreateChromatinAssay(
  counts     = atac_keep,
  sep        = c(":", "-"),
  fragments  = frag_file,
  annotation = annotation,
  min.cells  = 10
)
pbmc[["ATAC"]] <- chrom_assay
cat("ATAC assay attached.\n")

# --- 4. RNA 標準処理 ---
DefaultAssay(pbmc) <- "RNA"
pbmc <- NormalizeData(pbmc, verbose = FALSE)
pbmc <- FindVariableFeatures(pbmc, nfeatures = 3000, verbose = FALSE)
pbmc <- ScaleData(pbmc, verbose = FALSE)
pbmc <- RunPCA(pbmc, npcs = 30, verbose = FALSE)
pbmc <- FindNeighbors(pbmc, reduction = "pca", dims = 1:30, verbose = FALSE)
pbmc <- FindClusters(pbmc, resolution = 0.5, verbose = FALSE)
pbmc <- RunUMAP(pbmc, reduction = "pca", dims = 1:30, verbose = FALSE)
cat(sprintf("RNA clusters: %d\n", length(unique(Idents(pbmc)))))

# --- 5. marker-gene scoring + celltype 割り当て ---
# Canonical PBMC marker panels (Seurat tutorial / Azimuth level2 ベース)
markers <- list(
  `CD4 T Naive`   = c("CD3D","CD3E","CD4","CCR7","SELL","LEF1","TCF7"),
  `CD4 T Memory`  = c("CD3D","CD3E","CD4","IL7R","S100A4","KLRB1","ANXA1"),
  `CD8 T Naive`   = c("CD3D","CD3E","CD8A","CD8B","CCR7","LEF1"),
  `CD8 T Memory`  = c("CD3D","CD3E","CD8A","CD8B","GZMK","CCL5","NKG7"),
  `NK`            = c("GNLY","NKG7","KLRD1","KLRF1","NCAM1","FCGR3A"),
  `B Naive`       = c("MS4A1","CD79A","CD79B","IGHM","IGHD","TCL1A"),
  `B Memory`      = c("MS4A1","CD79A","CD27","TNFRSF13B","AIM2"),
  `CD14 Mono`     = c("CD14","LYZ","S100A8","S100A9","CSTA","VCAN"),
  `CD16 Mono`     = c("FCGR3A","MS4A7","CDKN1C","LST1","FCER1G"),
  `cDC`           = c("FCER1A","CST3","CLEC10A","CD1C","HLA-DQA1"),
  `pDC`           = c("LILRA4","IL3RA","CLEC4C","JCHAIN"),
  `Platelet`      = c("PPBP","PF4","NRGN","GP9","TUBB1")
)

for (ct in names(markers)) {
  gs <- intersect(markers[[ct]], rownames(pbmc))
  if (length(gs) == 0) next
  pbmc <- AddModuleScore(pbmc, features = list(gs),
    name = paste0("score_", gsub(" |/|\\+", "", ct)), nbin = 12)
}
score_cols <- grep("^score_", colnames(pbmc@meta.data), value = TRUE)

# Per cluster: which score is highest on average → celltype name
cluster_ids <- as.character(unique(Idents(pbmc)))
cluster_celltype <- sapply(cluster_ids, function(cl) {
  cells <- WhichCells(pbmc, idents = cl)
  m <- colMeans(pbmc@meta.data[cells, score_cols, drop = FALSE])
  best <- names(markers)[which.max(m)]
  best
})
names(cluster_celltype) <- cluster_ids
cat("\n=== Cluster → celltype assignment ===\n")
for (cl in sort(cluster_ids)) cat(sprintf("Cluster %s -> %s (n=%d)\n",
  cl, cluster_celltype[cl], sum(Idents(pbmc) == cl)))

# NOTE: cluster_celltype は names=cluster ID の named vector。indexed lookup
# の結果がそのまま入ると Seurat 内部 AddMetaData が「名前が cell ID と
# overlap しない」と判定して死ぬので unname() で剥がす必要がある。
pbmc$celltype <- unname(cluster_celltype[as.character(Idents(pbmc))])

# Broad lineage (guide Z 候補)
broad_map <- c(
  `CD4 T Naive`  = "T", `CD4 T Memory` = "T",
  `CD8 T Naive`  = "T", `CD8 T Memory` = "T",
  `NK`           = "NK",
  `B Naive`      = "B", `B Memory`     = "B",
  `CD14 Mono`    = "Myeloid", `CD16 Mono` = "Myeloid",
  `cDC`          = "Myeloid", `pDC` = "Myeloid",
  `Platelet`     = "Other"
)
pbmc$broad_lineage <- unname(broad_map[pbmc$celltype])

cat("\n=== celltype distribution ===\n"); print(table(pbmc$celltype))
cat("\n=== broad_lineage distribution ===\n"); print(table(pbmc$broad_lineage))

# --- 6. ATAC 処理: TF-IDF + 上位 features ---
DefaultAssay(pbmc) <- "ATAC"
pbmc <- RunTFIDF(pbmc, verbose = FALSE)
pbmc <- FindTopFeatures(pbmc, min.cutoff = "q5", verbose = FALSE)

# --- 7. GAM (Gene Activity Matrix) — baseline 用 ---
gam <- tryCatch({
  GeneActivity(pbmc, verbose = FALSE)
}, error = function(e) {
  cat("GeneActivity failed:", conditionMessage(e), "\n")
  NULL
})
if (!is.null(gam)) {
  pbmc[["GAM"]] <- CreateAssayObject(counts = gam)
  DefaultAssay(pbmc) <- "GAM"
  pbmc <- NormalizeData(pbmc, verbose = FALSE)
  cat(sprintf("GAM created: %d genes x %d cells\n", nrow(gam), ncol(gam)))
}

# --- 8. サブサンプル (任意; subsample_n=0 で skip) ---
if (subsample_n > 0 && subsample_n < ncol(pbmc)) {
  cat(sprintf("Subsampling to %d cells ...\n", subsample_n))
  set.seed(42)
  keep_cells <- sample(colnames(pbmc), subsample_n)
  pbmc <- subset(pbmc, cells = keep_cells)
}

# --- 9. RNA と ATAC を別ファイルとして保存 ---
# 同一細胞ペアなので、cell barcode に suffix を付けて両側で衝突しないようにする
# (Harmony の merge 等が duplicate cell name でコケるのを防ぐ)。
# 元 pbmc を deep-copy せず、必要な assay だけ残して名前変更する。

# RNA only object
rna_obj <- pbmc
DefaultAssay(rna_obj) <- "RNA"
rna_obj[["ATAC"]] <- NULL
if ("GAM" %in% Assays(rna_obj)) rna_obj[["GAM"]] <- NULL
rna_obj <- RenameCells(rna_obj, new.names = paste0(colnames(rna_obj), "_R"))
rna_obj <- NormalizeData(rna_obj, verbose = FALSE)
rna_obj <- FindVariableFeatures(rna_obj, nfeatures = 3000, verbose = FALSE)
saveRDS(rna_obj, out_rna_rds)
cat(sprintf("RNA saved: %s (%d cells x %d genes)\n",
            out_rna_rds, ncol(rna_obj), nrow(rna_obj)))

# ATAC only object
atac_obj <- pbmc
atac_obj[["RNA"]] <- NULL
atac_obj <- RenameCells(atac_obj, new.names = paste0(colnames(atac_obj), "_A"))
if ("GAM" %in% Assays(atac_obj)) {
  DefaultAssay(atac_obj) <- "GAM"
} else {
  DefaultAssay(atac_obj) <- "ATAC"
}
saveRDS(atac_obj, out_atac_rds)
cat(sprintf("ATAC saved: %s\n", out_atac_rds))

# --- 10. metadata 出力 (新しい cell ID で) ---
md_rna  <- rna_obj@meta.data
md_atac <- atac_obj@meta.data
md_rna$cell_id   <- rownames(md_rna)
md_atac$cell_id  <- rownames(md_atac)
# 元の paired barcode (PBMC ペア評価用に保持)
md_rna$paired_barcode  <- sub("_R$", "", rownames(md_rna))
md_atac$paired_barcode <- sub("_A$", "", rownames(md_atac))
# Snakemake / prepare_guide_labels が期待する列を埋める
md_rna$Stage <- "PBMC"
md_atac$Stage <- "PBMC"
md_rna$germ_layer    <- md_rna$broad_lineage
md_atac[["Germ.layer"]] <- md_atac$broad_lineage
write.csv(md_rna,  out_rna_meta,  row.names = TRUE)
write.csv(md_atac, out_atac_meta, row.names = TRUE)
cat(sprintf("Metadata saved: %s, %s\n", out_rna_meta, out_atac_meta))

cat("Done.\n")
