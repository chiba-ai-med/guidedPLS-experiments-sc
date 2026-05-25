#!/usr/bin/env python3
# =============================================================================
# scRNA-Seq E9.5-E10.5 前処理 (GSE186068, Loom形式)
# メモリ効率重視: チャンクごとにスパース行列構築
# =============================================================================
# Usage: python3 src/preprocess_rna_e95e105.py \
#   <loom.gz> <cell_annotate.csv.gz> <gene_annotate.csv.gz> \
#   <out_rds> <min_genes> <seed>
# =============================================================================

import sys
import os
import gzip
import shutil
import subprocess
import gc
import numpy as np
import pandas as pd
import h5py
from scipy import sparse
from scipy.io import mmwrite

loom_gz    = sys.argv[1]
cells_file = sys.argv[2]
genes_file = sys.argv[3]
out_rds    = sys.argv[4]
min_genes  = int(sys.argv[5])
seed       = int(sys.argv[6])

np.random.seed(seed)
print("=== scRNA-Seq E9.5-E10.5 preprocessing (Python, memory-efficient) ===")

# --- 1. Loom解凍 ---
loom_file = loom_gz.replace(".gz", "")
if not os.path.exists(loom_file):
    print("Decompressing loom file ...")
    with gzip.open(loom_gz, 'rb') as f_in:
        with open(loom_file, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out, length=64*1024*1024)

# --- 2. メタデータ読み込み・フィルタ ---
print("Reading cell annotations ...")
cell_ann = pd.read_csv(cells_file)
gene_ann = pd.read_csv(genes_file)

cell_ann["Stage"] = "E" + cell_ann["development_stage"].astype(str)

keep_stage = cell_ann["Stage"].isin(["E9.5", "E10.5"])
keep_qc = cell_ann["removed_by_low_quality_or_doublets"] == "No"
keep_ct = cell_ann["celltype"].notna() & (cell_ann["celltype"] != "")
keep = (keep_stage & keep_qc & keep_ct).values
keep_idx = np.where(keep)[0]

cell_ann_filtered = cell_ann.loc[keep].reset_index(drop=True)
print(f"After stage + QC filter: {len(cell_ann_filtered)} cells")

# --- 3. Loomからスパース行列をチャンク読み込み ---
print("Reading Loom file (chunked, memory-efficient) ...")
with h5py.File(loom_file, 'r') as f:
    matrix_ds = f['matrix']
    dim0, dim1 = matrix_ds.shape
    print(f"Loom matrix shape: {dim0} x {dim1}")

    # loom形式: matrix は genes × cells or cells × genes
    # dim0=49585(=genes), dim1=2452396(=cells) → genes × cells
    if dim0 == len(gene_ann) and dim1 == len(cell_ann):
        is_gene_by_cell = True
        n_genes_total = dim0
        n_cells_total = dim1
        print("Matrix orientation: genes x cells")
    else:
        is_gene_by_cell = False
        n_cells_total = dim0
        n_genes_total = dim1
        print("Matrix orientation: cells x genes")

    sorted_idx = np.sort(keep_idx)
    n_keep = len(sorted_idx)

    # チャンクごとにスパース行列を構築
    chunk_size = 2000
    sparse_chunks = []

    for start in range(0, n_keep, chunk_size):
        end = min(start + chunk_size, n_keep)
        batch_idx = sorted_idx[start:end]

        if is_gene_by_cell:
            # genes × cells → 列(cells)方向でスライス
            chunk = matrix_ds[:, batch_idx]  # genes × chunk_cells
            sparse_chunks.append(sparse.csc_matrix(chunk))
        else:
            # cells × genes → 行(cells)方向でスライス
            if len(batch_idx) > 0 and batch_idx[-1] - batch_idx[0] + 1 == len(batch_idx):
                chunk = matrix_ds[batch_idx[0]:batch_idx[-1]+1, :]
            else:
                chunk = matrix_ds[batch_idx, :]
            sparse_chunks.append(sparse.csr_matrix(chunk))

        if start % (chunk_size * 20) == 0:
            print(f"  {start}/{n_keep} cells processed ...")
            gc.collect()

    print("Stacking sparse chunks ...")
    if is_gene_by_cell:
        # genes × cells: 列方向に結合
        mat = sparse.hstack(sparse_chunks, format='csc')
    else:
        mat = sparse.vstack(sparse_chunks, format='csr')
        mat = mat.T.tocsc()

    del sparse_chunks
    gc.collect()

print(f"Sparse matrix: {mat.shape[0]} genes x {mat.shape[1]} cells")

# --- 4. 最小遺伝子数フィルタ ---
genes_per_cell = np.array((mat > 0).sum(axis=0)).ravel()
keep_genes_mask = genes_per_cell >= min_genes
mat = mat[:, keep_genes_mask]
cell_ann_filtered = cell_ann_filtered.loc[keep_genes_mask].reset_index(drop=True)
print(f"After min_genes filter: {mat.shape[1]} cells")

# ゼロ発現遺伝子除去
gene_sums = np.array(mat.sum(axis=1)).ravel()
gene_keep = gene_sums > 0
mat = mat[gene_keep, :]

# gene_annotateのカラム名を確認
if "gene_short_name" in gene_ann.columns:
    gene_col = "gene_short_name"
elif "gene_name" in gene_ann.columns:
    gene_col = "gene_name"
else:
    gene_col = gene_ann.columns[0]

gene_names = gene_ann[gene_col].values[gene_keep]
# 重複遺伝子名をユニーク化
from collections import Counter
name_count = Counter()
unique_names = []
for name in gene_names:
    name_count[name] += 1
    if name_count[name] > 1:
        unique_names.append(f"{name}.{name_count[name]-1}")
    else:
        unique_names.append(name)
gene_names = unique_names

print(f"Non-zero genes: {mat.shape[0]}")

# --- 5. 中間ファイル保存 ---
tmp_dir = "/tmp/rna_e95e105_preprocess"
os.makedirs(tmp_dir, exist_ok=True)

print("Saving intermediate mtx ...")
mmwrite(f"{tmp_dir}/matrix.mtx", mat.tocoo())
del mat
gc.collect()

pd.DataFrame({"gene": gene_names}).to_csv(f"{tmp_dir}/genes.csv", index=False)
cell_ann_filtered.to_csv(f"{tmp_dir}/cells.csv", index=False)

# --- 6. Rで Seuratオブジェクト作成 ---
print("Creating Seurat object via R ...")
r_script = f"""
suppressPackageStartupMessages({{
  library(Seurat)
  library(Matrix)
}})
set.seed({seed})

cat("Reading mtx ...\\n")
mat <- readMM("{tmp_dir}/matrix.mtx")
genes <- read.csv("{tmp_dir}/genes.csv", stringsAsFactors=FALSE)
cells <- read.csv("{tmp_dir}/cells.csv", stringsAsFactors=FALSE)

rownames(mat) <- genes$gene
colnames(mat) <- cells$sample

cat("Creating Seurat object ...\\n")
rna <- CreateSeuratObject(counts = mat, project = "RNA_E95_E105")
rownames(cells) <- cells$sample
rna <- AddMetaData(rna, cells)
rna <- NormalizeData(rna)

cat(sprintf("Final RNA E9.5-E10.5: %d genes x %d cells\\n", nrow(rna), ncol(rna)))
saveRDS(rna, "{out_rds}")
cat("Done.\\n")
"""

result = subprocess.run(["Rscript", "-e", r_script], capture_output=True, text=True)
print(result.stdout)
if result.returncode != 0:
    print("R stderr:", result.stderr, file=sys.stderr)
    sys.exit(1)

shutil.rmtree(tmp_dir, ignore_errors=True)
print("All done.")
