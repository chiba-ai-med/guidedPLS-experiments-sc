#!/usr/bin/env python3
# =============================================================================
# Scanorama + kNN Label Transfer
# =============================================================================
# Usage: python3 src/run_scanorama.py \
#   <rna.rds> <atac.rds> <out_labels.csv> <seed> <k>
# =============================================================================

import sys
import subprocess
import numpy as np
import pandas as pd
import scanorama
from sklearn.neighbors import KNeighborsClassifier

rna_rds    = sys.argv[1]
atac_rds   = sys.argv[2]
out_labels = sys.argv[3]
seed       = int(sys.argv[4])
k_nn       = int(sys.argv[5])

np.random.seed(seed)
print("=== Scanorama + kNN Label Transfer ===")

# --- 1. RDSからCSVを抽出 (R経由) ---
r_script = f"""
suppressPackageStartupMessages(library(Seurat))

rna  <- readRDS("{rna_rds}")
atac <- readRDS("{atac_rds}")
DefaultAssay(atac) <- "GAM"

common_genes <- intersect(VariableFeatures(rna), rownames(atac))
if (length(common_genes) < 100) {{
    common_genes <- intersect(rownames(rna), rownames(atac))
    common_genes <- head(sort(common_genes), 3000)
}}
cat(sprintf("Common genes: %d\\n", length(common_genes)))

if (inherits(rna[["RNA"]], "Assay5")) rna[["RNA"]] <- JoinLayers(rna[["RNA"]])
x_rna  <- as.matrix(LayerData(rna, layer = "data")[common_genes, ])
if (inherits(atac[["GAM"]], "Assay5")) atac[["GAM"]] <- JoinLayers(atac[["GAM"]])
x_atac <- as.matrix(LayerData(atac, layer = "data")[common_genes, ])

write.csv(t(x_rna),  "/tmp/scanorama_rna.csv",  row.names = FALSE)
write.csv(t(x_atac), "/tmp/scanorama_atac.csv", row.names = FALSE)

meta_rna <- data.frame(celltype = rna$celltype, stringsAsFactors = FALSE)
write.csv(meta_rna, "/tmp/scanorama_rna_meta.csv", row.names = FALSE)
"""

print("Extracting data via R ...")
result = subprocess.run(["Rscript", "-e", r_script], capture_output=True, text=True)
if result.returncode != 0:
    print("R stderr:", result.stderr, file=sys.stderr)
    sys.exit(1)
print(result.stdout)

# --- 2. データ読み込み ---
print("Loading matrices ...")
rna_mat  = pd.read_csv("/tmp/scanorama_rna.csv").values   # cells × genes
atac_mat = pd.read_csv("/tmp/scanorama_atac.csv").values
rna_meta = pd.read_csv("/tmp/scanorama_rna_meta.csv")

gene_names = pd.read_csv("/tmp/scanorama_rna.csv", nrows=0).columns.tolist()

print(f"RNA: {rna_mat.shape}, ATAC: {atac_mat.shape}")

# --- 3. Scanorama ---
print("Running Scanorama ...")
corrected, genes = scanorama.correct(
    [rna_mat, atac_mat],
    [gene_names, gene_names],
    return_dimred=False
)

# 補正後行列から埋め込みを取得
integrated, genes_int = scanorama.integrate(
    [rna_mat, atac_mat],
    [gene_names, gene_names]
)

rna_embed  = integrated[0]
atac_embed = integrated[1]

print(f"Embedding: RNA {rna_embed.shape}, ATAC {atac_embed.shape}")

# --- 4. kNN Label Transfer ---
print(f"kNN label transfer (k={k_nn}) ...")
knn = KNeighborsClassifier(n_neighbors=k_nn)
knn.fit(rna_embed, rna_meta["celltype"].values)

predicted = knn.predict(atac_embed)
proba = knn.predict_proba(atac_embed)
scores = proba.max(axis=1)

# --- 5. 保存 ---
result_df = pd.DataFrame({
    "predicted_celltype": predicted,
    "prediction_score": scores,
})
result_df.to_csv(out_labels, index=False)

print(f"Predicted labels saved: {len(result_df)} cells")
print("Done.")
