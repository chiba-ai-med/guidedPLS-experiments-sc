# =============================================================================
# PBMC scMultiome 前処理 (10x 10k)
# =============================================================================
# 同一細胞由来の RNA + ATAC を 1 つの h5 から読み、Seurat + Signac で前処理。
# RNA 側でクラスタリング → marker-gene scoring で celltype と broad_lineage を付与。
# =============================================================================

PBMC_DATA_DIR = "data/PBMC"
PBMC_H5       = f"{PBMC_DATA_DIR}/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5"
PBMC_FRAG     = f"{PBMC_DATA_DIR}/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"

rule preprocess_pbmc:
    input:
        h5   = PBMC_H5,
        frag = PBMC_FRAG,
    output:
        rna       = "output/PBMC/preprocess/rna_combined.rds",
        atac      = "output/PBMC/preprocess/atac_processed.rds",
        rna_meta  = "output/PBMC/preprocess/rna_metadata.csv",
        atac_meta = "output/PBMC/preprocess/atac_metadata.csv",
    params:
        subsample = config.get("subsample_n", 10000),
    threads: 4
    shell:
        """
        mkdir -p output/PBMC/preprocess
        Rscript src/preprocess_pbmc.R \
            {input.h5} {input.frag} \
            {output.rna} {output.atac} \
            {output.rna_meta} {output.atac_meta} \
            {params.subsample}
        """
