# =============================================================================
# 比較手法の実行 — 全データセット共通
# =============================================================================
# 使用する変数:
#   COMPARISON_METHODS — ["seurat", "harmony", "scanorama"]
# =============================================================================

rule run_seurat:
    input:
        rna="output/{dataset}/preprocess/rna_combined.rds",
        atac="output/{dataset}/preprocess/atac_processed.rds",
    output:
        labels="output/{dataset}/comparison/seurat/predicted_labels.csv",
    params:
        seed=config["random_seed"],
    threads: 4
    shell:
        """
        mkdir -p output/{wildcards.dataset}/comparison/seurat
        Rscript src/run_seurat_transfer.R \
            {input.rna} {input.atac} \
            {output.labels} {params.seed}
        """

rule run_harmony:
    input:
        rna="output/{dataset}/preprocess/rna_combined.rds",
        atac="output/{dataset}/preprocess/atac_processed.rds",
    output:
        labels="output/{dataset}/comparison/harmony/predicted_labels.csv",
    params:
        seed=config["random_seed"],
        k=config["guidedpls_k"],
    threads: 4
    shell:
        """
        mkdir -p output/{wildcards.dataset}/comparison/harmony
        Rscript src/run_harmony_knn.R \
            {input.rna} {input.atac} \
            {output.labels} {params.seed} {params.k}
        """

rule run_scanorama:
    input:
        rna="output/{dataset}/preprocess/rna_combined.rds",
        atac="output/{dataset}/preprocess/atac_processed.rds",
    output:
        labels="output/{dataset}/comparison/scanorama/predicted_labels.csv",
    params:
        seed=config["random_seed"],
        k=config["guidedpls_k"],
    threads: 4
    shell:
        """
        mkdir -p output/{wildcards.dataset}/comparison/scanorama
        python3 src/run_scanorama.py \
            {input.rna} {input.atac} \
            {output.labels} {params.seed} {params.k}
        """
