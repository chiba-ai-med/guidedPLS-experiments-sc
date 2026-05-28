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
        emb="output/{dataset}/comparison/seurat/embeddings.csv",
    benchmark: "output/{dataset}/benchmark/seurat.tsv"
    params:
        seed=config["random_seed"],
    threads: 4
    shell:
        """
        mkdir -p output/{wildcards.dataset}/comparison/seurat
        Rscript src/run_seurat_transfer.R \
            {input.rna} {input.atac} \
            {output.labels} {params.seed} {output.emb}
        """

rule run_harmony:
    input:
        rna="output/{dataset}/preprocess/rna_combined.rds",
        atac="output/{dataset}/preprocess/atac_processed.rds",
    output:
        labels="output/{dataset}/comparison/harmony/predicted_labels.csv",
        emb="output/{dataset}/comparison/harmony/embeddings.csv",
    benchmark: "output/{dataset}/benchmark/harmony.tsv"
    params:
        seed=config["random_seed"],
        k=config["guidedpls_k"],
    threads: 4
    shell:
        """
        mkdir -p output/{wildcards.dataset}/comparison/harmony
        Rscript src/run_harmony_knn.R \
            {input.rna} {input.atac} \
            {output.labels} {params.seed} {params.k} {output.emb}
        """

rule run_scanorama:
    input:
        rna="output/{dataset}/preprocess/rna_combined.rds",
        atac="output/{dataset}/preprocess/atac_processed.rds",
    output:
        labels="output/{dataset}/comparison/scanorama/predicted_labels.csv",
        emb="output/{dataset}/comparison/scanorama/embeddings.csv",
    benchmark: "output/{dataset}/benchmark/scanorama.tsv"
    params:
        seed=config["random_seed"],
        k=config["guidedpls_k"],
    threads: 4
    shell:
        """
        mkdir -p output/{wildcards.dataset}/comparison/scanorama
        python3 src/run_scanorama.py \
            {input.rna} {input.atac} \
            {output.labels} {params.seed} {params.k} {output.emb}
        """
