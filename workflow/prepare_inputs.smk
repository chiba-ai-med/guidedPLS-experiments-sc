# =============================================================================
# guided-PLS 入力データ準備 (X1, X2, Y1, Y2) — 全データセット共通
# =============================================================================

rule prepare_inputs:
    input:
        rna_rds="output/{dataset}/preprocess/rna_combined.rds",
        rna_meta="output/{dataset}/preprocess/rna_metadata.csv",
        atac_rds="output/{dataset}/preprocess/atac_processed.rds",
        atac_meta="output/{dataset}/preprocess/atac_metadata.csv",
    output:
        x1="output/{dataset}/guidedpls_input/{condition}/X1.mtx",
        x2="output/{dataset}/guidedpls_input/{condition}/X2.mtx",
        y1="output/{dataset}/guidedpls_input/{condition}/Y1.csv",
        y2="output/{dataset}/guidedpls_input/{condition}/Y2.csv",
        rna_meta_out="output/{dataset}/guidedpls_input/{condition}/rna_metadata.csv",
        atac_meta_out="output/{dataset}/guidedpls_input/{condition}/atac_metadata.csv",
    params:
        seed=config["random_seed"],
    shell:
        """
        mkdir -p output/{wildcards.dataset}/guidedpls_input/{wildcards.condition}
        python3 src/prepare_guide_labels.py \
            --rna-rds {input.rna_rds} \
            --rna-meta {input.rna_meta} \
            --atac-rds {input.atac_rds} \
            --atac-meta {input.atac_meta} \
            --condition {wildcards.condition} \
            --dataset {wildcards.dataset} \
            --outdir output/{wildcards.dataset}/guidedpls_input/{wildcards.condition} \
            --seed {params.seed}
        """
