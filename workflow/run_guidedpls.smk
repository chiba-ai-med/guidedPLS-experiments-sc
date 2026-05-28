# =============================================================================
# guided-PLS 実行 — 全データセット共通
# =============================================================================

rule run_guidedpls:
    input:
        x1="output/{dataset}/guidedpls_input/{condition}/X1.csv",
        x2="output/{dataset}/guidedpls_input/{condition}/X2.csv",
        y1="output/{dataset}/guidedpls_input/{condition}/Y1.csv",
        y2="output/{dataset}/guidedpls_input/{condition}/Y2.csv",
        rna_meta="output/{dataset}/guidedpls_input/{condition}/rna_metadata.csv",
    output:
        labels="output/{dataset}/guidedpls/{condition}/predicted_labels.csv",
        rdata="output/{dataset}/guidedpls/{condition}/guidedpls.RData",
    benchmark: "output/{dataset}/benchmark/guidedpls_{condition}.tsv"
    params:
        r=config["guidedpls_r"],
        k=config["guidedpls_k"],
    threads: 4
    shell:
        """
        mkdir -p output/{wildcards.dataset}/guidedpls/{wildcards.condition}
        Rscript src/run_guidedpls.R \
            {input.x1} {input.x2} {input.y1} {input.y2} \
            {input.rna_meta} \
            {output.labels} {output.rdata} \
            {params.r} {params.k}
        """

rule save_gpls_embeddings:
    input:
        rdata="output/{dataset}/guidedpls/{condition}/guidedpls.RData",
        rna_meta="output/{dataset}/preprocess/rna_metadata.csv",
        atac_meta="output/{dataset}/preprocess/atac_metadata.csv",
    output:
        emb="output/{dataset}/guidedpls/{condition}/embeddings.csv",
    shell:
        """
        Rscript src/save_gpls_embeddings.R \
            {input.rdata} {input.rna_meta} {input.atac_meta} {output.emb}
        """
