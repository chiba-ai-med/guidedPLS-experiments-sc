# =============================================================================
# Organogenesis データセット前処理
# =============================================================================
# data/Organogenesis/
#   CNP0003941_scATAC/  — scATAC-Seq fragments (E8.5/E9.5/E10.5)
#   GSE186069_scRNA/    — scRNA-Seq E8.5 (mtx)
#   GSE186068_scRNA/    — scRNA-Seq E9.5-E13.5 (loom)
# =============================================================================

rule preprocess_organogenesis_rna_e85:
    input:
        "data/Organogenesis/GSE186069_scRNA/GSE186069_gene_count.mtx.gz",
    output:
        rds="output/Organogenesis/preprocess/rna_GSE186069.rds",
    params:
        mtx="data/Organogenesis/GSE186069_scRNA/GSE186069_gene_count.mtx.gz",
        cells="data/Organogenesis/GSE186069_scRNA/GSE186069_cell_annotate.csv.gz",
        genes="data/Organogenesis/GSE186069_scRNA/GSE186069_gene_annotate.csv.gz",
        min_genes=config["min_genes_rna"],
        seed=config["random_seed"],
    threads: 4
    shell:
        """
        mkdir -p output/Organogenesis/preprocess
        Rscript src/preprocess_rna_e85.R \
            {params.mtx} {params.cells} {params.genes} \
            {output.rds} {params.min_genes} {params.seed}
        """

rule preprocess_organogenesis_rna_e95e105:
    input:
        "data/Organogenesis/GSE186068_scRNA/GSE186068_gene_count.loom.gz",
    output:
        rds="output/Organogenesis/preprocess/rna_GSE186068.rds",
    params:
        loom="data/Organogenesis/GSE186068_scRNA/GSE186068_gene_count.loom.gz",
        cells="data/Organogenesis/GSE186068_scRNA/GSE186068_cell_annotate.csv.gz",
        genes="data/Organogenesis/GSE186068_scRNA/GSE186068_gene_annotate.csv.gz",
        min_genes=config["min_genes_rna"],
        seed=config["random_seed"],
    threads: 4
    shell:
        """
        mkdir -p output/Organogenesis/preprocess
        python3 src/preprocess_rna_e95e105.py \
            {params.loom} {params.cells} {params.genes} \
            {output.rds} {params.min_genes} {params.seed}
        """

rule preprocess_organogenesis_merge_rna:
    input:
        rds1="output/Organogenesis/preprocess/rna_GSE186069.rds",
        rds2="output/Organogenesis/preprocess/rna_GSE186068.rds",
    output:
        rds="output/Organogenesis/preprocess/rna_combined.rds",
        meta="output/Organogenesis/preprocess/rna_metadata.csv",
    params:
        subsample_n=config["subsample_n"],
        n_hvg=config["n_hvg"],
        seed=config["random_seed"],
    threads: 4
    shell:
        """
        Rscript src/merge_rna.R \
            {input.rds1} {input.rds2} \
            {output.rds} {output.meta} \
            {params.subsample_n} {params.n_hvg} {params.seed}
        """

rule preprocess_organogenesis_atac:
    input:
        "data/Organogenesis/CNP0003941_scATAC/mouse.embryo.E8.5.E9.5-E10.5.metadata.upload.txt",
    output:
        rds="output/Organogenesis/preprocess/atac_processed.rds",
        meta="output/Organogenesis/preprocess/atac_metadata.csv",
    params:
        metadata="data/Organogenesis/CNP0003941_scATAC/mouse.embryo.E8.5.E9.5-E10.5.metadata.upload.txt",
        fragments_dir="data/Organogenesis/CNP0003941_scATAC",
        min_peaks=config["min_peaks_atac"],
        subsample_n=config["subsample_n"],
        seed=config["random_seed"],
    threads: 8
    shell:
        """
        mkdir -p output/Organogenesis/preprocess
        Rscript src/preprocess_atac.R \
            {params.metadata} {params.fragments_dir} \
            {output.rds} {output.meta} \
            {params.min_peaks} {params.subsample_n} {params.seed}
        """
