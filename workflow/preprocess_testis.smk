# =============================================================================
# Testis データセット前処理
# =============================================================================
# data/Testis/
#   Save-ArchR-Project.rds        — scATAC-Seq (ArchR project, peakSet + metadata)
#   testis.allcell.scRNAseq.rds   — scRNA-Seq (Seurat object)
#   E18/, P0/, P3/, P6/           — fragment files
# =============================================================================

rule preprocess_testis_rna:
    input:
        "data/Testis/testis.allcell.scRNAseq.rds",
    output:
        rds="output/Testis/preprocess/rna_combined.rds",
        meta="output/Testis/preprocess/rna_metadata.csv",
    params:
        subsample_n=config["subsample_n"],
        n_hvg=config["n_hvg"],
        seed=config["random_seed"],
    threads: 4
    shell:
        """
        mkdir -p output/Testis/preprocess
        Rscript src/preprocess_testis_rna.R \
            {input} {output.rds} {output.meta} \
            {params.subsample_n} {params.n_hvg} {params.seed}
        """

rule preprocess_testis_atac:
    input:
        archr="data/Testis/Save-ArchR-Project.rds",
        rna="output/Testis/preprocess/rna_combined.rds",
        e18="data/Testis/E18/alignments.possorted.tagged.bap.fragments.tsv.gz",
        p0="data/Testis/P0/alignments.possorted.tagged.bap.fragments.tsv.gz",
        p3="data/Testis/P3/alignments.possorted.tagged.bap.fragments.tsv.gz",
        p6="data/Testis/P6/alignments.possorted.tagged.bap.fragments.tsv.gz",
    output:
        rds="output/Testis/preprocess/atac_processed.rds",
        meta="output/Testis/preprocess/atac_metadata.csv",
    params:
        fragments_dir="data/Testis",
        min_peaks=config["min_peaks_atac"],
        subsample_n=config["subsample_n"],
        seed=config["random_seed"],
    threads: 8
    shell:
        """
        mkdir -p output/Testis/preprocess
        Rscript src/preprocess_testis_atac.R \
            {input.archr} {input.rna} \
            {params.fragments_dir} \
            {output.rds} {output.meta} \
            {params.min_peaks} {params.subsample_n} {params.seed}
        """
