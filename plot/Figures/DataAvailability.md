# Data Availability — guidedPLS-experiments-sc

Datasets used in this repository, with provenance for the manuscript Data Availability statement.

## Primary dataset (Fig 3 main): 10x Genomics 10k PBMC scMultiome

Paired single-cell RNA-Seq + ATAC-Seq from a healthy human donor, ~11k cells after QC.

- **Publisher**: 10x Genomics, Inc. (public demonstration dataset, Chromium X)
- **Catalog page**: https://www.10xgenomics.com/datasets/10k-Human-PBMCs-Multiome-v1-0-Chromium-X
- **Files used**:
  - `pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5` (RNA + ATAC counts)
  - `pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz` (+ `.tbi` index)
  - `pbmc_granulocyte_sorted_10k_per_barcode_metrics.csv`
- **Direct download URLs** (Cell Ranger ARC v2.0.0 release):
  - https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5
  - https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz
  - https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz.tbi
  - https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_per_barcode_metrics.csv
- **Licence**: 10x Genomics demonstration data; freely available for academic use under the 10x Genomics Software License Agreement.
- **Cell-type annotation in this repository**: derived in-repo by Seurat clustering of the RNA modality + marker-gene scoring (`src/preprocess_pbmc.R`); no external annotation file is required.

## Scale-up dataset (Supp): mouse organogenesis

Unpaired scRNA-Seq and scATAC-Seq of mouse embryonic development (E8.5–E10.5).

### scRNA-Seq (GEO)
- **GSE186069** — E8.5 (sci-RNA-seq3)
  - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186069
  - Files: `GSE186069_gene_count.mtx.gz`, `GSE186069_cell_annotate.csv.gz`, `GSE186069_gene_annotate.csv.gz`
- **GSE186068** — E9.5–E13.5 (sci-RNA-seq3)
  - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186068
  - Files: `GSE186068_gene_count.loom.gz`, `GSE186068_cell_annotate.csv.gz`, `GSE186068_gene_annotate.csv.gz`
- **Source paper**: Qiu, C. *et al.* "Systematic reconstruction of the cellular trajectories of mammalian embryogenesis." *Nature Genetics* 54, 328–341 (2022).
- This repository uses only E8.5 / E9.5 / E10.5 timepoints for parity with the ATAC side.

### scATAC-Seq (CNGB)
- **CNP0003941** (CNGB Nucleotide Sequence Archive)
  - FTP: https://ftp.cngb.org/pub/CNSA/data5/CNP0003941/Single_Cell/CSE0000217/
  - Files: per-stage `*.fragments.tsv.gz` (E8.5 × 9 libraries, E9.5 × 10, E10.5 × 10) + `mouse.embryo.E8.5.E9.5-E10.5.metadata.upload.txt`
- **Source paper**: "Single-cell chromatin accessibility profiling of cell-state-specific gene regulatory programs during mouse organogenesis." *Frontiers in Neuroscience* (2023).
- Cell-type and germ-layer annotations are taken directly from the published metadata file.

## Dropped from manuscript scope: mouse testis

This dataset is preserved in `data/Testis/` and `output/Testis/` for internal reference but is **not part of the manuscript**. It was excluded because the upstream ArchR scATAC representation uses a `sample#barcode` cell-ID format that is not directly compatible with the Signac Gene Activity Matrix construction the baseline integration methods require; the baseline comparison could therefore not be run fairly, and the explanation cost outweighed the narrative value (see commit `ca38b90`).

- **scRNA-Seq**: stages E18, P0, P3, P6.
- **scATAC-Seq**: ArchR project (`Save-ArchR-Project.rds`), stages E18, P0, P3, P6.
- **Origin**: not from a public archive. To the best of our knowledge this dataset reached the repository as a pre-processed ArchR project (RDS) rather than from a GEO / SRA / CNGB accession. **If kept in any future iteration of the manuscript, please confirm the original publication / collaborator attribution before listing it in Data Availability.** Treat the rows below as a placeholder until confirmed:
  - Original collaborator / lab: *to be filled in by the corresponding author*.
  - Public deposition (if any): *to be filled in*.
  - Re-use permissions: *to be confirmed*.

## Code availability

This Snakemake workflow (preprocessing, integration, baselines, evaluation, figure generation) is published in this repository under the MIT licence:

- https://github.com/chiba-ai-med/guidedPLS-experiments-sc

The `guidedPLS` R package itself is in a separate repository (see the manuscript's main Code Availability statement).
