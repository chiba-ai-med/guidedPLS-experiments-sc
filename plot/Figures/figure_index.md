# Figure Index — guidedPLS-experiments-sc

Single-cell multi-omics experiments contributing to **Fig. 3** of the guided-PLS manuscript.

**Important convention used throughout this repo:**

- **Guide Z** is one of `{none, stage, germlayer, stage_germlayer, stage_binned}` — defined in `config/config.yaml`.
- **Evaluation label** is `celltype` (cell type annotation from the original publications).
- **Cell type is NEVER used as a guide / training input for guided-PLS.** It is only used post-hoc to score the alignment quality. All baseline methods (Seurat, Harmony, Scanorama) likewise produce predictions that are scored against `celltype`.

## Main Figures (`plot/Figures/main/`)

| File | Panel | Dataset | Guide Z | Eval label | Method comparison | What it shows | Main/Supp | Notes |
|---|---|---|---|---|---|---|---|---|
| `Fig3A_task_overview.{png,pdf}` *(TBD — schematic)* | A | Both | — | celltype | — | Task schematic: unpaired scRNA-Seq / scATAC-Seq integration, role of guide Z, distinction between guide Z and evaluation label, Organogenesis vs Testis data structure | Main | Conceptual figure, not data-driven. To be drawn (Illustrator/Inkscape). |
| `Fig3B_organogenesis_performance.png` | B | Organogenesis (E8.5/E9.5/E10.5) | none / stage / germlayer / stage_germlayer | celltype | guidedPLS (4 guide conditions) vs Seurat / Harmony / Scanorama | Bar plot of accuracy and balanced accuracy across all methods. Currently Seurat leads (acc ≈ 0.48), best gPLS = stage_germlayer (acc ≈ 0.16). | Main | Source: `output/Organogenesis/figures/barplot_accuracy.png`. The `none` condition has predicted labels but is currently absent from `metrics.csv` — investigate. |
| `Fig3C_testis_performance.png` | C | Testis (E18/P0/P3/P6) | none / stage_binned | celltype | guidedPLS only (Seurat/Harmony/Scanorama not run; fragment-barcode format prevents GAM) | Bar plot of accuracy and balanced accuracy for the two gPLS guide conditions. `none` (acc ≈ 0.29) outperforms `stage_binned` (acc ≈ 0.22) — interpret with stage mismatch in mind. | Main | Source: `output/Testis/figures/barplot_accuracy.png`. Used as **robustness / secondary benchmark**, not primary claim. |
| `Fig3D_trueZ_vs_noZ_or_shuffledZ.{png,pdf}` *(TBD — needs experiment)* | D | Both | true Z / no Z / shuffled Z | celltype | guidedPLS only | Effect of guide Z informativeness: true Z vs no Z (already partially in B/C) vs shuffled Z (not yet implemented). | Main | Shuffled-Z baseline not yet run. Until then, the `none` vs guided comparison in Fig3B/C is the closest available signal. |
| `Fig3E_organogenesis_confusion_best_gPLS.png` | E | Organogenesis | stage_germlayer | celltype | guidedPLS (best gPLS variant) | Confusion matrix of best-performing gPLS condition. Rows = true cell type, columns = predicted cell type, normalized within row. | Main | Source: `output/Organogenesis/figures/confusion_guidedpls_stage_germlayer.png`. Alternative: per-class F1 heatmap (not yet generated; CSVs in supplementary). |

## Supplementary Figures and Tables (`plot/Figures/supplementary/`)

| File | Type | Dataset | Methods / Conditions | What it shows | Notes |
|---|---|---|---|---|---|
| `SuppFig_organogenesis_ari_nmi.png` | Bar plot | Organogenesis | All (gPLS × 4 + 3 baselines) | ARI and NMI per method | Complement to Fig3B |
| `SuppFig_testis_ari_nmi.png` | Bar plot | Testis | gPLS × 2 | ARI and NMI per method | Complement to Fig3C |
| `SuppFig_summary_both_datasets.png` | Cross-dataset | Both | All available | Side-by-side comparison of headline metrics across datasets | — |
| `SuppFig_organogenesis_confusion_guidedpls_none.png` | Confusion matrix | Organogenesis | guidedPLS, Z=none | Per-method confusion matrix | — |
| `SuppFig_organogenesis_confusion_guidedpls_stage.png` | Confusion matrix | Organogenesis | guidedPLS, Z=stage | — | — |
| `SuppFig_organogenesis_confusion_guidedpls_germlayer.png` | Confusion matrix | Organogenesis | guidedPLS, Z=germlayer | — | — |
| `SuppFig_organogenesis_confusion_guidedpls_stage_germlayer.png` | Confusion matrix | Organogenesis | guidedPLS, Z=stage_germlayer | Same image as Fig3E (kept here for completeness) | — |
| `SuppFig_organogenesis_confusion_seurat.png` | Confusion matrix | Organogenesis | Seurat | — | — |
| `SuppFig_organogenesis_confusion_harmony.png` | Confusion matrix | Organogenesis | Harmony | — | — |
| `SuppFig_organogenesis_confusion_scanorama.png` | Confusion matrix | Organogenesis | Scanorama | — | — |
| `SuppFig_testis_confusion_guidedpls_none.png` | Confusion matrix | Testis | guidedPLS, Z=none | — | — |
| `SuppFig_testis_confusion_guidedpls_stage_binned.png` | Confusion matrix | Testis | guidedPLS, Z=stage_binned | — | — |
| `SuppFig_workflow_dag.png` | DAG | — | — | Snakemake rule graph for the entire pipeline | Regenerable via `bash workflow/dag.sh` |
| `SuppTable_organogenesis_metrics.csv` | Table | Organogenesis | All | accuracy / balanced_accuracy / ARI / NMI / macro_f1 / weighted_f1 / n_cells per method | — |
| `SuppTable_testis_metrics.csv` | Table | Testis | gPLS × 2 | same columns | — |
| `SuppTable_organogenesis_per_class_f1.csv` | Table | Organogenesis | All | Per-cell-type F1 score for every method | Candidate source for a per-class F1 heatmap figure |
| `SuppTable_testis_per_class_f1.csv` | Table | Testis | gPLS × 2 | Per-cell-type F1 score | — |

## Known gaps

- **Fig3A schematic** has not been drawn yet.
- **Fig3D shuffled-Z baseline** has not been run.
- The `guidedpls_none` row is missing from `output/Organogenesis/evaluation/metrics.csv` even though predicted labels and confusion matrix exist — needs investigation before Fig3D can be assembled cleanly.
- **UMAP / latent score plots** are listed in the task brief but not yet produced by the workflow.

## Regenerating figures

```bash
# Full pipeline (regenerates output/{dataset}/figures/)
conda activate snakemake
snakemake -s workflow/Snakefile --cores 8 -p

# Then refresh the curated copies
cp output/Organogenesis/figures/barplot_accuracy.png plot/Figures/main/Fig3B_organogenesis_performance.png
# ... etc. (see this file as the authoritative mapping)
```
