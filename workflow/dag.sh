#!/bin/bash
# DAG graph (rulegraph) for the unified guidedPLS scRNA/scATAC workflow
mkdir -p plot
snakemake -s workflow/Snakefile --rulegraph | dot -Tpng > plot/Snakefile.png
