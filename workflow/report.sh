#!/bin/bash
# HTML report for the unified guidedPLS scRNA/scATAC workflow
mkdir -p report
snakemake -s workflow/Snakefile --report report/Snakefile.html
