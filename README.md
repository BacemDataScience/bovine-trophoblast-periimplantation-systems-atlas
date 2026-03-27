# bovine-trophoblast-periimplantation-systems-atlas

This repository contains the analysis code used for the manuscript:

**Integrative single-cell systems analysis reveals coordinated IFNT signaling and extracellular matrix remodeling during bovine peri-implantation**

## Data source

The workflow analyzes the publicly available bovine peri-implantation single-cell RNA-seq dataset **GSE234335**.

Raw input files should be placed in:

`data/raw/`

Expected input files include:

- `*_matrix.mtx.gz`
- `*_features.tsv.gz`
- `*_barcodes.tsv.gz`

## Main analysis file

The full workflow is contained in:

`bovine_periimplantation_pipeline.R`

This script performs:
- data loading
- quality control
- integration and clustering
- marker analysis
- stage composition analysis
- differential expression
- module scoring
- monotonic gene analysis
- pseudobulk analysis
- pseudotime inference
- export of manuscript figures and supplementary tables

## Output

Analysis outputs are written to:

`results/`

## Software

The analysis was performed in R using Seurat and related packages.

Session information is provided in:

`sessionInfo.txt`

## Reproducibility note

This repository provides the code used to generate the main analytical results, figures, and supplementary tables. Large intermediate objects are not included because of file size limitations.
