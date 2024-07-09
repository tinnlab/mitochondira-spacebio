# Gene set enrichment analysis for mitochondrial pathways in space biology data

This folder contains R scripts to download and process gene expression data from the NASA GeneLab repository for gene set enrichment analysis (GSEA) using mitochondrial pathways.
The following packages are required to run the scripts:

- `fgsea`, version 1.20.0
- `tidyverse`, version >=2.0.0
- `HarmonicMeanP`, version >=3.0.1
- `jsonlite`, version >1.7.2

The analysis can be replicated using the following steps:

1. Run the `download-data.R` script to download the data from the NASA GeneLab repository.
2. Run the `enrichment-analysis.R` script to perform the gene set enrichment analysis using the downloaded data.