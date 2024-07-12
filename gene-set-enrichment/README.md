# Gene set enrichment analysis and Meta analysis for mitochondrial pathways in space biology data

The `Codes` folder contains R scripts to download and process gene expression data from the NASA GeneLab repository for gene set enrichment analysis (GSEA) and Meta analysis using mitochondrial pathways. The `gmt` folder contains the defined genesets for different organisms. The `MetadataFiles` folder includes the required metadata files to collect the results. These metadata files are manually prepared utilizing the metadata provided by NASA GeneLab repository for each study.  The `Plots` folder contains all generated plots for the manuscript. The `Results` folder contains all the collected enrichment and meta analyses results, which are discussed in the manuscript. 

The following packages are required to run the scripts:

- `fgsea`, version 1.20.0
- `tidyverse`, version >=2.0.0
- `jsonlite`, version >1.7.2
- `metafor`, version 3.0-2
- `meta`, version 4.19-2
- `ggplot2`, version 3.3.3
- `gridExtra`, version 2.3
- `cowplot`, version 1.1.1
- `ggpubr`, version 0.4.0

The analysis can be replicated using the following steps:

1. Run the `download-data.R` script to download the data from the NASA GeneLab repository.
2. Run the `enrichment-analysis.R` script to perform the gene set enrichment analysis using the downloaded data.
3. Run the `meta-analysis.R` script to perform the meta analysis on mouse data using the gene set enrichment analysis results.
4. Run the `forest-plot.R` script to generate the forest plots for mice studies utilizing both gene set enrichment analysis and meta analysis results.