# Set working directory
setwd("./gene-set-enrichment")

# Load necessary libraries (add these at the beginning)
library(tidyverse)
library(fgsea)
library(dplyr)

# Read comparisons data
comparisons <- read.csv("./MetadataFiles/MouseComparisions.csv", header = TRUE)

# Process comparisons data
comparisons <- comparisons %>%
  mutate(
    Log2fc = sapply(comparisons$comparison, function(comp) paste0("Log2fc_", comp)),
    Pval = sapply(comparisons$comparison, function(comp) paste0("P.value_", comp)),
    AdjPval = sapply(comparisons$comparison, function(comp) paste0("Adj.p.value_", comp)),
    firstSD = sapply(comparisons$comparison, function(comp) paste0("Group.Stdev_", strsplit(comp, "\\.v\\.")[[1]][1], ".")),
    secondSD = sapply(comparisons$comparison, function(comp) paste0("Group.Stdev_.", strsplit(comp, "\\.v\\.")[[1]][2])),
    group = comparisons$comparison
  )

list.data <- unique(comparisons$dataset)

# Load all data
allData <- lapply(list.data, function(id) {
  message(id)
  f <- list.files(paste0("./data/", id), pattern = ".rds$", full.names = TRUE)
  if (file.exists(f)) {
    return(list(
      data = readRDS(f),
      metadata = read.table(list.files(paste0("./data/", id), pattern = ".txt$", full.names = TRUE)[[1]], 
                            header = TRUE, sep = "\t", fill = NA)
    ))
  }
  return(NULL)
})

allData <- allData[!sapply(allData, is.null)]
names(allData) <- list.data

# Function to process GMT files
gmt2geneset <- function(path) {
  genesets <- read_tsv(path, col_names = FALSE) %>% 
    apply(MARGIN = 1, function(r) {
      genes = unique(r[-(1:2)])
      path_name <- strsplit(r[[1]], ":")[[1]]
      
      list(
        id = paste(path_name, collapse = "_"),
        description = r[2],
        genes = genes[!is.na(genes)]
      )
    })
  
  gs <- lapply(genesets, function(g) g$genes %>% as.character())
  names(gs) <- sapply(genesets, function(g) g$id)
  gs
}

# Load gene set
geneset <- gmt2geneset("./gmt/MitoPathways3.gmt")

# Process statistical data
statDat <- lapply(names(allData), function(d) {
  message(d)
  data <- allData[[d]]$data
  cols_to_extract <- comparisons[comparisons$dataset == d,]
  
  lapply(1:nrow(cols_to_extract), function(rnum) {
    message(rnum)
    cur_row <- cols_to_extract[rnum,]
    
    if(any(sapply(c(cur_row$Log2fc, cur_row$Pval, cur_row$AdjPval), function(x) length(data[[x]]) == 0))) return(NULL)
    
    data_res <- data.frame(
      SYMBOL = str_to_upper(data$SYMBOL),
      log2FC = data[[cur_row$Log2fc]],
      pVal = data[[cur_row$Pval]],
      pAdj = data[[cur_row$AdjPval]],
      sdSF = data[[cur_row$firstSD]],
      sdGC = data[[cur_row$secondSD]],
      stringsAsFactors = FALSE
    ) %>% drop_na()
    
    list(
      id = d,
      name = cur_row$comparison,
      data = data_res
    )
  })
}) %>% unlist(recursive = FALSE)

statDat <- statDat[!sapply(statDat, is.null)]

# Perform GSEA
allRes <- lapply(statDat, function(data) {
  set.seed(2)
  
  data$data$pVal[data$data$pVal == 0] <- 1e-300
  geneRanks <- (-log10(data$data$pVal)) * sign(data$data$log2FC)
  names(geneRanks) <- data$data$SYMBOL
  geneRanks <- geneRanks[!is.na(geneRanks)]
  
  range_geneRanks <- range(geneRanks)
  if(all(is.finite(range_geneRanks))) {
    res <- fgseaSimple(pathways = geneset, stats = geneRanks, nperm = 1e+4, nproc = 1)
    
    final_res <- data.frame(
      pathway = names(geneset),
      stringsAsFactors = FALSE
    ) %>%
      left_join(res$pvals, by = "pathway") %>%
      mutate(
        pval = ifelse(is.na(pval), 1, pval),
        padj = ifelse(is.na(padj), 1, padj),
        NES = ifelse(is.na(NES), 0, NES),
        ES = ifelse(is.na(ES), 0, ES),
        pval = ifelse(pval == 0, 1e-6, ifelse(pval == 1, 1 - 1e-6, pval))
      )
    
    #final_res$HMDPval <- pharmonicmeanp(final_res$pval, L = length(final_res$pval), lower.tail = TRUE)
    
    list(
      id = data$id,
      name = data$name,
      data = data$data,
      res = final_res
    )
  }
})

allRes <- allRes[!sapply(allRes, is.null)]