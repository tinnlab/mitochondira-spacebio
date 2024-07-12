library(dplyr)
library(stringr)
library(tidyverse)
library(fgsea)
library(meta)

mouse_studies <- c("OSD-162", "OSD-164", "OSD-173", "OSD-194", "OSD-238", "OSD-239", "OSD-240", "OSD-241", "OSD-242",
                   "OSD-243", "OSD-244", "OSD-245", "OSD-246", "OSD-248", "OSD-253", "OSD-254", "OSD-255", "OSD-270",
                   "OSD-288", "OSD-289", "OSD-326", "OSD-379", "OSD-401", "OSD-419", "OSD-420", "OSD-421", "OSD-576",
                   "OSD-580", "OSD-352", "OSD-98", "OSD-99", "OSD-163", "OSD-161", "OSD-104", "OSD-102", "OSD-101",
                   "OSD-100", "OSD-137", "OSD-47", "OSD-48", "OSD-103", "OSD-105", "OSD-21", "OSD-25", "OSD-111",
                   "OSD-135", "OSD-4", "OSD-87")

setwd("./gene-set-enrichment")

mouse_comparisons <- read.csv("./MetadataFiles/MouseComparisions.csv", stringsAsFactors = FALSE)

# Process comparisons data
comparisons <- mouse_comparisons %>%
  mutate(
    Log2fc = sapply(mouse_comparisons$comparison, function(comp) paste0("Log2fc_", comp)),
    Pval = sapply(mouse_comparisons$comparison, function(comp) paste0("P.value_", comp)),
    AdjPval = sapply(mouse_comparisons$comparison, function(comp) paste0("Adj.p.value_", comp)),
    firstSD = sapply(mouse_comparisons$comparison, function(comp) paste0("Group.Stdev_", strsplit(comp, "\\.v\\.")[[1]][1], ".")),
    secondSD = sapply(mouse_comparisons$comparison, function(comp) paste0("Group.Stdev_.", strsplit(comp, "\\.v\\.")[[1]][2])),
    group = mouse_comparisons$comparison
  )

# Load all data
allData <- lapply(mouse_studies, function(id) {
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
names(allData) <- mouse_studies

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
geneset <- gmt2geneset("./gmt/24_1.22_custom_Mitochondria_SpaceBio_Pathways-Space6.gmt")

# Process statistical data
statDat <- lapply(names(allData), function(d) {
  message(d)
  data <- allData[[d]]$data
  cols_to_extract <- comparisons[comparisons$dataset == d,]

  lapply(1:nrow(cols_to_extract), function(rnum) {
    message(rnum)
    cur_row <- cols_to_extract[rnum,]

    if (any(sapply(c(cur_row$Log2fc, cur_row$Pval, cur_row$AdjPval), function(x) length(data[[x]]) == 0))) return(NULL)

    data <- data.frame(
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
      data = data
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
  if (all(is.finite(range_geneRanks))) {
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

    list(
      id = data$id,
      name = data$name,
      data = data$data,
      res = final_res
    )
  }
})

allRes <- allRes[!sapply(allRes, is.null)]

#Meta analysis
#Due to the variance in studies metadata, we meticulously collected the required information by thoroughly reviewing each study.
mouse_metadata <- read.csv("./MetadataFiles/MouseMetadata.csv", stringsAsFactors = FALSE)

resDf <- allRes %>% lapply(function(res) {
  cur_meta <- mouse_metadata %>% filter(dataset == res$id, comparision == res$name)
    xx <- data.frame(
      dataset = res$id,
      summarized.name = cur_meta$summarized.name,
      tissue = tolower(cur_meta$tissue),
      strain = cur_meta$strain,
      sample.size = cur_meta$n.sample,
      res$res,
      stringsAsFactors = F
    )
  }) %>% do.call(what = rbind)


# tissue map
tissueMap = list(
  "leg muscle" = c("calf muscle", "extensor digitorum longus- both sides", "gastrocnemius-left", "quadriceps-left", "soleus-both sides", "tibialis anterior-left", "right tibialis anterior", "gastrocnemius-right", "quadriceps femoris", "soleus"),
  "bone marrow" = c("bone marrow cells"),
  "kidney" = c("kidney-left", "left kidney"),
  "eye" = c("eye-left"),
  "skin" = c("dorsal skin", "femoral lateral skin", "femoral skin"),
  "adrenal gland" = c("adrenal glands- both sides"),
  "retina" = c("left retina", "right retina"),
  "back muscle" = c("longissimus dorsi muscle"),
  "spleen" = c("spleen-distal"),
  "thymus" = c("thymus gland"),
  "lung" = c("left lung lobe")
)

for (tissue in names(tissueMap)) {
  resDf$tissue[resDf$tissue %in% tissueMap[[tissue]]] <- tissue
}

# strain map
strainMap = list(
  "C57BL/6" = c("C57BL/6NTac", "C57BL/6NCrl", "C57BL/6Tac", "C57BL/6J", "C57BL/6CR", "C57/BL6", "C57BL/6N"),
  "BALB/c" = c("BALB/cAnNTac")
)

for (strain in names(strainMap)) {
  resDf$strain[resDf$strain %in% strainMap[[strain]]] <- strain
}

metaRes <- resDf %>%
  group_by(tissue, pathway) %>%
  group_split() %>%
  lapply(function(resDfsm) {
    nes_var <- var(resDfsm$NES)
    if (is.na(nes_var) || nes_var == 0) return(NA)
    resDfsm$NES.sd <- abs((resDfsm$NES - ifelse(resDfsm$NES > 0, 1, -1)) / qnorm(resDfsm$pval))
    message(paste0(resDfsm$dataset[1], "  ", resDfsm$pathway[1]))
    res <- meta::metagen(data = resDfsm, studlab = pathway, TE = NES, seTE = NES.sd, sm = "SMD",
                         n.e = sample.size,
                         method.tau = "REML",
                         hakn = TRUE)

    NES.combined <- res$TE.fixed
    NES.combined.sd <- res$seTE.fixed

    pval <- pnorm((ifelse(NES.combined > 0, 1, -1) - NES.combined) / NES.combined.sd)
    if (NES.combined < 0) pval <- 1 - pval

    if (abs(NES.combined) < 1 & pval < 0.5) {
      xxx <<- resDfsm
    }

    data.frame(
    tissue = resDfsm$tissue[1],
    pathway = resDfsm$pathway[1],
    NES.combined = NES.combined,
    NES.combined.sd = NES.combined.sd,
    pval.combined = pval,
    count = nrow(resDfsm),
    stringsAsFactors = F
  )
  })

metaRes <- metaRes[!is.na(metaRes)] %>% do.call(what = rbind) %>% as.data.frame()

metaRes <- metaRes %>% group_by(tissue) %>% group_split() %>% lapply(function(r){
  r$p.fdr <- p.adjust(r$pval.combined, method = "fdr")
  r
}) %>% do.call(what = rbind)

