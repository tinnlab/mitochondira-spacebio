library(tidyverse)
library(stringr)
library(ggplot2)
library(dplyr)
library(grid)
library(gridExtra)
library(cowplot)
library(ggpubr)

setwd("./gene-set-enrichment")

resDf <- readRDS("./Results/CustomGSResults/MousePARes.rds")
metaRes <- readRDS("/Results/CustomGSResults/MetaAnalysis.rds")

mitoPathwayGroups <- read.table("./MetadataFiles/MitoPathwaysGroups.csv", sep = ",", header = TRUE, stringsAsFactors = F)
mitoPathwayGroups <- mitoPathwayGroups[mitoPathwayGroups$MitoPathway != "",]
mitoPathwayGroups$category <- lapply(mitoPathwayGroups$MitoPathway.Hierarchy, function(hierarchy) {
  hierarchy <- as.character(hierarchy)
  if (grepl(">", hierarchy)) {
    return(str_trim(strsplit(hierarchy, ">")[[1]][1]))
  } else return(hierarchy)
})

second.level <- mitoPathwayGroups %>% filter(category %in% c("Mitochondrial central dogma", "Metabolism"))
second.level$category.second <- lapply(second.level$MitoPathway.Hierarchy, function(hierarchy) {
  hierarchy <- as.character(hierarchy)
  if (grepl(">", hierarchy)) {
    return(str_trim(strsplit(hierarchy, ">")[[1]][2]))
  } else return(hierarchy)
})

mitoPathwayGroups$category[1:21] <- second.level$category.second[1:21]
mitoPathwayGroups$category[55:126] <- second.level$category.second[22:93]

resDf <- resDf %>%
  dplyr::mutate(category = mitoPathwayGroups$category[match(resDf$pathway, mitoPathwayGroups$MitoPathway)])

backupResDf <- resDf

for (i in 1:nrow(resDf)) {
  plotDat <- resDf[i,]
  print(i)

  plotDat$NES.sd <- abs((plotDat$NES - ifelse(plotDat$NES > 0, 1, -1)) / qnorm(plotDat$pval))

  sd <- plotDat$NES.sd

  sd[sd > 0.5] <- 0.5

  plotDat$min <- plotDat$NES - sd * 2
  plotDat$max <- plotDat$NES + sd * 2

  plotDat$min[plotDat$min < -2.5] <- -2.5
  plotDat$max[plotDat$max > 2.5] <- 2.5

  resDf[i, c("min")] <- plotDat$min
  resDf[i, c("max")] <- plotDat$max
  resDf[i, c("NES.sd")] <- plotDat$NES.sd
}

resDf$category <- as.vector(unlist(resDf$category))

categories <- unique(resDf$category)

sorted_datasets <- distinct(resDf[, c("summarized.name", "strain")])
sorted_datasets <- sorted_datasets[order(sorted_datasets$strain),]

resDf$summarized.name <- factor(resDf$summarized.name, levels = sorted_datasets$summarized.name[length(sorted_datasets$summarized.name):1])

strain_shapes <- c("BALB/c" = 21, "C57BL/6" = 22, "C3H/HeJ" = 23)

legend_data <- data.frame(
  Strain = c("BALB/c", "C57BL/6", "C3H/HeJ"),
  x = c(1, 2, 3),
  y = c(1, 2, 3)
)

legend_plot <- ggplot(legend_data, aes(x = x, y = y, shape = Strain)) +
  geom_point(size = 20) +
  scale_shape_manual(values = strain_shapes) +
  theme_minimal() +
  theme(legend.position = "top",
        legend.text = element_text(size = 40, face = "bold"),
        legend.title = element_text(size = 40, face = "bold"))

legend <- cowplot::get_legend(legend_plot)

plotsCategory <- lapply(categories, function(cat) {
  filteredDf <- resDf[resDf$category == cat,]
  sortedPathwaysCat <- filteredDf %>%
    group_by(pathway) %>%
    summarise(c = length(pathway)) %>%
    arrange(desc(c)) %>%
    `$`("pathway")

  plots <- filteredDf %>%
    group_by(tissue) %>%
    group_split() %>%
    lapply(function(plotDat) {
      plotDat <- data.frame(plotDat)
      if (length(unique(plotDat$summarized.name)) <= 1) return(NA)
      dsColor <- ggsci::pal_npg("nrc", alpha = 1)(length(unique(plotDat$summarized.name)))
      names(dsColor) <- unique(plotDat$summarized.name)
      plotDat$pathway <- factor(plotDat$pathway, levels = sortedPathwaysCat[length(sortedPathwaysCat):1])

      ggplot(plotDat, aes(x = NES, y = pathway, group = summarized.name, col = summarized.name, shape = strain)) +
        theme_minimal() +
        geom_point(size = 7, position = position_dodge(0.7), aes(fill = summarized.name), show.legend = FALSE) +
        geom_rect(
          aes(
            xmin = -Inf, xmax = Inf, ymin = as.numeric(pathway) - 0.5, ymax = as.numeric(pathway) + 0.5
          ),
          fill = ifelse((as.numeric(plotDat$pathway) %% 2 == 0), "white", "#eeeeee"),
          color = "white"
        ) +
        geom_point(size = 7, position = position_dodge(0.7), aes(fill = summarized.name)) +
        geom_vline(xintercept = c(-1, 1), colour = "#FA8072", linetype = "longdash") +
        geom_vline(xintercept = c(0), colour = "grey", linetype = "solid") +
        coord_cartesian(clip = "off", xlim = c(-2.5, 2.5)) +
        geom_errorbarh(height = .3, width = 1, aes(xmin = min, xmax = max), position = position_dodge(0.7)) +
        theme_bw() +
        theme(
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.margin = unit(c(5, 5, 7, 5), "pt"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text = element_text(size = 40, face = "bold"),
          legend.position = "top",
          legend.title = element_blank(),
          legend.text = element_text(size = 40, face = "bold"),
          axis.text = element_text(size = 40, face = "bold"),
          axis.text.x = element_text(vjust = -0.75),
          legend.key = element_rect(colour = NA, fill = NA),
          legend.box.background = element_blank()
        ) +
        labs(x = "NES") +
        guides(
          colour = guide_legend(ncol = 2),
          shape = FALSE,
          fill = FALSE
        ) +
        scale_shape_manual(values = strain_shapes)
    })
  plots <- plots[!is.na(plots)]
})

names(plotsCategory) <- categories

metaRes <- metaRes %>%
  dplyr::mutate(category = mitoPathwayGroups$category[match(metaRes$pathway, mitoPathwayGroups$MitoPathway)])

plotsMetaCategory <- lapply(categories, function(cat) {
  filtered.MetaRes <- metaRes[metaRes$category == cat,]
  sortedPathways <- filtered.MetaRes %>%
    group_by(pathway) %>%
    summarise(c = length(pathway)) %>%
    arrange(desc(c)) %>%
    `$`("pathway")

  plts.f <- filtered.MetaRes %>%
    group_by(tissue) %>%
    group_split() %>%
    lapply(function(plotDat) {
      plotDat <- data.frame(plotDat)
      plotDat$pathway <- factor(plotDat$pathway, levels = sortedPathways[length(sortedPathways):1])

      xxx <<- plotDat

      sd <- plotDat$NES.combined.sd

      sd[sd > 0.5] <- 0.5

      plotDat$min <- plotDat$NES.combined - sd * 2
      plotDat$max <- plotDat$NES.combined + sd * 2

      plotDat$min[plotDat$min < -2.5] <- -2.5
      plotDat$max[plotDat$max > 2.5] <- 2.5

      ggplot(plotDat, aes(y = pathway, x = NES.combined, xmin = min, xmax = max)) +
        theme_minimal() +
        geom_point(size = 7, color = "red") +
        geom_rect(
          aes(
            xmin = -Inf, xmax = Inf, ymin = as.numeric(pathway) - 0.5, ymax = as.numeric(pathway) + 0.5
          ),
          fill = ifelse((as.numeric(plotDat$pathway) %% 2 == 0), "white", "#eeeeee"),
          color = "white"
        ) +
        geom_point(size = 7, color = "red") +
        geom_vline(xintercept = c(-1, 1), colour = "#FA8072", linetype = "longdash") +
        geom_vline(xintercept = c(0), colour = "grey", linetype = "solid") +
        coord_cartesian(clip = "off", xlim = c(-2.5, 2.5)) +
        geom_errorbarh(height = .2, width = 1) +
        theme_bw() +
        theme(
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.margin = unit(c(5, 5, 7, 5), "pt"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text = element_text(size = 40, face = "bold"),
          plot.title = element_text(size = 50, face = "bold", hjust = 0.5),
          legend.position = "top",
          legend.title = element_blank(),
          legend.text = element_text(size = 40, face = "bold"),
          axis.text = element_text(size = 40, face = "bold", hjust = 0.5),
          axis.text.x = element_text(vjust = -0.75),
          legend.key = element_rect(colour = NA, fill = NA),
          legend.box.background = element_blank()
        ) +
        labs(x = "NES", title = str_to_title(na.omit(plotDat$tissue)[[1]]))
    })
})

names(plotsMetaCategory) <- categories

plotsAllCategory <- list()

for (cat in categories) {
  print(cat)
  plots1 <- plotsCategory[names(plotsCategory) == cat][[1]]
  plots2 <- plotsMetaCategory[names(plotsMetaCategory) == cat][[1]]
  allPlots.f <- list()
  j <- 1
  for (i in 1:length(plots1)) {
    print(j)
    allPlots.f[[j]] <- plots1[[i]]
    j <- j + 1
    allPlots.f[[j]] <- plots2[[i]]
    j <- j + 1
  }
  plotsAllCategory[[cat]] <- allPlots.f
}

names(plotsAllCategory) <- categories

plotsAllLegends <- list()

for (i in 1:length(plotsAllCategory)) {
  currentCategory <- plotsAllCategory[[i]]
  plotTitle <- names(plotsAllCategory)[[i]]

  tempList <- list()
  for (j in 1:length(currentCategory)) {
    p <- currentCategory[[j]]
    legendP <- cowplot::get_legend(p)
    titleP <- cowplot::get_title(p)

    l <- list(legend = legendP,
              title = titleP)
    tempList[[j]] <- l
  }

  plotsAllLegends[[plotTitle]] <- tempList
}

plotsAllCategoryWithoutLegend <- lapply(plotsAllCategory, function(cat) {
  lapply(cat, function(plot) {
    plot <- plot + theme(legend.position = "none", title = element_blank(), plot.title = element_blank())
  })
})

legendIDs <- c(1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23)
titleIDs <- c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24)
totalLength <- length(legendIDs) + length(titleIDs)

## Due to the differences in length of geneset names the dimensions should be fixed for each plot manually.
## The below setting is for i=2 where category is Mtdna Maintenance
#for (i in 1:length(plotsAllCategoryWithoutLegend)) {
i <- 2
currentCategory <- plotsAllCategoryWithoutLegend[[i]]
plotTitle <- names(plotsAllCategoryWithoutLegend)[i]
plotTitle <- str_to_title(plotTitle)
curLegend <- plotsAllLegends[[i]]

titleGrob <- grid::textGrob(plotTitle, gp = grid::gpar(fontsize = 60, fontface = 'bold', color = "black"))

plot1 <- list()
for (j in 1:totalLength) {
  if ((j %% 2) == 0) {
    plot1[[j]] <- curLegend[[j]]$title
  }else {
    plot1[[j]] <- curLegend[[j]]$legend
  }
}

##
onlyLegends <- lapply(legendIDs, function(id) curLegend[[id]]$legend)
onlyTitles <- lapply(titleIDs, function(id) curLegend[[id]]$title)

gs_legends <- lapply(onlyLegends, function(p)
  ggpubr::as_ggplot(p) + theme(
    legend.margin = margin(c(0, 0, 0, 0))))

gs_titles <- lapply(onlyTitles, function(p)
  ggpubr::as_ggplot(p) + theme(
    legend.margin = margin(c(0, 0, 0, 0))))

listAll <- c(gs_titles[1:6], c(gs_legends[1:6],
                               list(
                                 ggplot(currentCategory[[1]]$data, aes(y = pathway, x = "")) +
                                   labs(y = "", title = "", x = "") +
                                   theme_minimal() +
                                   theme(panel.grid.major.x = element_blank(), plot.margin = unit(c(5, 5, 5, 5), "pt"), text = element_text(size = 50, face = "bold", colour = "black"))
                               ),
                               currentCategory[1:12], gs_titles[7:12], gs_legends[7:12],
                               list(
                                 ggplot(currentCategory[[1]]$data, aes(y = pathway, x = "")) +
                                   labs(y = "", title = "", x = "") +
                                   theme_minimal() +
                                   theme(panel.grid.major.x = element_blank(), plot.margin = unit(c(5, 5, 5, 5), "pt"), text = element_text(size = 50, face = "bold", colour = "black"))
                               ),
                               currentCategory[13:24]
))

g <- c(list(titleGrob), list(legend), listAll)

##
layout_matrix <- rbind(
  c(rep(1, 13)),
  c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  c(rep(2, 13)),
  c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  c(NA, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8),
  c(NA, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 14),
  c(rep(15:27)),
  c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  c(NA, 28, 28, 29, 29, 30, 30, 31, 31, 32, 32, 33, 33),
  c(NA, 34, 34, 35, 35, 36, 36, 37, 37, 38, 38, 39, 39),
  c(rep(40:52)),
  c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
)

pdf(file = paste0("./Plots/", plotTitle, ".pdf"), width = 110, height = 50)

xx <- gridExtra::grid.arrange(
  grobs = g,
  widths = c(6, rep(5, 12)),
  heights = c(0.2, 0.2, 0.2, 0.3, 0.08, 0.5, 2.5, 0.3, 0.08, 0.7, 2.5, 0.1),
  layout_matrix = layout_matrix
)

dev.off()
#}
