library(chipmine)
library(org.Anidulans.eg.db)
library(esquisse)
library(summarytools)
library(here)
library(ggpubr)

## 1) plots SM cluster wise genes binding signal plot as geom_tiles
## 2) plots cluster wise binding signal and polII signal fold change plot using geom_tile
## 3) plots the binding signal and polII fold change values in two groups: bound and unbound
rm(list = ls())

##################################################################################
analysisName <- "kdmB_del"
outPrefix <- here::here("kdmB_analysis/SM_analysis", analysisName)

tfIds <- c("An_kdmB_20h_HA_1", "An_kdmB_48h_HA_1")
polII1 <- "An_untagged_20h_polII_1"
polII2 <- "An_untagged_48h_polII_1"

otherPolIIs <- c("An_kdmB_del_20h_polII_1", "An_kdmB_del_48h_polII_1")

## polII signal fold change pairs
polIIDiffPairs <- list(
  p1 = list(
    name = "48h_vs_20h_untagged_polII",
    title = "polII log2(untagged_48h \n vs untagged_20h)",
    samples = c(polII2, polII1)
  ),
  p2 = list(
    name = "48h_vs_20h_kdmB_del_polII",
    title = "polII log2(kdmB_del_48h \n vs kdmB_del_20h)",
    samples = c(otherPolIIs[2], otherPolIIs[1])
  ),
  p3 = list(
    name = "kdmB_del_vs_untagged_20h",
    title = "polII log2(kdmB_del \n vs kdmB_untagged 20h)",
    samples = c("An_kdmB_del_20h_polII_1", "An_untagged_20h_polII_1")
  ),
  p4 = list(
    name = "kdmB_del_vs_untagged_48h",
    title = "polII log2(kdmB_del \n vs kdmB_untagged 48h)",
    samples = c("An_kdmB_del_48h_polII_1", "An_untagged_48h_polII_1")
  )
)

polIIPairId <- "p1"

## genes to read
file_exptInfo <- here::here("data", "referenceData/sampleInfo.txt")
file_genes <- here::here("data", "referenceData/AN_genesForPolII.bed")
file_topGoMap <- "E:/Chris_UM/Database/A_Nidulans/ANidulans_OrgDb/geneid2go.ANidulans.topGO.map"
file_geneInfo <- "E:/Chris_UM/Database/A_Nidulans/A_nidulans_FGSC_A4_geneClasses.txt"

TF_dataPath <- here::here("data", "TF_data")
polII_dataPath <- here::here("data", "polII_data")
hist_dataPath <- here::here("data", "histone_data")


orgDb <- org.Anidulans.eg.db

##################################################################################

## genes to read
geneSet <- data.table::fread(file = file_genes, header = F,
                             col.names = c("chr", "start", "end", "gene", "score", "strand")) %>% 
  dplyr::mutate(length = end - start)

geneDesc <- AnnotationDbi::select(x = orgDb, keys = geneSet$gene, columns = c("SM_CLUSTER"), keytype = "GID")

geneSet <- dplyr::left_join(x = geneSet, y = geneDesc, by = c("gene" = "GID"))


##################################################################################


tfInfo <- get_sample_information(exptInfoFile = file_exptInfo,
                                 samples = tfIds,
                                 dataPath = TF_dataPath,
                                 matrixSource = "normalizedmatrix")

polII_info <- get_sample_information(exptInfoFile = file_exptInfo,
                                     samples = c(polII1, polII2, otherPolIIs),
                                     dataPath = polII_dataPath,
                                     matrixSource = "normalizedmatrix")

exptData <- dplyr::bind_rows(tfInfo, polII_info)

polII_ids <- exptData$sampleId[which(exptData$IP_tag == "polII")]
tfIds <- exptData$sampleId[which(exptData$IP_tag %in% c("HA", "MYC", "TAP") & exptData$TF != "untagged")]


polIICols <- list(
  exp = structure(polII_ids, names = polII_ids),
  is_expressed = structure(paste("is_expressed", ".", polII_ids, sep = ""), names = polII_ids)
)

tfCols <- sapply(c("hasPeak", "pval", "peakType", "tesPeakType", "peakDist", "summitDist", "upstreamExpr", "peakExpr", "relativeDist"),
                 FUN = function(x){ structure(paste(x, ".", tfIds, sep = ""), names = tfIds) },
                 simplify = F, USE.NAMES = T)


##################################################################################

chipData <- get_polII_expressions(genesDf = geneSet, exptInfo = polII_info)

## add fold change columns
for (i in names(polIIDiffPairs)) {
  chipData <- get_fold_change(df = chipData,
                              s1 = polIIDiffPairs[[i]]$samples[1],
                              s2 = polIIDiffPairs[[i]]$samples[2],
                              newCol = polIIDiffPairs[[i]]$name,
                              isExpressedCols = polIICols$is_expressed)
}


chipData <- get_TF_binding_data(genesDf = chipData, exptInfo = tfInfo, allColumns = FALSE) %>% 
  dplyr::filter(! is.na(SM_CLUSTER)) %>% 
  dplyr::mutate(SM_CLUSTER = gsub(pattern = "SM_cluster_", replacement = "", x = SM_CLUSTER, fixed = TRUE)) %>% 
  dplyr::group_by(SM_CLUSTER) %>% 
  dplyr::arrange(start, .by_group = TRUE) %>% 
  dplyr::mutate(index = 1:n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(gene, SM_CLUSTER, starts_with("hasPeak."), index, !! unname(purrr::map_chr(polIIDiffPairs, "name"))) %>% 
  as.data.frame()

view(dfSummary(chipData))


##################################################################################


## SM cluster binding plot
pltDf1 <- tidyr::gather(chipData, key = "sample", value = "hasPeak",
                        -gene, -SM_CLUSTER, -index, - ends_with("_polII"))

pltDf1$sample <- factor(pltDf1$sample, levels = c(unname(tfCols$hasPeak)))

pltTitle <- paste("SM cluster binding comparison:", analysisName)

pt1 <- ggplot(data = pltDf1) +
  geom_tile(mapping = aes(x = index, y = sample, fill = hasPeak, color = sample), size = 0.75, height = 0.9) +
  scale_fill_manual(
    values = c("TRUE" = "black", "FALSE" = "white"),
    guide = FALSE
  ) +
  scale_color_discrete(
    breaks = unname(tfCols$hasPeak),
    # values = structure(c("#B35806", "#542788"), names = unname(tfCols$hasPeak))
    name = ""
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  facet_wrap(facets = SM_CLUSTER ~ ., scales = "free_y", ncol = 2, strip.position = "left", dir = "v") +
  ggtitle(pltTitle) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.spacing = unit(0.15, "lines"),
        strip.text.y = element_text(hjust = 0.5, size = 14, face = "bold", angle = 180),
        strip.background = element_rect(fill="white"),
        legend.text = element_text(size = 13),
        legend.position = "bottom",
        legend.title = element_text(size = 13, face = "bold"),
        plot.margin = unit(rep(0.5, 4), "cm"))


pdf(file = paste(outPrefix, "_cluster_binding_cmp.pdf", sep = ""), width = 10, height = 10)
pt1
dev.off()

##################################################################################
## SM cluster polII signal plot
pltDf2 <- tidyr::gather(chipData, key = "sample", value = "value",
                        -gene, -SM_CLUSTER, -index)


pltDf2$sample <- factor(pltDf2$sample,
                        levels = unname(c(tfCols$hasPeak[1], purrr::map_chr(polIIDiffPairs, "name"), tfCols$hasPeak[2])))


pltTitle <- paste(c("SM cluster polII fold change:", analysisName), collapse = " ")

pt2 <- ggplot() +
  geom_segment(data = dplyr::filter(pltDf2, sample %in% unname(tfCols$hasPeak)) %>% 
                 dplyr::mutate(value = as.character(value)),
               mapping = aes(x = index - 0.5, xend = index + 0.5, y = sample, yend = sample, color = value),
               size = 5) +
  geom_tile(data = dplyr::filter(pltDf2, sample %in% unname(purrr::map_chr(polIIDiffPairs, "name"))),
            mapping = aes(x = index, y = sample, fill = value),
            color = "black", size = 0.5, height = 1) +
  scale_fill_gradient2(
    name = paste("log2(", "polII fold change", ")", sep = ""),
    low = "#B35806", mid = "#F7F7F7", high = "#542788", midpoint = 0
  ) +
  scale_colour_manual(
    name = "",
    values = c("1" = "#006600", "0" = "white"),
    breaks = c("1"),
    labels = c("peak detected")) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(limits = unname(c(tfCols$hasPeak[1], purrr::map_chr(polIIDiffPairs, "name"), tfCols$hasPeak[2])),
                   expand = expand_scale(add = 0.0)) +
  facet_wrap(facets = SM_CLUSTER ~ ., scales = "free_y", ncol = 4, strip.position = "left", dir = "v",
             labeller = labeller(vs = label_both, am = label_value)) +
  ggtitle(pltTitle) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.spacing = unit(0.15, "lines"),
        strip.text.y = element_text(hjust = 0.5, size = 14, face = "bold", angle = 180),
        strip.background = element_rect(fill="white"),
        legend.text = element_text(size = 13),
        legend.position = "bottom",
        legend.title = element_text(size = 13, face = "bold"),
        plot.margin = unit(rep(0.5, 4), "cm"))


# png(filename = paste(outPrefix, "_SM_cluster_polII_diff.png"), width = 4000, height = 6000, res = 350)
pdf(file = paste(outPrefix, "_SM_cluster_polII_diff.pdf", sep = ""), width = 15, height = 8)
pt2
dev.off()


##################################################################################
## plot tiles in two groups, showing kdmB peak and not showing peak
pt3Df <- dplyr::mutate(.data = chipData,
                       hasPeak = !! as.name(tfCols$hasPeak[1]) | !! as.name(tfCols$hasPeak[2]) ) %>% 
  dplyr::arrange_at(.vars = unname(tfCols$hasPeak), .funs = funs(desc(.)))

## convert gene column to factor 
pt3Df$gene <- factor(pt3Df$gene, levels = pt3Df$gene)

dplyr::group_by(pt3Df, hasPeak) %>% 
  dplyr::summarise(n = n())

pt3 <- ggplot() +
  geom_segment(
    data = tidyr::gather(pt3Df, key = "tf", value = "tfPeak",
                         -gene, -SM_CLUSTER, -index, -hasPeak, -ends_with("_polII")) %>% 
      dplyr::mutate(tfPeak = as.logical(tfPeak)),
    mapping = aes(x = 0.5, xend = 1.5, y = tf, yend = tf,
                  color = tfPeak),
    size = 3
  ) +
  geom_tile(
    data = tidyr::gather(pt3Df, key = "polII", value = "lfc",
                         -gene, -SM_CLUSTER, -index, -hasPeak, -starts_with("hasPeak.")),
    mapping = aes(x = 1, y = polII, fill = lfc),
    color = "black", size = 0.5, height = 1
  ) +
  scale_fill_gradient2(
    name = paste("log2(", "polII fold change", ")", sep = ""),
    low = "#B35806", mid = "#F7F7F7", high = "#542788", midpoint = 0
  ) +
  scale_colour_manual(
    name = "",
    values = c("TRUE" = "#006600", "FALSE" = "white"),
    breaks = c("TRUE"),
    labels = c("peak detected")) +
  scale_y_discrete(limits = unname( c(tfCols$hasPeak[1], purrr::map_chr(polIIDiffPairs, "name"), tfCols$hasPeak[2]) ),
                   expand = expand_scale(add = 0.0)) +
  scale_x_continuous(expand = c(0, 0)) +
  facet_wrap(facets = gene ~ ., scales = "free", nrow = 16, strip.position = "left", dir = "h") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.spacing = unit(0.15, "lines"),
        strip.text.y = element_blank(),
        strip.background = element_rect(fill="white"),
        legend.text = element_text(size = 13),
        legend.position = "bottom",
        legend.title = element_text(size = 13, face = "bold"),
        plot.margin = unit(rep(0.5, 4), "cm"))


## plot list for arranging side by side
ptList <- dplyr::group_by(pt3Df, hasPeak) %>%
  dplyr::do(
    n = as.numeric(nrow(.)),
    pt = ggplot() +
      geom_segment(
      data = tidyr::gather(., key = "tf", value = "tfPeak",
                           -gene, -SM_CLUSTER, -index, -hasPeak, -ends_with("_polII")) %>% 
        dplyr::mutate(tfPeak = as.logical(tfPeak)),
      mapping = aes(x = 0.5, xend = 1.5, y = tf, yend = tf,
                    color = tfPeak),
      size = 3
    ) +
      geom_tile(
        data = tidyr::gather(., key = "polII", value = "lfc",
                             -gene, -SM_CLUSTER, -index, -hasPeak, -starts_with("hasPeak.")),
        mapping = aes(x = 1, y = polII, fill = lfc),
        color = "black", size = 0.5, height = 1
      ) +
      scale_fill_gradient2(
        name = paste("log2(", "polII fold change", ")", sep = ""),
        low = "#B35806", mid = "#F7F7F7", high = "#542788", midpoint = 0
      ) +
      scale_colour_manual(
        name = "",
        values = c("TRUE" = "#006600", "FALSE" = "white"),
        breaks = c("TRUE"),
        labels = c("peak detected")) +
      scale_y_discrete(limits = unname( c(tfCols$hasPeak[1], purrr::map_chr(polIIDiffPairs, "name"), tfCols$hasPeak[2]) ),
                       expand = expand_scale(add = 0.0)) +
      scale_x_continuous(expand = c(0, 0)) +
      facet_wrap(facets = gene ~ ., scales = "free", nrow = 12, strip.position = "left", dir = "h") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
            axis.text = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            panel.grid = element_blank(),
            panel.spacing = unit(0.15, "lines"),
            strip.text.y = element_blank(),
            strip.background = element_rect(fill="white"),
            legend.text = element_text(size = 13),
            legend.position = "bottom",
            legend.title = element_text(size = 13, face = "bold"),
            plot.margin = unit(rep(1, 4), "cm"))
  ) %>% 
  dplyr::arrange(desc(hasPeak))



pt3 <- ggpubr::ggarrange(plotlist = ptList$pt,
                         ncol = 2, nrow = 1,
                         labels = paste("hasPeak:", ptList$hasPeak),
                         common.legend = TRUE, legend = "bottom",
                         widths = unlist(ptList$n)
)


# png(filename = paste(outPrefix, "_SM_cluster_hasPeak_pairs.png", sep = ""), width = 6000, height = 6000, res = 550)
pdf(file = paste(outPrefix, "_SM_cluster_hasPeak_pairs.pdf", sep = ""), width = 14, height = 8)
pt3
dev.off()






