library(chipmine)
library(org.Anidulans.eg.db)
library(esquisse)
library(summarytools)
library(here)



rm(list = ls())

##################################################################################
analysisName <- "kdmB_SM"
outPrefix <- here::here("kdmB_analysis/SM_analysis", analysisName)

tfIds <- c("An_kdmB_20h_HA_1", "An_kdmB_48h_HA_1")
polII1 <- "An_untagged_20h_polII_1"
polII2 <- "An_untagged_48h_polII_1"


## polII signal fold change pairs
polIIDiffPairs <- list(
  p1 = list(
    name = "48h_vs_20h_untagged_polII",
    title = "polII log2(untagged_48h \n vs untagged_20h)",
    samples = c(polII2, polII1)
  )
)



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
                                     samples = c(polII1, polII2),
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
  dplyr::group_by(SM_CLUSTER) %>% 
  dplyr::arrange(start, .by_group = TRUE) %>% 
  dplyr::mutate(index = 1:n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(gene, SM_CLUSTER, starts_with("hasPeak."), index, polIIDiffPairs$p1$name) %>% 
  as.data.frame()

view(dfSummary(chipData))

pltDf <- tidyr::gather(chipData, key = "sample", value = "hasPeak",
                       -gene, -SM_CLUSTER, -index, - !!as.name(polIIDiffPairs$p1$name)) %>% 
  dplyr::mutate(newId = paste(SM_CLUSTER, sample, sep = "_"))


pltDf$sample <- factor(pltDf$sample, levels = rev(unname(tfCols$hasPeak)))

pt <- ggplot(data = pltDf) +
  geom_tile(mapping = aes(x = index, y = sample, fill = hasPeak, color = sample), size = 0.75, height = 0.9) +
  scale_fill_manual(
    values = c("TRUE" = "black", "FALSE" = "white"),
    guide = FALSE
  ) +
  scale_color_discrete(
    guide = guide_legend(reverse=TRUE)
  ) +
  # facet_grid(SM_CLUSTER ~ ., scales = "free_y", switch = "y") +
  facet_wrap(facets = SM_CLUSTER ~ ., scales = "free_y", ncol = 2, strip.position = "left") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 1, size = 14, face = "bold"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.spacing = unit(0.25, "lines"),
        strip.text.y = element_text(hjust = 0.5, size = 14, face = "bold", angle = 180),
        strip.background = element_rect(fill="white"),
        legend.text = element_text(size = 13),
        legend.position = "bottom",
        legend.title = element_text(size = 13, face = "bold"),
        plot.margin = unit(rep(0.5, 4), "cm"))



# png(filename = paste(outPrefix, "_cluster_binding_cmp.png"), width = 4000, height = 6000, res = 350)
pdf(file = paste(outPrefix, "_cluster_binding_cmp.pdf"), width = 10, height = 10)
pt
dev.off()















