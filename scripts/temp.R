library(chipmine)
library(org.Anidulans.eg.db)
library(esquisse)
library(summarytools)

rm(list = ls())



path <- "E:/Chris_UM/Analysis/21_CL2017_ChIPmix_ULAS_MIX/ULAS_AN/kdmB_analysis/SM_analysis/"
setwd(path)


## genes to read
file_exptInfo <-"E:/Chris_UM/Analysis/21_CL2017_ChIPmix_ULAS_MIX/ULAS_AN/data/referenceData/sampleInfo.txt"
file_genes <- "E:/Chris_UM/Analysis/21_CL2017_ChIPmix_ULAS_MIX/ULAS_AN/data/referenceData/AN_genesForPolII.bed"
file_topGoMap <- "E:/Chris_UM/Database/A_Nidulans/ANidulans_OrgDb/geneid2go.ANidulans.topGO.map"
file_topGoMap2 <- "E:/Chris_UM/Database/A_Nidulans/ANidulans_OrgDb/geneid2go.ANidulans.topGO.map2"
file_geneInfo <- "E:/Chris_UM/Database/A_Nidulans/A_nidulans_FGSC_A4_geneClasses.txt"

TF_dataPath <- "E:/Chris_UM/Analysis/21_CL2017_ChIPmix_ULAS_MIX/ULAS_AN/data/TF_data"
polII_dataPath <- "E:/Chris_UM/Analysis/21_CL2017_ChIPmix_ULAS_MIX/ULAS_AN/data/polII_data"
hist_dataPath <- "E:/Chris_UM/Analysis/21_CL2017_ChIPmix_ULAS_MIX/ULAS_AN/data/histone_data"

orgDb <- org.Anidulans.eg.db

tfIds <- c("An_kdmB_20h_HA_1", "An_kdmB_48h_HA_1")
##################################################################################

## genes to read
geneSet <- data.table::fread(file = file_genes, header = F,
                             col.names = c("chr", "start", "end", "gene", "score", "strand")) %>% 
  dplyr::mutate(length = end - start)

geneDesc <- AnnotationDbi::select(x = orgDb, keys = geneSet$gene, columns = c("SM_CLUSTER"), keytype = "GID")

geneSet <- dplyr::left_join(x = geneSet, y = geneDesc, by = c("gene" = "GID")) %>% 
  dplyr::filter(! is.na(SM_CLUSTER))


##################################################################################

tfInfo <- get_sample_information(exptInfoFile = file_exptInfo,
                                 samples = tfIds,
                                 dataPath = TF_dataPath,
                                 matrixSource = "normalizedmatrix")


tfData <- get_TF_binding_data(genesDf = geneSet, exptInfo = tfInfo, allColumns = FALSE) %>% 
  dplyr::group_by(SM_CLUSTER) %>% 
  dplyr::arrange(start, .by_group = TRUE) %>% 
  dplyr::mutate(index = 1:n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(gene, SM_CLUSTER, starts_with("hasPeak."), index) %>% 
  as.data.frame()

view(dfSummary(tfData))

pltDf <- tidyr::gather(tfData, key = "sample", value = "hasPeak", -gene, -SM_CLUSTER, -index) %>% 
  dplyr::mutate(newId = paste(SM_CLUSTER, sample, sep = "_"))
  

pt <- ggplot(data = pltDf, mapping = aes(x = index, y = sample)) +
  geom_tile(mapping = aes(fill = hasPeak, color = sample), size = 0.75, height = 0.9) +
  scale_fill_manual(
    values = c("TRUE" = "black", "FALSE" = "white")
  ) +
  facet_grid(SM_CLUSTER ~ ., scales = "free_y", switch = "y") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 1, size = 14, face = "bold"),
        axis.text = element_blank(),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        # panel.spacing = unit(1.5, "lines"),
        strip.text.y = element_text(hjust = 0.5, size = 14, face = "bold", angle = 180),
        strip.background = element_rect(fill="white"),
        legend.text = element_text(size = 13),
        legend.position = "bottom",
        legend.title = element_text(size = 13, face = "bold"),
        plot.margin = unit(rep(0.5, 4), "cm"))



png(filename = "temp.png", width = 4000, height = 6000, res = 320)
pt
dev.off()















