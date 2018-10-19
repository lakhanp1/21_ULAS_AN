library(chipmine)
library(org.Anidulans.eg.db)

## this script plots the polII signal heatmap and polII fold change heatmap
## for the kdmB complex members

rm(list = ls())

path <- "E:/Chris_UM/Analysis/21_CL2017_ChIPmix_ULAS_MIX/ULAS_AN/kdmB_analysis/polII_signal"
setwd(path)


file_exptInfo <-"E:/Chris_UM/Analysis/21_CL2017_ChIPmix_ULAS_MIX/ULAS_AN/data/referenceData/sampleInfo.txt"
TF_dataPath <- "E:/Chris_UM/Analysis/21_CL2017_ChIPmix_ULAS_MIX/ULAS_AN/data/TF_data"
polII_dataPath <- "E:/Chris_UM/Analysis/21_CL2017_ChIPmix_ULAS_MIX/ULAS_AN/data/polII_data"
file_genes <- "E:/Chris_UM/Analysis/21_CL2017_ChIPmix_ULAS_MIX/ULAS_AN/data/referenceData/AN_genesForPolII.bed"

file_geneInfo <- "E:/Chris_UM/Database/A_Nidulans/A_nidulans_FGSC_A4_geneClasses.txt"

file_factorSignal <- paste(polII_dataPath, "/factors_polII_signal.tab", sep = "")
file_factorLfc <- paste(polII_dataPath, "/factors_polII_signal_LFC.tab", sep = "")

geneCdsFile <- "E:/Chris_UM/Database/A_Nidulans/A_nidulans_FGSC_A4_version_s10-m04-r03_CDS_Unique.bed"

orgDb <- org.Anidulans.eg.db

file_polIIsamples <- "kdmB_polII_sample.list"


##################################################################################
exptInfo <- readr::read_tsv(file = file_exptInfo, col_names = T, na = character())

polIISamples <- fread(file = file_polIIsamples, sep = "\t", header = F,
                      stringsAsFactors = F, col.names = c("id"), data.table = F) %>% 
  dplyr::left_join(y = exptInfo, by = c("id" = "sampleId"))


factorSignal <- readr::read_tsv(file = file_factorSignal, col_names = T) %>% 
  dplyr::select(sample, kdmB, rpdA, sntB, ecoA, sudA, laeA)

dt <- dplyr::left_join(x = polIISamples, y = factorSignal, by = c("id" = "sample")) %>% 
  dplyr::mutate(sample = gsub(pattern = "_[12]$", replacement = "", x = id, perl = T))


##################################################################################
## replicate 1 polII signal heatmap

rep1Df <- dplyr::filter(dt, rep == 1) %>% 
  dplyr::select(id, colnames(factorSignal)[-1])


rep1Mat <- log2(as.matrix(tibble::column_to_rownames(rep1Df, "id")) + 1)


quantile(rep1Mat, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)

rep1Col <- colorRamp2(breaks = quantile(rep1Mat, c(0, 0.1, 0.9, 0.99, 0.999)),
                    colors = c("white", "#f7f4f9", "#ce1256", "#980043", "#67001f"))


rep1Ht <- Heatmap(matrix = rep1Mat,
              col = rep1Col,
              column_title = "polII signal for deletion mutants of kdmB complex",
              column_title_gp = gpar(fontsize = 24, fontface = "bold"),
              # cell_fun = function(j, i, x, y, width, height, fill) {
              #   grid.text(sprintf("%.1f", rep1Mat[i, j]), x, y, gp = gpar(fontsize = 10))
              # },
              heatmap_legend_param = list(
                title = "log2(polII signal)",
                # at = colAt,
                legend_height = unit(4, "cm"),
                labels_gp = gpar(fontsize = 14),
                title_gp = gpar(fontsize = 18, fontface = "bold"),
                title_position = "topcenter"
              ),
              cluster_rows = FALSE, cluster_columns = FALSE,
              row_names_side = "left",
              column_names_gp = gpar(fontsize = 20, fontface = "bold"),
              row_names_gp = gpar(fontsize = 18),
              row_names_max_width = unit(20, "cm")
)


png(filename = "kdmB_rep1_member_signal_heatmap.png", width = 5000, height = 5000, res = 450)
draw(rep1Ht)
dev.off()

##################################################################################
## polII signal fold change heatmap

lfcDataAll <- readr::read_tsv(file = file_factorLfc, col_names = T) %>% 
  dplyr::select(lfcCol, sampleId, control, kdmB, rpdA, sntB, ecoA, sudA, laeA)


lfcDf <- dplyr::left_join(x = polIISamples, y = lfcDataAll, by = c("id" = "sampleId")) %>% 
  dplyr::filter(! is.na(lfcCol)) %>% 
  dplyr::select(!! colnames(lfcDataAll)[-2:-3])


lfcMat <- tibble::column_to_rownames(lfcDf, var = "lfcCol") %>% 
  as.matrix()

rownames(lfcMat) <- gsub(pattern = "(lfc.|An_|_polII)", replacement = "", x = rownames(lfcMat)) %>% 
  gsub(pattern = "_vs_", replacement = " / ", x = .)


lfc_color <- colorRamp2(breaks = c(-2, -1, -0.5, 0, 0.5, 1, 2),
                        # colors = RColorBrewer::brewer.pal(n = 7, name = "PuOr"),
                        colors = c("#b35806", "#f1a340", "#f7f7f7", "#f7f7f7", "#f7f7f7", "#998ec3", "#542788"))


htLfc <- Heatmap(matrix = lfcMat,
                 col = lfc_color,
                 column_title = "polII signal fold change for deletion mutants of kdmB complex",
                 column_title_gp = gpar(fontsize = 20, fontface = "bold"),
                 # cell_fun = function(j, i, x, y, width, height, fill) {
                 #   grid.text(sprintf("%.1f", lfcMat[i, j]), x, y, gp = gpar(fontsize = 10))
                 # },
                 heatmap_legend_param = list(
                   title = "log2(fold change)",
                   legend_height = unit(4, "cm"),
                   labels_gp = gpar(fontsize = 14),
                   title_gp = gpar(fontsize = 18, fontface = "bold"),
                   title_position = "topcenter"
                 ),
                 cluster_rows = FALSE, cluster_columns = FALSE,
                 row_names_side = "left",
                 column_names_gp = gpar(fontsize = 16, fontface = "bold"),
                 row_names_gp = gpar(fontsize = 14),
                 row_names_max_width = unit(10, "cm")
)


png(filename = "kdmB_member_lfc_heatmap.png", width = 4000, height = 4000, res = 380)
draw(htLfc)
dev.off()





