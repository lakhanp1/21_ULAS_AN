library(chipmine)
library(org.Anidulans.eg.db)


## this script performs following tasks
## 1) extract the polII signal for all samples
## 2) prepare polII signal matrix for gene of interest (kdmB complex members) and plot heatmap
## 3) calculate the log2 fold change for each polII sample using appropriate control data


rm(list = ls())

path <- "E:/Chris_UM/Analysis/21_CL2017_ChIPmix_ULAS_MIX/ULAS_AN/data"
setwd(path)


file_exptInfo <-"E:/Chris_UM/Analysis/21_CL2017_ChIPmix_ULAS_MIX/ULAS_AN/data/referenceData/sampleInfo.txt"
TF_dataPath <- "E:/Chris_UM/Analysis/21_CL2017_ChIPmix_ULAS_MIX/ULAS_AN/data/TF_data"
polII_dataPath <- "E:/Chris_UM/Analysis/21_CL2017_ChIPmix_ULAS_MIX/ULAS_AN/data/polII_data"
file_genes <- "E:/Chris_UM/Analysis/21_CL2017_ChIPmix_ULAS_MIX/ULAS_AN/data/referenceData/AN_genesForPolII.bed"

file_geneInfo <- "E:/Chris_UM/Database/A_Nidulans/A_nidulans_FGSC_A4_geneClasses.txt"

file_factors <- "E:/Chris_UM/Analysis/21_CL2017_ChIPmix_ULAS_MIX/ULAS_AN/data/referenceData/factor_names.list"

orgDb <- org.Anidulans.eg.db

file_polIIsamples <- paste(polII_dataPath, "/polII_sample.list", sep = "")
file_polIICtrlPairs <- paste(polII_dataPath, "/polII_sample_control_pairs.txt", sep = "")


##################################################################################
## extract polII signal matrix for all the genes
geneSet <- data.table::fread(file = file_genes, header = F,
                             col.names = c("chr", "start", "end", "gene", "score", "strand")) %>% 
  dplyr::select(-score) %>% 
  dplyr::mutate(length = end - start)


geneDesc <- select(x = orgDb, keys = geneSet$gene, columns = "DESCRIPTION", keytype = "GID")

geneSet <- dplyr::left_join(x = geneSet, y = geneDesc, by = c("gene" = "GID"))


polIISamples <- fread(file = file_polIIsamples, sep = "\t", header = F,
                      stringsAsFactors = F, col.names = c("id"), data.table = F)


polII_info <- get_sample_information(exptInfoFile = file_exptInfo,
                                     samples = polIISamples$id,
                                     dataPath = polII_dataPath)

polIICols <- list(
  exp = structure(polIISamples$id, names = polIISamples$id),
  is_expressed = structure(paste("is_expressed", ".", polIISamples$id, sep = ""), names = polIISamples$id)
)



polIIMat <- get_polII_expressions(genesDf = geneSet, exptInfo = polII_info) %>% 
  # dplyr::select(-starts_with("is_expressed")) %>%
  dplyr::select(chr, start, end, gene, strand, length, DESCRIPTION, everything())


fwrite(x = polIIMat, file = paste(polII_dataPath, "/polII_signal_matrix.tab", sep = ""),
       sep = "\t", col.names = T, quote = F, row.names = F)

##################################################################################
## polII signal percentile matrix
quantile(x = polIIMat$An_untagged_20h_polII_1,
         c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.997, 0.999, 0.9999, 1), na.rm = T)


polIIQuantiles <- lapply(
  X = polIICols$exp,
  FUN = function(x){
    quantile(
      x = polIIMat[[x]],
      c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.997, 0.999, 0.9999, 1), na.rm = T)
  }) %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "quantile") %>% 
  # dplyr::mutate(quantile = paste("quantile_", quantile, sep = "")) %>% 
  tidyr::gather(key = sample, value = signal, -quantile) %>% 
  tidyr::spread(key = quantile, value = signal)


fwrite(x = polIIQuantiles, file = paste(polII_dataPath, "/polII_signal_quantiles.tab", sep = ""),
       sep = "\t", col.names = T, quote = F, row.names = F)


##################################################################################
## polII signal matrix for genes of interest

goi <- read_tsv(file = file_factors, col_names = T)

goiPolII <- dplyr::left_join(goi, polIIMat, by = c("id" = "gene")) 


exprDf <- dplyr::select(goiPolII, -id, -chr, -start, -end, -strand, -length, -DESCRIPTION,
                        -starts_with("is_expressed")) %>% 
  tidyr::gather(key = sample, value = expression, -name, factor_key = TRUE) %>% 
  tidyr::spread(key = name, value = expression) %>% 
  dplyr::mutate_if(.predicate = is.factor, .funs = as.character) %>% 
  dplyr::select(sample, goi$name) %>% 
  as.data.frame()


fwrite(x = exprDf, file = paste(polII_dataPath, "/factors_polII_signal.tab", sep = ""),
       sep = "\t", col.names = T, quote = F, row.names = F)



exprMat <- log2(as.matrix(tibble::column_to_rownames(exprDf, "sample")) + 1)


quantile(exprMat, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)

htCol <- colorRamp2(breaks = quantile(exprMat, c(0, 0.1, 0.9, 0.99, 0.999)),
                    colors = c("white", "#f7f4f9", "#ce1256", "#980043", "#67001f"))

colAt <- as.numeric(sprintf("%.0f", c(seq(quantile(exprMat, 0.1), quantile(exprMat, 0.99), length = 4),
                                      quantile(exprMat, 0.995))))

ht <- Heatmap(matrix = exprMat,
              col = htCol,
              column_title = "polII signal for various factors in mutants",
              column_title_gp = gpar(fontsize = 20, fontface = "bold"),
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf("%.1f", exprMat[i, j]), x, y, gp = gpar(fontsize = 10))
              },
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
              column_names_gp = gpar(fontsize = 16, fontface = "bold"),
              row_names_gp = gpar(fontsize = 14),
              row_names_max_width = unit(10, "cm")
)


png(filename = paste(polII_dataPath, "/factors_polII_signal.png", sep = ""), width = 4000, height = 5000, res = 400)
draw(ht)
dev.off()


##################################################################################
## polII data fold change results

polIICtrlPairs <- readr::read_tsv(file = file_polIICtrlPairs, col_names = T) %>% 
  as.data.frame() %>% 
  dplyr::mutate(lfcCol = paste("lfc.", sampleId, "_vs_", control, sep = ""))

lfcMat <- polIIMat

# i <- 1

for (i in 1:nrow(polIICtrlPairs)) {
  
  lfcMat <- get_fold_change(df = lfcMat,
                            nmt = polIICtrlPairs[i, 1],
                            dmt = polIICtrlPairs[i, 2],
                            newCol = polIICtrlPairs[i, 3],
                            isExpressedCols = polIICols$is_expressed)
  
}


lfcMat <- dplyr::select(lfcMat, gene, chr, start, end, strand, length, DESCRIPTION, starts_with("lfc."))


fwrite(x = lfcMat, file = paste(polII_dataPath, "/polII_LFC_matrix.tab", sep = ""),
       sep = "\t", col.names = T, quote = F, row.names = F)


##################################################################################
## fold change heatmap for genes of interest
goi <- read_tsv(file = file_factors, col_names = T)

goiLfc <- dplyr::left_join(goi, lfcMat, by = c("id" = "gene")) 


lfcDf <- dplyr::select(goiLfc, name, starts_with("lfc.")) %>% 
  tidyr::gather(key = lfcCol, value = expression, -name, factor_key = TRUE) %>% 
  tidyr::spread(key = name, value = expression) %>% 
  dplyr::mutate_if(.predicate = is.factor, .funs = as.character) %>% 
  dplyr::left_join(y = polIICtrlPairs, by = c("lfcCol" = "lfcCol")) %>% 
  dplyr::select(lfcCol, sampleId, control, goi$name) %>% 
  as.data.frame()


fwrite(x = lfcDf, file = paste(polII_dataPath, "/factors_polII_signal_LFC.tab", sep = ""),
       sep = "\t", col.names = T, quote = F, row.names = F)


lfcMat <- dplyr::select(lfcDf, -sampleId, -control) %>% 
  tibble::column_to_rownames(var = "lfcCol") %>% 
  as.matrix()

rownames(lfcMat) <- gsub(pattern = "(lfc.|An_|_polII_\\d)", replacement = "", x = rownames(lfcMat)) %>% 
  gsub(pattern = "_vs_", replacement = " / ", x = .)

lfc_color <- colorRamp2(breaks = c(-2, -1, -0.5, 0, 0.5, 1, 2),
                        # colors = RColorBrewer::brewer.pal(n = 7, name = "PuOr"),
                        colors = c("#b35806", "#f1a340", "#f7f7f7", "#f7f7f7", "#f7f7f7", "#998ec3", "#542788"))


htLfc <- Heatmap(matrix = lfcMat,
                 col = lfc_color,
                 column_title = "polII signal fold change for various factors in mutants",
                 column_title_gp = gpar(fontsize = 20, fontface = "bold"),
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


png(filename = paste(polII_dataPath, "/factors_polII_signal_LFC.png", sep = ""), width = 4000, height = 5000, res = 400)
draw(htLfc)
dev.off()




