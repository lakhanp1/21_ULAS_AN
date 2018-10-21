library(chipmine)
library(org.Anidulans.eg.db)
library(scales)
library(ggplot2)


## This script compares the binding of two TF


rm(list = ls())

source(file = "E:/Chris_UM/Codes/GO_enrichment/topGO_functions.R")

# path <- "E:/Chris_UM/Analysis/21_CL2017_ChIPmix_ULAS_MIX/ULAS_AN/kdmB_analysis/kdmB_48h_vs_20h"
# setwd(path)

##################################################################################

## IMP: the first sampleID will be treated primary and clustering will be done/used for/of this sample
comparisonName <- "kdmB_48h_vs_20h"
outPrefix <- here::here("kdmB_analysis/kdmB_48h_vs_20h", comparisonName)

tf1 <- "An_kdmB_20h_HA_1"
tf2 <- "An_kdmB_48h_HA_1"
polII1 <- "An_untagged_20h_polII_1"
polII2 <- "An_untagged_48h_polII_1"

otherTfs <- c("An_untagged_20h_HA_1", "An_untagged_48h_HA_1")
# "An_kdmB_20h_HA_2", "An_kdmB_48h_HA_2", 
otherPolII <- c("An_kdmB_del_20h_polII_1", "An_kdmB_del_48h_polII_1")
otherHist <- c("An_H3_20h_HIST_1", "An_H3_48h_HIST_1")

# "deeptools", "miao", "normalizedmatrix", "normalizedmatrix_5kb"
matrixType <- "normalizedmatrix_5kb"
matrixDim = c(500, 200, 100, 10)


## polII signal fold change pairs
polIIDiffPairs <- list(
  p1 = list(
    name = "untagged_48h_vs_untagged_20h_polII",
    title = "polII log2(untagged_48h \n vs untagged_20h)",
    samples = c(polII2, polII1)
  ),
  p2 = list(
    name = "kdmB_del_48h_vs_kdmB_del_20h_polII",
    title = "polII log2(kdmB_del_48h \n vs kdmB_del_20h)",
    samples = otherPolII[c(2, 1)]
  ),
  p3 = list(
    name = "kdmB_del_20h_vs_untagged_20h_polII",
    title = "polII log2(kdmB_del_20h \n vs untagged_20h_polII)",
    samples = c(otherPolII[1], polII1)
  ),
  p4 = list(
    name = "kdmB_del_48h_vs_untagged_48h_polII",
    title = "polII log2(kdmB_del_48h \n vs untagged_48h_polII)",
    samples = c(otherPolII[2], polII2)
  )
)


sampleList <- c(tf1, tf2, polII1, polII2, otherTfs, otherPolII, otherHist)



clusterStorePath <- paste(outPrefix, "_profile.kmeans.clusters.txt", sep = "")

## genes to read
file_exptInfo <- here::here("data", "referenceData/sampleInfo.txt")
file_genes <- here::here("data", "referenceData/AN_genesForPolII.bed")
file_topGoMap <- "E:/Chris_UM/Database/A_Nidulans/ANidulans_OrgDb/geneid2go.ANidulans.topGO.map"
file_geneInfo <- "E:/Chris_UM/Database/A_Nidulans/A_nidulans_FGSC_A4_geneClasses.txt"


TF_dataPath <- here::here("data", "TF_data")
polII_dataPath <- here::here("data", "polII_data")
hist_dataPath <- here::here("data", "histone_data")


orgDb <- org.Anidulans.eg.db

showExpressionHeatmap <- FALSE


outPrefix_tfPolII <- paste(outPrefix, "_tfPolII", sep = "")
outPrefix_tfSpecific <- paste(outPrefix, "_specific_binding", sep = "")
outPrefix_peakExp <- paste(outPrefix, "_pkExpGenes", sep = "")


anLables <- list()
anLables[["is_SM_gene"]] <- "SM gene"
anLables[["is_TF"]] <- "Transcription Factor"
anLables[["gene_length"]] <- "Gene Length"


##################################################################################

## genes to read
geneSet <- data.table::fread(file = file_genes, header = F,
                             col.names = c("chr", "start", "end", "gene", "score", "strand")) %>% 
  dplyr::mutate(length = end - start)

geneDesc <- AnnotationDbi::select(x = orgDb, keys = geneSet$gene, columns = "DESCRIPTION", keytype = "GID")

geneSet <- dplyr::left_join(x = geneSet, y = geneDesc, by = c("gene" = "GID"))

## gene information annotations: cluster and TF and polII expression values
geneInfo <- add_gene_info(file = file_geneInfo, clusterDf = geneSet)

head(geneInfo)

##################################################################################

## read the experiment sample details and select only those which are to be plotted

tfData <- get_sample_information(exptInfoFile = file_exptInfo,
                                 samples = c(tf1, tf2, otherTfs),
                                 dataPath = TF_dataPath,
                                 matrixSource = matrixType)

polIIData <- get_sample_information(exptInfoFile = file_exptInfo,
                                    samples = c(polII1, polII2, otherPolII),
                                    dataPath = polII_dataPath,
                                    matrixSource = matrixType)

histData <- get_sample_information(exptInfoFile = file_exptInfo,
                                   samples = otherHist,
                                   dataPath = hist_dataPath,
                                   matrixSource = matrixType)

exptData <- dplyr::bind_rows(tfData, polIIData, histData)


dfToList <- function(x, n){
  return(structure(x, names = n))
}


sampleList <- lapply(exptData, dfToList, n = exptData$sampleId)


polII_ids <- exptData$sampleId[which(exptData$IP_tag == "polII")]
tfIds <- exptData$sampleId[which(exptData$IP_tag %in% c("HA", "MYC", "TAP") & exptData$TF != "untagged")]


polIICols <- list(
  exp = structure(polII_ids, names = polII_ids),
  is_expressed = structure(paste("is_expressed", ".", polII_ids, sep = ""), names = polII_ids)
)


tfCols <- sapply(c("hasPeak", "pval", "peakType", "tesPeakType", "peakDist", "summitDist", "upstreamExpr", "peakExpr", "relativeDist"),
                 FUN = function(x){ structure(paste(x, ".", tfIds, sep = ""), names = tfIds) },
                 simplify = F, USE.NAMES = T)


expressionData <- get_polII_expressions(exptInfo = exptData,
                                        genesDf = geneInfo)

## add fold change columns
for (i in names(polIIDiffPairs)) {
  expressionData <- get_fold_change(df = expressionData,
                                    s1 = polIIDiffPairs[[i]]$samples[1],
                                    s2 = polIIDiffPairs[[i]]$samples[2],
                                    newCol = polIIDiffPairs[[i]]$name,
                                    isExpressedCols = polIICols$is_expressed)
}


expressionData <- get_TF_binding_data(exptInfo = tfData,
                                      genesDf = expressionData)


expressionData %>% 
  dplyr::select(gene, starts_with("hasPeak")) %>% 
  dplyr::group_by_at(.vars = vars(starts_with("hasPeak"))) %>% 
  dplyr::summarise(n = n())


expressionData$group <- dplyr::group_by_at(expressionData, unname(tfCols$hasPeak[c(tf1, tf2)])) %>%
  dplyr::group_indices()

rownames(expressionData) = expressionData$gene

fwrite(x = expressionData, file = paste(outPrefix, "_data.tab", sep = ""),
       sep = "\t", col.names = T, quote = F, na = "NA")

## peak data
hasPeakDf <- dplyr::filter_at(.tbl = expressionData,
                              .vars = vars(unname(tfCols$hasPeak[c(tf1, tf2)])),
                              .vars_predicate = any_vars(. == "TRUE"))

rownames(hasPeakDf) <- hasPeakDf$gene


dplyr::group_by_at(.tbl = hasPeakDf, .vars = vars(unname(tfCols$hasPeak[c(tf1, tf2)]), group)) %>%
  dplyr::summarise(n = n())

## check the MA plot
plot_MA_gg(df = expressionData, s1 = polII2, s2 = polII1, title = "MA plot", colorCol = "group")
# plot_MA_gg(df = expressionData, s1 = polII1, s2 = polII2, title = "MA plot", colorCol = "group")
plot_MA_gg(df = expressionData, s1 = otherPolII[2], s2 = otherPolII[1], title = "MA plot", colorCol = "group")
# plot_MA_gg(df = expressionData, s1 = otherPolII[1], s2 = otherPolII[2], title = "MA plot", colorCol = "group")


plot_scatter(df = expressionData, s1 = polII1, s2 = polII2, title = "Scatter plot", colorCol = "group")
plot_scatter(df = expressionData, s1 = polII2, s2 = polII1, title = "Scatter plot", colorCol = "group")


##################################################################################

## join the profile matrices and do clustering
# mergedKm = merged_profile_matrix_cluster(name = comparisonName,
#                                          exptInfo = exptData,
#                                          genes = geneSet$gene,
#                                          clusterStorePath = clusterStorePath,
#                                          k = 12)

# clusterInfo <- fread(file = clusterStorePath, sep = "\t", header = T, stringsAsFactors = F)
# 
# expressionData <- dplyr::left_join(x = expressionData, y = clusterInfo, by = c("gene" = "gene"))




##################################################################################
## topGO enrichment
goEnrich <- dplyr::group_by_at(.tbl = hasPeakDf, .vars = vars(unname(tfCols$hasPeak[c(tf1, tf2)]), group)) %>%
  do(topGO_enrichment(goMapFile = file_topGoMap, genesOfInterest = .$gene, goNodeSize = 5))


fwrite(x = goEnrich,
       file = paste(outPrefix, "_GO.tab", sep = ""), sep = "\t", col.names = T, quote = F)


## clusterProfiler groupGO assignment
grpGo <- dplyr::group_by_at(.tbl = hasPeakDf, .vars = vars(unname(tfCols$hasPeak[c(tf1, tf2)]), group)) %>% 
  do(clusterProfiler_groupGO(genes = .$gene, org = orgDb, goLevel = 3, type = "BP", keyType ="GID"))

fwrite(x = grpGo,
       file = paste(outPrefix, "_groupGO.tab", sep = ""), sep = "\t", col.names = T, quote = F)


## pathway enrichment
keggEnr <- dplyr::group_by_at(.tbl = hasPeakDf, .vars = vars(unname(tfCols$hasPeak[c(tf1, tf2)]), group)) %>% 
  do(keggprofile_enrichment(genes = .$gene, orgdb = orgDb, keytype = "GID", keggOrg = "ani", pvalCut = 0.05))


fwrite(x = keggEnr,
       file = paste(outPrefix, "_KEGG.tab", sep = ""), sep = "\t", col.names = T, quote = F)


##################################################################################
## GO enrichment and KEGG testing
gl <- hasPeakDf$gene[ which( !hasPeakDf[[tfCols$hasPeak[tf1]]] & hasPeakDf[[tfCols$hasPeak[tf2]]] ) ]


## KEGG pathway enrichment using enrichKEGG
keggIds <- AnnotationDbi::select(x = orgDb, keys = gl, columns = c("GID", "KEGG_ID"), keytype = "GID") %>% 
  dplyr::filter(!is.na(KEGG_ID))

kk <- clusterProfiler::enrichKEGG(gene = keggIds$KEGG_ID, organism = 'ani', keyType = "kegg",
                                  pvalueCutoff = 0.05, pAdjustMethod = "none", qvalueCutoff = 1,
                                  minGSSize = 1)


kkRes <- dplyr::filter(kk@result, pvalue <= 0.05)

## using KEGGprofile package
kp <- KEGGprofile::find_enriched_pathway(gene = keggIds$KEGG_ID, species = 'ani',
                                         returned_pvalue = 0.05, returned_adjpvalue = 1,
                                         returned_genenumber = 1, download_latest = TRUE)

kp$stastic

##################################################################################
## clusterProfiler enrichment: testing

# ## groupGO
# grpGo <- clusterProfiler::groupGO(gene = gl,
#                                   OrgDb = orgDb,
#                                   keyType = "GID", ont = "BP", level = 3)
# 
# grpGoDf <- dplyr::filter(grpGo@result, Count > 0)
# 
# grpGoDf[, 1:4]
# 
# ## enricher
# t2g <- AnnotationDbi::select(x = orgDb, keys = keys(orgDb), columns = c("GO", "GID","ONTOLOGY"), keytype = "GID") %>% 
#   dplyr::filter(ONTOLOGY == "BP") %>% 
#   dplyr::select(GO, GID)
# 
# t2n <- AnnotationDbi::select(x = GO.db, keys = unique(t2g$GO), columns = c("GOID", "TERM"), keytype = "GOID")
# 
# enrGo <- clusterProfiler::enricher(gene = gl, TERM2GENE = t2g, TERM2NAME = t2n)
# enrGoDf <- dplyr::filter(enrGo@result, pvalue <= 0.05)
# 
# enrGoDf[, 1:7]
# 
# 
# ## enricher
# enrGo2 <- clusterProfiler::enrichGO(gene = gl, OrgDb = orgDb, keyType = "GID", ont = "BP",
#                                     pvalueCutoff = 0.05)
# 
# enrGo2Df <- dplyr::filter(enrGo2@result, pvalue <= 0.05)
# enrGo2Df[, 1:7]
# 
# write.table(enrGo2Df, "clipboard", sep = "\t", quote = F)
# 
# 
# ## formula based enrichment
# keytypes(orgDb)
# reformulate(termlabels = c(tfCols$hasPeak, "group"), response = "gene")
# 
# formulaRes <- clusterProfiler::compareCluster(
#   geneClusters = gene ~ group,
#   fun = "enrichGO",
#   data = hasPeakDf,
#   OrgDb = orgDb,
#   keyType = "GID",
#   ont = "BP", minGSSize = 2)

# formSimple <- clusterProfiler::simplify(formulaRes, cutoff = 0.7)
# 
# plt <- dotplot(formSimple, showCategory = NULL) +
#   scale_y_discrete(labels = wrap_format(60)) +
#   ggtitle(paste(titleName, ": GO enrichment for DEGs")) +
#   theme_bw() +
#   theme(plot.title = element_text(hjust = 0.8, size = 16, face = "bold"),
#         axis.text.x = element_text(size = 14),
#         axis.text.y = element_text(size = 14),
#         axis.title = element_text(face = "bold"),
#         panel.grid = element_blank(),
#         legend.text = element_text(size = 14),
#         legend.title = element_text(size = 14, face = "bold"))
# 

##################################################################################

## polII signal matrix
polIIMat <- data.matrix(log2(expressionData[, polII_ids] + 1))

quantile(polIIMat, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.997, 0.999, 0.9999, 1), na.rm = T)


polII_color <- colorRamp2(breaks = c(0, quantile(polIIMat, c(0.5, 0.8, 0.9, 0.92, 0.95, 0.97, 0.99, 0.995, 0.999))),
                          colors = c("white", RColorBrewer::brewer.pal(n = 9, name = "RdPu")))


## polII signal fold change matrix
lfcMat <- as.matrix(expressionData[, purrr::map_chr(polIIDiffPairs, "name"), drop = FALSE])

lfc_color <- colorRamp2(breaks = c(-2, -1, -0.5, 0, 0.5, 1, 2),
                        colors = RColorBrewer::brewer.pal(n = 7, name = "PuOr"))

##################################################################################

# ## peak pval heatmap
# quantile(hasPeakDf[, tfCols$pval], c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)
# 
# peakPvalCol <- colorRamp2(breaks = c(0, 1, quantile(hasPeakDf[, tfCols$pval], c(0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 0.9, 0.95), na.rm = T)),
#                           colors = c("grey", "white", RColorBrewer::brewer.pal(n = 9, name = "YlOrRd")))
# 
# 
# macs2Ht <- signal_heatmap(log2_matrix = hasPeakDf[, tfCols$pval],
#                           htName = "macs2_peak",
#                           col_title = tfCols$pval,
#                           legend_title = "log10(macs2_peak)",
#                           color = peakPvalCol,
#                           htSplit = hasPeakDf["group"],
#                           cluster_rows = TRUE, cluster_columns = FALSE)
# 
# 
# polIIMat_peaks <- polIIMat[hasPeakDf$gene , ]
# 
# polIIht_peaks <- signal_heatmap(log2_matrix = polIIMat,
#                                 htName = "polII_exp",
#                                 col_title = polII_ids,
#                                 legend_title = "log2(polII_singal)",
#                                 color = polII_color)

##################################################################################
## colors for profile matrix
matList <- profile_matrix_list(exptInfo = exptData,
                               geneList = geneInfo$gene,
                               source = matrixType,
                               up = matrixDim[1], target = matrixDim[2], down = matrixDim[3])



## tf colors
tfMeanProfile <- getSignalsFromList(lt = matList[tfData$sampleId])
quantile(tfMeanProfile, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)
tfMeanColor <- colorRamp2(quantile(tfMeanProfile, c(0.50, 0.995), na.rm = T), c("white", "red"))

tfColorList <- sapply(X = tfData$sampleId, FUN = function(x){return(tfMeanColor)})
tfWiseColors <- sapply(X = tfData$sampleId,
                       FUN = function(x){
                         cat(x, "\n")
                         print(quantile(
                           matList[[x]],
                           c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T))
                         
                         return(
                           colorRamp2(quantile(matList[[x]], c(0.50, 0.99), na.rm = T), c("white", "red"))
                         )
                       }
)

## histone colors
histMeanProfile <- getSignalsFromList(lt = matList[histData$sampleId])
quantile(histMeanProfile, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)
histMeanColor <- colorRamp2(quantile(histMeanProfile, c(0.20, 0.995), na.rm = T), c("black", "yellow"))

histColorList <- sapply(X = histData$sampleId, FUN = function(x){return(histMeanColor)})

## LFC(TF2/TF1) matrix
## adding Tf1's median to TF2 matrix and vice a versa. Doing this will remove the background skewness of fold change
## matrix and the median will be around 0
tfLfcMat <- log2((matList[[tf2]] + quantile(matList[[tf1]], 0.5)) / (matList[[tf1]] + quantile(matList[[tf2]], 0.5)))

quantile(log2((matList[[tf2]] + quantile(matList[[tf2]], 0.5)) / (matList[[tf1]] + quantile(matList[[tf1]], 0.5))),
         c(0, 0.01, 0.05, seq(0.1, 0.9, by = 0.1), 0.95, 0.99, 1))

plot(density(log2((matList[[tf2]] + quantile(matList[[tf1]], 0.5)) / (matList[[tf1]] + quantile(matList[[tf2]], 0.5))),
             na.rm = TRUE)
)

lfcProfileCol <- colorRamp2(breaks = c(-2, -1.5, -1, -0.75, 0, 0.75, 1, 1.5, 2),
                            colors = RColorBrewer::brewer.pal(n = 9, name = "PuOr"))

## scalled LFC matrix
tfScalledLfcMat <- scale(x = tfLfcMat, center = TRUE, scale = TRUE)

quantile(tfScalledLfcMat,
         c(0, 0.005, 0.01, 0.05, seq(0.1, 0.9, by = 0.1), 0.95, 0.99, 0.995, 1))

## TF1 scalled matrix
tf1ScalledMat <- scale(x = matList[[tf1]], center = TRUE, scale = TRUE)
quantile(tf1ScalledMat, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)

## TF2 scalled matrix
tf2ScalledMat <- scale(x = matList[[tf2]], center = TRUE, scale = TRUE)
quantile(tf2ScalledMat, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)

## Difference between TF2 and TF1 scalled matrix
scalledTfDiffMat <- tf2ScalledMat - tf1ScalledMat
plot(density(scalledTfDiffMat))
quantile(scalledTfDiffMat, c(0, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.02, 0.05, seq(0.1, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)


scalledTfDiffColor <- colorRamp2(breaks = c(-3, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 3),
                                 colors = RColorBrewer::brewer.pal(n = 11, name = "PuOr"))

##################################################################################
## profile plot with genes showing peak specific and common between the conditions
# tfSpecificDf <- dplyr::filter(hasPeakDf,
#                               !! as.name(tfCols$hasPeak[tf1]) == FALSE | !! as.name(tfCols$hasPeak[tf2]) == FALSE )

tfSpecificDf <- hasPeakDf

rownames(tfSpecificDf) <- tfSpecificDf$gene

clusters_tfSpecific <- dplyr::select(tfSpecificDf, gene, group) %>% 
  dplyr::rename(cluster = group)

ylimList <- sapply(X = tfData$sampleId, FUN = function(x){c(0, 30)}, USE.NAMES = T, simplify = F)

## TF profile plot
multiProf_tfSpecific <- multi_profile_plots(exptInfo = tfData,
                                            genesToPlot = tfSpecificDf$gene,
                                            matSource = matrixType,
                                            matBins = matrixDim,
                                            profileColors = tfColorList,
                                            clusters = clusters_tfSpecific,
                                            column_title_gp = gpar(fontsize = 12),
                                            ylimFraction = ylimList)


## Scalled TF diff profile
scalledTfDiffProf_tfSpecific <- profile_heatmap(profileMat = scalledTfDiffMat[tfSpecificDf$gene, ],
                                                signalName = comparisonName,
                                                profileColor = scalledTfDiffColor,
                                                geneGroups = clusters_tfSpecific)

## histone profile plots
histProfiles <- multi_profile_plots(exptInfo = histData,
                                    genesToPlot = tfSpecificDf$gene,
                                    matSource = matrixType,
                                    matBins = matrixDim,
                                    profileColors = histColorList,
                                    clusters = clusters_tfSpecific,
                                    column_title_gp = gpar(fontsize = 12),
                                    drawClusterAn = FALSE)

## polII signal heatmap
polIIMat_tfSpecific <- polIIMat[tfSpecificDf$gene, ]

polIIHt_tfSpecific <- signal_heatmap(log2_matrix = polIIMat_tfSpecific,
                                     htName = "polII_exp",
                                     col_title = polII_ids,
                                     legend_title = "log2(polII_singal)",
                                     color = polII_color,
                                     cluster_columns = FALSE)

## polII signal fold change heatmap
lfc_tfSpecific <- lfcMat[tfSpecificDf$gene, ]

lfc_heatmap <- signal_heatmap(log2_matrix = lfc_tfSpecific,
                              htName = "polII_lfc",
                              col_title = purrr::map_chr(polIIDiffPairs, "title"),
                              legend_title = "log2(fold change)",
                              color = lfc_color,
                              cluster_columns = FALSE)


gl_tfSpecific <- gene_length_heatmap_annotation(bedFile = file_genes, genes = tfSpecificDf$gene)

htlist_tfSpecific <- gl_tfSpecific$an + 
  multiProf_tfSpecific$heatmapList + 
  scalledTfDiffProf_tfSpecific$heatmap +
  polIIHt_tfSpecific + 
  lfc_heatmap
# histProfiles$heatmapList




if( all(rownames(multiProf_tfSpecific$profileHeatmaps[[1]]$heatmap@matrix) == tfSpecificDf$gene) ){
  rowOrd <- order(tfSpecificDf[[ tfCols$peakDist[tf1] ]], tfSpecificDf[[ tfCols$peakDist[tf2] ]], decreasing = TRUE)
}


title_tfSpecific <- "Differential binding of kdmB at 20h and 48h"


# draw Heatmap and add the annotation name decoration
png(filename = paste0(outPrefix_tfSpecific, "_profile.png", collapse = ""), width=5000, height=3500, res = 250)

draw(htlist_tfSpecific,
     main_heatmap = tfData$profileName[1],
     # annotation_legend_list = list(profile1$legend),
     column_title = title_tfSpecific,
     column_title_gp = gpar(fontsize = 14, fontface = "bold"),
     row_sub_title_side = "left",
     gap = unit(6, "mm"),
     row_order = rowOrd,
     padding = unit(rep(0.5, times = 4), "cm")
)

## decorate the annotations
add_annotation_titles(annotations = c("gene_length"), anTitle = anLables, fontSize = 12)

row_annotation_axis(an = "gene_length",
                    at = c(0, 2000, 4000),
                    labels = c("0kb", "2kb", ">4kb"),
                    slice = length(unique(tfSpecificDf$group)))

dev.off()


##################################################################################
## genes which have peak in both TFs and polII signal in either one or both TFs

peakExpDf <- dplyr::filter_at(.tbl = expressionData,
                              .vars = unname(tfCols$hasPeak[c(tf1, tf2)]),
                              .vars_predicate = all_vars(. == TRUE)) %>% 
  dplyr::filter_at(.vars = unname(polIICols$is_expressed[c(polII1, polII2)]),
                   .vars_predicate = any_vars(. == TRUE)) 

## group genes into polII fold change bins
peakExpDf <- dplyr::mutate(peakExpDf,
                           group = case_when(
                             !! as.name(polIIDiffPairs$p1$name) >= 2 ~ "G1: LFC >= 2",
                             !! as.name(polIIDiffPairs$p1$name) >= 1 ~ "G2: 1 <= LFC < 2",
                             !! as.name(polIIDiffPairs$p1$name) >= 0.5 ~ "G3: 0.5 <= LFC < 1",
                             !! as.name(polIIDiffPairs$p1$name) >= 0 ~ "G4: 0 <= LFC < 0.5",
                             !! as.name(polIIDiffPairs$p1$name) <= -2 ~ "G8: -2 >= LFC",
                             !! as.name(polIIDiffPairs$p1$name) <= -1 ~ "G7: -2 < LFC <= -1",
                             !! as.name(polIIDiffPairs$p1$name) <= -0.5 ~ "G6: -1 < LFC <= -0.5",
                             !! as.name(polIIDiffPairs$p1$name) < 0 ~ "G5: -0.5 < LFC < 0",
                             TRUE ~ "0"
                           ))

dplyr::group_by(peakExpDf, group) %>% 
  dplyr::summarise(n = n())


peakExp_clusters <- dplyr::select(peakExpDf, gene, group) %>% 
  dplyr::rename(cluster = group)


## profile heatmap: no need to include input
multiProf_peakExp <- multi_profile_plots(exptInfo = tfData[tfData$sampleId %in% tfIds, ],
                                         genesToPlot = peakExpDf$gene,
                                         matSource = matrixType,
                                         matBins = matrixDim,
                                         clusters = peakExp_clusters,
                                         profileColors = tfWiseColors,
                                         # ylimFraction = ylimList,
                                         column_title_gp = gpar(fontsize = 14))


## Scalled TF diff profile
scalledTfDiffProf_peakExp <- profile_heatmap(profileMat = scalledTfDiffMat[peakExpDf$gene, ],
                                             signalName = comparisonName,
                                             profileColor = scalledTfDiffColor,
                                             geneGroups = peakExp_clusters,
                                             ylimFraction = c(-2.3, 1.7))
## polII signal heatmap
polIIMat_peakExp <- polIIMat[peakExpDf$gene, ]

polIIht_peakExp <- signal_heatmap(log2_matrix = polIIMat_peakExp,
                                  htName = "polII_exp",
                                  col_title = polII_ids,
                                  legend_title = "log2(polII_singal)",
                                  color = polII_color,
                                  cluster_columns = FALSE)

## polII fold change heatmap
lfc_peakExp <- lfcMat[peakExpDf$gene, ]

lfcHt_peakExp <- signal_heatmap(log2_matrix = lfc_peakExp,
                                htName = "polII_lfc",
                                col_title = purrr::map_chr(polIIDiffPairs, "title"),
                                legend_title = "log2(fold change)",
                                color = lfc_color,
                                cluster_columns = FALSE)

## gene length annotation
gl_peakExp <- gene_length_heatmap_annotation(bedFile = file_genes, genes = peakExpDf$gene)

htlist_peakExp <- gl_peakExp$an +
  multiProf_peakExp$heatmapList +
  scalledTfDiffProf_peakExp$heatmap +
  polIIht_peakExp +
  lfcHt_peakExp



if( all(rownames(multiProf_peakExp$profileHeatmaps[[1]]$heatmap@matrix) == peakExpDf$gene) ){
  ## row order by polII LFC
  rowOrd <- order(peakExpDf[[ polIIDiffPairs$p1$name ]], decreasing = TRUE)
}


title_peakExp <- "kdmB binding vs polII transcription: genes bound both at 20h and 48h and atleast one timepoint showing polII signal "


# draw Heatmap and add the annotation name decoration
png(filename = paste0(outPrefix_peakExp, "_profile.png", collapse = ""), width=4500, height=3500, res = 270)

draw(htlist_peakExp,
     main_heatmap = tfData$profileName[1],
     # annotation_legend_list = list(profile1$legend),
     column_title = title_peakExp,
     column_title_gp = gpar(fontsize = 14, fontface = "bold"),
     row_sub_title_side = "left",
     gap = unit(6, "mm"),
     row_order = rowOrd,
     padding = unit(rep(0.5, times = 4), "cm")
)

## decorate the annotations
add_annotation_titles(annotations = c("gene_length"), anTitle = anLables, fontSize = 12)

row_annotation_axis(an = "gene_length",
                    at = c(0, 2000, 4000),
                    labels = c("0kb", "2kb", ">4kb"),
                    slice = length(unique(peakExp_clusters$cluster)))

dev.off()


##################################################################################
## genes with either peak detected or with top 10% polII signal 

expressionData %>% 
  dplyr::group_by_at(.vars = vars(unname(c(tfCols$hasPeak[c(tf1, tf2)], polIICols$is_expressed[c(polII1, polII2)])))) %>% 
  dplyr::summarise(n = n())

tfPolIIDf <- dplyr::filter_at(.tbl = expressionData,
                              .vars = unname(c(tfCols$hasPeak[c(tf1, tf2)], polIICols$is_expressed[c(polII1, polII2)])),
                              .vars_predicate = any_vars(. == TRUE))

tfPolIIDf$group <- tfPolIIDf %>% 
  dplyr::group_by_at(.vars = vars(unname(c(tfCols$hasPeak[c(tf1, tf2)], polIICols$is_expressed[c(polII1, polII2)])))) %>% 
  dplyr::group_indices()

tfpolII_clusters <- dplyr::select(tfPolIIDf, gene, group) %>% 
  dplyr::rename(cluster = group)

tfPolIIDf %>% 
  dplyr::group_by_at(
    .vars = vars(unname(c(tfCols$hasPeak[c(tf1, tf2)], polIICols$is_expressed[c(polII1, polII2)])), "group")) %>% 
  dplyr::summarise(n = n())

## profile heatmap
multiProf_tfPolII <- multi_profile_plots(exptInfo = tfData[tfData$sampleId %in% tfIds, ],
                                         genesToPlot = tfPolIIDf$gene,
                                         matSource = matrixType,
                                         matBins = matrixDim,
                                         profileColors = tfColorList,
                                         column_title_gp = gpar(fontsize = 12),
                                         clusters = tfpolII_clusters)


## Scalled TF diff profile
scalledTfDiffProf_tfPolII <- profile_heatmap(profileMat = tfLfcMat[tfPolIIDf$gene, ],
                                             signalName = comparisonName,
                                             profileColor = scalledTfDiffColor,
                                             geneGroups = tfpolII_clusters)

## polII signal heatmap
polIIMat_tfPolII <- polIIMat[tfPolIIDf$gene, ]

polIIht_tfPolII <- signal_heatmap(log2_matrix = polIIMat_tfPolII,
                                  htName = "polII_exp",
                                  col_title = polII_ids,
                                  legend_title = "log2(polII_singal)",
                                  color = polII_color,
                                  cluster_columns = FALSE)

## polII fold change heatmap
lfc_tfPolII <- lfcMat[tfPolIIDf$gene, ]

lfcHt_tfPolII <- signal_heatmap(log2_matrix = lfc_tfPolII,
                                htName = "polII_lfc",
                                col_title = purrr::map_chr(polIIDiffPairs, "title"),
                                legend_title = "log2(fold change)",
                                color = lfc_color,
                                cluster_columns = FALSE)

## gene length annotation
gl_tfPolII <- gene_length_heatmap_annotation(bedFile = file_genes, genes = tfPolIIDf$gene)

htlist_tfPolII <- gl_tfPolII$an +
  multiProf_tfPolII$heatmapList +
  scalledTfDiffProf_tfPolII$heatmap +
  polIIht_tfPolII + lfcHt_tfPolII



if( all(rownames(multiProf_tfPolII$profileHeatmaps[[1]]$heatmap@matrix) == tfPolIIDf$gene) ){
  ## order by peak distance
  rowOrd <- order(tfPolIIDf[[ tfCols$peakDist[tf1] ]], tfPolIIDf[[ polIIDiffPairs$p1$name ]], decreasing = TRUE)
}


title_tfPolII <- "Differential binding of kdmB at 20h and 48h: genes with macs2 peaks or top 10% polII signal"


# draw Heatmap and add the annotation name decoration
png(filename = paste0(outPrefix_tfPolII, "_profile2.png", collapse = ""), width=5000, height=3500, res = 270)

draw(htlist_tfPolII,
     main_heatmap = tfData$profileName[1],
     # annotation_legend_list = list(profile1$legend),
     column_title = title_tfPolII,
     column_title_gp = gpar(fontsize = 14, fontface = "bold"),
     row_sub_title_side = "left",
     gap = unit(6, "mm"),
     row_order = rowOrd,
     padding = unit(rep(0.5, times = 4), "cm")
)

## decorate the annotations
add_annotation_titles(annotations = c("gene_length"), anTitle = anLables, fontSize = 12)

row_annotation_axis(an = "gene_length",
                    at = c(0, 2000, 4000),
                    labels = c("0kb", "2kb", ">4kb"),
                    slice = length(unique(tfpolII_clusters$cluster)))

dev.off()


##################################################################################

## polII expressed genes
polIIDf <- dplyr::filter_at(.tbl = expressionData,
                            .vars = unname(c(polIICols$is_expressed[c(polII1, polII2)])),
                            .vars_predicate = any_vars(. == TRUE))

plot(frequency(polIIDf[[polIIDiffPairs$p1$name]]))









p <- plot_scatter(df = tfPolIIDf,
                  s1 = "enrichment.An_kdmB_20h_HA_1",
                  s2 = "enrichment.An_kdmB_48h_HA_1",
                  colorCol = "black",
                  transformation = "log2")

ggplot(data = tfPolIIDf, mapping = aes(x = enrichment.An_kdmB_20h_HA_1, y = enrichment.An_kdmB_20h_HA_2)) +
  geom_point() +
  stat_smooth(method = "lm") +
  geom_quantile() +
  scale_x_continuous(trans = "log2") +
  scale_y_continuous(trans = "log2")


cor(x = tfPolIIDf$enrichment.An_kdmB_48h_HA_2,
    y = tfPolIIDf$enrichment.An_kdmB_48h_HA_1,
    method = "spearman", use = "complete.obs")




