library(chipmine)
library(org.Anidulans.eg.db)
library(here)

## generate profile plots for kdmB complex members
## group the genes into different categories based on kdmB complex member binding status
## plot the profile plots for groups
## plot average signal for the groups

rm(list = ls())

source(file = "E:/Chris_UM/Codes/GO_enrichment/topGO_functions.R")


## IMP: the first sampleID will be treated primary and clustering will be done/used for/of this sample
comparisonName <- "kdmB_complex_20h"
outPrefix <- here::here("kdmB_analysis", comparisonName, comparisonName)

file_plotSamples <- here::here("kdmB_analysis", comparisonName, "samples.txt")

# "deeptools", "miao", "normalizedmatrix", "normalizedmatrix_5kb"
matrixType <- "normalizedmatrix_5kb"
matrixDim = c(500, 200, 100, 10)
showExpressionHeatmap = TRUE

## genes to read
file_exptInfo <- here::here("data", "referenceData/sampleInfo.txt")
file_genes <- here::here("data", "referenceData/AN_genesForPolII.bed")
file_topGoMap <- "E:/Chris_UM/Database/A_Nidulans/ANidulans_OrgDb/geneid2go.ANidulans.topGO.map"
file_geneInfo <- "E:/Chris_UM/Database/A_Nidulans/A_nidulans_FGSC_A4_geneClasses.txt"

TF_dataPath <- here::here("data", "TF_data")
polII_dataPath <- here::here("data", "polII_data")
hist_dataPath <- here::here("data", "histone_data")

orgDb <- org.Anidulans.eg.db

anLables <- list()
outPrefix_all <- paste0(outPrefix, "_allGenes", collapse = "")
outPrefix_expressed <- paste0(outPrefix, "_expressedGenes", collapse = "")
outPrefix_sm <- paste0(outPrefix, "_SM_genes", collapse = "")
outPrefix_peaks <- paste0(outPrefix, "_peaksGenes", collapse = "")
outPrefix_pkExp <- paste0(outPrefix, "_pkExpGenes", collapse = "")

##################################################################################

sampleList <- suppressMessages(readr::read_tsv(file = file_plotSamples, col_names = T, comment = "#"))

## genes to read
geneSet <- suppressMessages(
  readr::read_tsv(file = file_genes, col_names = c("chr", "start", "end", "gene", "score", "strand"))
) %>% 
  dplyr::mutate(length = end - start)

geneDesc <- AnnotationDbi::select(x = orgDb, keys = geneSet$gene, columns = "DESCRIPTION", keytype = "GID")

geneSet <- dplyr::left_join(x = geneSet, y = geneDesc, by = c("gene" = "GID"))

## gene information annotations: cluster and TF and polII expression values
geneInfo <- add_gene_info(file = file_geneInfo, clusterDf = geneSet)

head(geneInfo)

##################################################################################

## read the experiment sample details and select only those which are to be plotted
exptData <- get_sample_information(exptInfoFile = file_exptInfo,
                                   samples = sampleList$sampleId,
                                   dataPath = TF_dataPath,
                                   matrixSource = matrixType)


polII_ids <- exptData$sampleId[which(exptData$IP_tag == "polII")]
tfIds <- exptData$sampleId[which(exptData$IP_tag %in% c("HA", "MYC", "TAP") & exptData$TF != "untagged")]


polIICols <- list(
  exp = structure(polII_ids, names = polII_ids),
  is_expressed = structure(paste("is_expressed", ".", polII_ids, sep = ""), names = polII_ids)
)


tfCols <- sapply(c("hasPeak", "pval", "peakType", "tesPeakType", "peakDist", "summitDist", "upstreamExpr", "peakExpr", "relativeDist"),
                 FUN = function(x){ structure(paste(x, ".", tfIds, sep = ""), names = tfIds) },
                 simplify = F, USE.NAMES = T)

expressionData <- get_TF_binding_data(exptInfo = exptData,
                                      genesDf = geneInfo)

dplyr::group_by_at(expressionData, .vars = vars(starts_with("hasPeak."))) %>% 
  dplyr::summarise(n = n())


hasPeakDf <- expressionData %>% 
  dplyr::filter_at(.vars = vars(starts_with("hasPeak")), .vars_predicate = any_vars(. == TRUE)) %>% 
  dplyr::mutate(group = group_indices(., !!! lapply(unname(tfCols$hasPeak), as.name))) %>% 
  dplyr::mutate(group = sprintf(fmt = "%02d", group)) %>% 
  as.data.frame()

readr::write_tsv(x = hasPeakDf, path = paste(outPrefix, "_peaks_data.tab", sep = ""))

## group label data mean profile facet plot
groupLabelDf <- dplyr::group_by_at(.tbl = hasPeakDf, .vars = vars(starts_with("hasPeak."), group)) %>%
  dplyr::summarise(n = n()) %>% 
  dplyr::mutate(groupLabels = paste(group, ": ", n, " genes", sep = ""))

groupLabels <- structure(groupLabelDf$groupLabels, names = groupLabelDf$group)

##################################################################################
## topGO enrichment
goEnrich <- dplyr::group_by_at(.tbl = hasPeakDf, .vars = vars(starts_with("hasPeak."), group)) %>%
  do(topGO_enrichment(goMapFile = file_topGoMap, genesOfInterest = .$gene, goNodeSize = 5))


fwrite(x = goEnrich,
       file = paste(outPrefix, "_peakGroups_GO_enrichment.tab", sep = ""), sep = "\t", col.names = T, quote = F)



## clusterProfiler groupGO assignment
grpGo <- dplyr::group_by_at(.tbl = hasPeakDf, .vars = vars(starts_with("hasPeak."), group)) %>%
  do(clusterProfiler_groupGO(genes = .$gene, org = orgDb, goLevel = 3, type = "BP", keyType ="GID"))

fwrite(x = grpGo,
       file = paste(outPrefix, "_peakGroups_GO_assignment.tab", sep = ""), sep = "\t", col.names = T, quote = F)



## pathway enrichment
keggEnr <- dplyr::group_by_at(.tbl = hasPeakDf, .vars = vars(starts_with("hasPeak."), group)) %>%
  do(keggprofile_enrichment(genes = .$gene, orgdb = orgDb, keytype = "GID", keggOrg = "ani", pvalCut = 0.05))


fwrite(x = keggEnr,
       file = paste(outPrefix, "_peakGroups_KEGG_enrichment.tab", sep = ""), sep = "\t", col.names = T, quote = F)


##################################################################################

## profile matrix of the genes which show binding
matList <- profile_matrix_list(exptInfo = exptData, geneList = geneInfo$gene, source = matrixType,
                               up = matrixDim[1], target = matrixDim[2], down = matrixDim[3])

## get average signal over all factors to select color
meanSignal <- getSignalsFromList(lt = matList[tfIds])
quantile(meanSignal, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)
meanCol <- colorRamp2(quantile(meanSignal, c(0.50, 0.99), na.rm = T), c("white", "red"))

colorList <- sapply(X = exptData$sampleId, FUN = function(x){return(meanCol)})

##################################################################################
ylimList <- sapply(exptData$sampleId, function(x){return(c(0,25))}, simplify = FALSE)

profiles_peaks <- multi_profile_plots(exptInfo = exptData,
                                      genesToPlot = hasPeakDf$gene,
                                      clusters = NULL,
                                      profileColors = colorList,
                                      matSource = matrixType,
                                      matBins = matrixDim,
                                      ylimFraction = ylimList,
                                      column_title_gp = gpar(fontsize = 12)
)

anGl_peaks <- gene_length_heatmap_annotation(bedFile = file_genes, genes = hasPeakDf$gene)

peaks_htlist <- anGl_peaks$an + profiles_peaks$heatmapList

pdfWd <- 2 + (length(profiles_peaks$heatmapList@ht_list) * 2) +
  (length(polII_ids) * 0.25 * showExpressionHeatmap) + 2

title_peak= "Transcription factor binding profile: kdmB complex TF bound genes"

# draw Heatmap and add the annotation name decoration
pdf(file = paste(outPrefix, "_peaks_unclustered.pdf", sep = ""), width = pdfWd, height = 12)

draw(peaks_htlist,
     main_heatmap = exptData$profileName[1],
     # annotation_legend_list = list(profile1$legend),
     column_title = title_peak,
     column_title_gp = gpar(fontsize = 14, fontface = "bold"),
     row_sub_title_side = "left",
     heatmap_legend_side = "right",
     gap = unit(5, "mm"),
     padding = unit(rep(0.5, times = 4), "cm")
)


row_annotation_axis(an = "gene_length",
                    at = c(0, 2000, 4000),
                    labels = c("0kb", "2kb", ">4kb"),
                    slice = 1)

dev.off()



##################################################################################

## profile matrix with selected genes only
newClusters <- dplyr::select(hasPeakDf, gene, group) %>% 
  dplyr::rename(cluster = group)

profiles_peakGroups <- multi_profile_plots(exptInfo = exptData,
                                          genesToPlot = hasPeakDf$gene,
                                          clusters = newClusters,
                                          clustOrd = sort(unique(newClusters$cluster)),
                                          profileColors = colorList,
                                          matSource = matrixType,
                                          matBins = matrixDim,
                                          column_title_gp = gpar(fontsize = 12)
)



## heatmap of binary assignment of samples to different group
htMat <- dplyr::select(hasPeakDf, gene, starts_with("hasPeak")) %>% 
  dplyr::mutate_if(.predicate = is.logical, .funs = as.character) %>% 
  tibble::column_to_rownames("gene")

## column name as annotation for Heatmap
colNameAnn <- HeatmapAnnotation(
  colName = anno_text(
    x = unname(tfCols$hasPeak),
    rot = 90, just = "left",
    offset = unit(1, "mm"),
    gp = gpar(fontsize = 10)),
  annotation_height = unit.c(max_text_width(unname(tfCols$hasPeak)))
)

grp_ht <- Heatmap(htMat,
                  col = c("TRUE" = "black", "FALSE" = "white"),
                  heatmap_legend_param = list(title = "Peak detected"),
                  # column_names_side = "top",
                  show_column_names = FALSE,
                  top_annotation = colNameAnn,
                  cluster_columns = FALSE, cluster_rows = FALSE,
                  width = unit(3, "cm"),
                  show_row_names = FALSE
)

## gene length annotation
anGl_peakGroups <- gene_length_heatmap_annotation(bedFile = file_genes, genes = hasPeakDf$gene)


peakGroups_htlist <- anGl_peakGroups$an + profiles_peakGroups$heatmapList + grp_ht

pdfWd <- 2 + 
  (length(profiles_peakGroups$heatmapList@ht_list) * 2) +
  (length(polII_ids) * 0.25 * showExpressionHeatmap) + 1

title_peak= "Transcription factor binding profile: kdmB complex TF bound genes"

# draw Heatmap and add the annotation name decoration
pdf(file = paste(outPrefix, "_peak_comparison.pdf", sep = ""), width = pdfWd, height = 12)

draw(peakGroups_htlist,
     main_heatmap = exptData$profileName[1],
     # annotation_legend_list = list(profile1$legend),
     column_title = title_peak,
     column_title_gp = gpar(fontsize = 14, fontface = "bold"),
     row_sub_title_side = "left",
     heatmap_legend_side = "right",
     gap = unit(5, "mm"),
     padding = unit(rep(0.5, times = 4), "cm")
)


row_annotation_axis(an = "gene_length",
                    at = c(0, 2000, 4000),
                    labels = c("0kb", "2kb", ">4kb"),
                    slice = length(unique(newClusters$cluster)))

dev.off()


##################################################################################
## average line plots for peak groups
testDf <- dplyr::filter(hasPeakDf, group == "07")

lineColors = structure(.Data = c(RColorBrewer::brewer.pal(n = 4, name = "Set1"), "black"),
                       names = exptData$sampleId)

lineShape = structure(.Data = c(1, 1, 1, 1, 6),
                      names = exptData$sampleId)


draw_avg_profile_plot(exptInfo = exptData,
                      profileMats = matList,
                      genes = testDf$gene,
                      lineColors = lineColors,
                      lineShape = lineShape)


ap <- geneset_average_profile(exptInfo = exptData,
                              profileMats = matList,
                              genes = testDf$gene,
                              cluster = "group_1")


## average profile data for each group
groupMeanProfiles <- dplyr::group_by_at(.tbl = hasPeakDf, .vars = vars(starts_with("hasPeak."), group)) %>%
  dplyr::do(
    geneset_average_profile(exptInfo = exptData,
                            profileMats = matList,
                            genes = .$gene,
                            cluster = unique(.$group))
  ) %>% 
  dplyr::ungroup() %>% 
  as.data.frame()

## decide the axis labels
axisBrk <- NULL
axisLab <- NULL
targetEnd <- NULL
if(attributes(matList[[1]])$target_is_single_point){
  axisBrk <- c(
    attributes(matList[[1]])$upstream_index[1],
    attributes(matList[[1]])$target_index[1],
    tail(attributes(matList[[1]])$downstream_index, 1)
  )
  
  axisLab <- c(
    -attributes(matList[[1]])$extend[1],
    attributes(matList[[1]])$target_name,
    attributes(matList[[1]])$extend[2]
  )
  
  targetEnd <- axisBrk[2] + 1
  
} else if(! attributes(matList[[1]])$target_is_single_point){
  axisBrk <- c(
    attributes(matList[[1]])$upstream_index[1],
    attributes(matList[[1]])$target_index[1],
    tail(attributes(matList[[1]])$target_index, 1),
    tail(attributes(matList[[1]])$downstream_index, 1)
  )
  
  axisLab <- c(
    -attributes(matList[[1]])$extend[1],
    "START", "END",
    attributes(matList[[1]])$extend[2]
  )
  
  targetEnd <- axisBrk[3]
}

## plot the groups using facets
p = ggplot(data = groupMeanProfiles) +
  geom_line(mapping = aes(x = bin, y = mean, group = sample, color = sample, linetype = sample),
            size = 0.8, alpha = 0.8) +
  geom_hline(yintercept = 0, size = 2, color = "grey70") +
  geom_segment(mapping = aes(x = axisBrk[2], y = 0, xend = targetEnd, yend = 0),
               size = 10, lineend = "butt", color = "grey50") +
  scale_color_manual(values = lineColors) +
  scale_linetype_manual(values = lineShape) +
  scale_x_continuous(breaks = axisBrk, labels = axisLab) +
  ggtitle(paste("Average profile for different binding patterns of kdmB complex members: ", comparisonName)) +
  ylab("Read coverage") +
  facet_wrap(group ~ ., ncol = 5, scales = "free_y", labeller = labeller(group = groupLabels)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = element_text(size = 13, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 15),
        strip.text = element_text(size = 12, face = "bold"),
        legend.position = "bottom",
        legend.spacing.x = unit(0.5, 'cm'),
        # legend.justification = c(1.1, 1.1),
        legend.title = element_blank(),
        legend.text = element_text(size = 13),
        plot.margin = unit(rep(0.5, 4), "cm")
  ) +
  guides(
    color = guide_legend(
      label.position = "bottom")
  )


pdf(file = paste(outPrefix, "_mean_profiles.pdf", sep = ""), width = 18, height = 10)
p
dev.off()






