library(chipmine)
library(org.Anidulans.eg.db)
library(here)

rm(list = ls())

source(file = "E:/Chris_UM/Codes/GO_enrichment/topGO_functions.R")


## IMP: the first sampleID will be treated primary and clustering will be done/used for/of this sample
comparisonName <- "kdmB_complex_48h"
sampleList <- c("An_kdmB_48h_HA_1",
                "An_sntB_48h_HA_1",
                "An_ecoA_48h_HA_1",
                "An_rpdA_48h_HA_1",
                "An_untagged_48h_HA_1")

# path <- "E:/Chris_UM/Analysis/21_CL2017_ChIPmix_ULAS_MIX/ULAS_AN/kdmB_analysis/kdmB_complex_48h"
# setwd(path)


outPrefix <- here::here("kdmB_analysis", "kdmB_complex_48h", comparisonName)

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
anLables[["gene_length"]] = "Gene Length"
##################################################################################

## genes to read
geneSet <- data.table::fread(file = file_genes, header = F,
                             col.names = c("chr", "start", "end", "gene", "score", "strand")) %>% 
  dplyr::mutate(length = end - start)

geneDesc <- select(x = orgDb, keys = geneSet$gene, columns = "DESCRIPTION", keytype = "GID")

geneSet <- dplyr::left_join(x = geneSet, y = geneDesc, by = c("gene" = "GID"))

## gene information annotations: cluster and TF and polII expression values
geneInfo <- add_gene_info(file = file_geneInfo, clusterDf = geneSet)

head(geneInfo)

##################################################################################

## read the experiment sample details and select only those which are to be plotted
exptData <- get_sample_information(exptInfoFile = file_exptInfo,
                                   samples = sampleList,
                                   dataPath = TF_dataPath,
                                   matrixSource = "deeptools")


polII_ids <- exptData$sampleId[which(exptData$IP_tag == "polII")]
tfIds <- exptData$sampleId[which(exptData$IP_tag %in% c("HA", "MYC", "TAP") & exptData$TF != "untagged")]

polII_expIds <- paste("is_expressed.", polII_ids, sep = "")
hasPeakCols <- paste("hasPeak.", tfIds, sep = "")


expressionData <- get_TF_binding_data(exptInfo = exptData,
                                      genesDf = geneInfo)

dplyr::group_by_at(expressionData, .vars = vars(starts_with("hasPeak."))) %>% 
  dplyr::summarise(n = n())


dt <- expressionData %>% 
  dplyr::filter_at(.vars = vars(starts_with("hasPeak")), .vars_predicate = any_vars(. == TRUE)) %>% 
  dplyr::mutate(group = group_indices(., !!! lapply(hasPeakCols, as.name)))


## read the cluster information for the first sample
# clusterData <- data.table::fread(file = exptData$clusterFile[1], sep = "\t", header = T, stringsAsFactors = F)
# 
# dt <- dplyr::left_join(dt, clusterData, by = c("gene" = "gene"))

rownames(dt) <- dt$gene

## group label data mean profile facet plot
groupLabelDf <- dplyr::group_by_at(.tbl = dt, .vars = vars(starts_with("hasPeak."), group)) %>%
  dplyr::summarise(n = n()) %>% 
  dplyr::mutate(groupLabels = paste(group, ": ", n, " genes", sep = ""))

groupLabels <- structure(groupLabelDf$groupLabels, names = groupLabelDf$group)

##################################################################################
## topGO enrichment
goEnrich <- dplyr::group_by_at(.tbl = dt, .vars = vars(starts_with("hasPeak."), group)) %>%
  do(topGO_enrichment(goMapFile = file_topGoMap, genesOfInterest = .$gene, goNodeSize = 5))


fwrite(x = goEnrich,
       file = paste(outPrefix, "_GO_enrichment.tab", sep = ""), sep = "\t", col.names = T, quote = F)



## clusterProfiler groupGO assignment
grpGo <- dplyr::group_by_at(.tbl = dt, .vars = vars(starts_with("hasPeak."), group)) %>%
  do(clusterProfiler_groupGO(genes = .$gene, org = orgDb, goLevel = 3, type = "BP", keyType ="GID"))

fwrite(x = grpGo,
       file = paste(outPrefix, "_GO_assignment.tab", sep = ""), sep = "\t", col.names = T, quote = F)



## pathway enrichment
keggEnr <- dplyr::group_by_at(.tbl = dt, .vars = vars(starts_with("hasPeak."), group)) %>%
  do(keggprofile_enrichment(genes = .$gene, orgdb = orgDb, keytype = "GID", keggOrg = "ani", pvalCut = 0.05))


fwrite(x = keggEnr,
       file = paste(outPrefix, "_KEGG_enrichment.tab", sep = ""), sep = "\t", col.names = T, quote = F)


##################################################################################

## heatmap of binary assignment of samples to different group
htMat <- dplyr::select(dt, gene, starts_with("hasPeak")) %>% 
  dplyr::mutate_if(.predicate = is.logical, .funs = as.character) %>% 
  tibble::column_to_rownames("gene")

## column name as annotation for Heatmap
colNameAnn <- HeatmapAnnotation(colName = anno_text(x = hasPeakCols,
                                                    rot = 90, just = "left",
                                                    offset = unit(1, "mm"),
                                                    gp = gpar(fontsize = 10))
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

##################################################################################
## profile matrix
## generate profile plot for all the genes
multiProfiles_all <- multi_profile_plots(exptInfo = exptData,
                                         genesToPlot = geneInfo$gene,
                                         clusters = NULL)

matList <- profile_matrix_list(exptInfo = exptData, geneList = geneInfo$gene, source = "deeptools")

## get average signal over all factors to select color
meanSignal <- getSignalsFromList(lt = matList[tfIds])
quantile(meanSignal, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)
meanCol <- colorRamp2(quantile(meanSignal, c(0.50, 0.995), na.rm = T), c("white", "red"))

colorList <- sapply(X = exptData$sampleId, FUN = function(x){return(meanCol)})

## profile matrix with selected genes only
newClusters <- dplyr::select(dt, gene, group) %>% 
  dplyr::rename(cluster = group)

multiProfiles_peak <- multi_profile_plots(exptInfo = exptData,
                                          genesToPlot = dt$gene,
                                          clusters = newClusters,
                                          clustOrd = sort(unique(newClusters$cluster)),
                                          profileColors = colorList
                                          )


## gene length annotation
anGl_peaks <- gene_length_heatmap_annotation(bedFile = file_genes, genes = dt$gene)


peaks_htlist <- anGl_peaks$an + multiProfiles_peak$heatmapList + grp_ht

wd <- 500 + (nrow(exptData) * 700) + (length(polII_ids)) + 200
title_peak= "Transcription factor binding profile: kdmB complex TF bound genes"

# draw Heatmap and add the annotation name decoration
png(filename = paste(outPrefix, "_peak_comparison.png", sep = ""), width=wd, height=3500, res = 270)

draw(peaks_htlist,
     main_heatmap = exptData$profileName[1],
     # annotation_legend_list = list(profile1$legend),
     column_title = title_peak,
     column_title_gp = gpar(fontsize = 14, fontface = "bold"),
     row_sub_title_side = "left",
     heatmap_legend_side = "bottom",
     gap = unit(c(5, rep(10, length(exptData$sampleId)), 5), "mm"),
     padding = unit(rep(0.5, times = 4), "cm")
)

add_annotation_titles(annotations = c("gene_length"), anTitle = anLables)

row_annotation_axis(an = "gene_length",
                    at = c(0, 2000, 4000),
                    labels = c("0kb", "2kb", ">4kb"),
                    slice = length(unique(newClusters$cluster)))

dev.off()


##################################################################################
## line plots
testDf <- dplyr::filter(dt, group == 7)

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
groupMeanProfiles <- dplyr::group_by_at(.tbl = dt, .vars = vars(starts_with("hasPeak."), group)) %>%
  dplyr::do(
    geneset_average_profile(exptInfo = exptData,
                            profileMats = matList,
                            genes = .$gene,
                            cluster = unique(.$group))
  ) %>% 
  dplyr::ungroup() %>% 
  as.data.frame()



## plot the groups using facets
p = ggplot(data = groupMeanProfiles) +
  geom_line(mapping = aes(x = bin, y = mean, group = sample, color = sample, linetype = sample),
            size = 1) +
  geom_segment(mapping = aes(x = 200, y = 0, xend = 400, yend = 0), size = 10, lineend = "butt", color = "grey70") +
  geom_hline(yintercept = 0, size = 2, color = "grey70") +
  annotate(geom = "text", x = 230, y = 0, label = "gene", hjust = 0, vjust = 0.3) +
  scale_color_manual(values = lineColors) +
  scale_linetype_manual(values = lineShape) +
  scale_x_continuous(breaks = c(0, 100, 200, 400, 500),
                     labels = c("-2kb", "-1kb", "ATG", "STOP", "+1kb")) +
  ggtitle(paste("Average profile for different binding patterns of kdmB complex members: ", outPrefix)) +
  ylab("Read coverage") +
  facet_wrap(group ~ ., ncol = 5, scales = "free_y", labeller = labeller(group = groupLabels)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 15),
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


png(filename = paste(outPrefix, "_mean_profiles.png", sep = ""), width = 5500, height = 3000, res = 320)
p
dev.off()






