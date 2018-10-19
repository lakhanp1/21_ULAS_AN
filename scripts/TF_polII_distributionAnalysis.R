library(dplyr)
library(data.table)
library(tibble)
library(ggplot2)
library(ggpubr)
library(ltm)
library(VennDiagram)
library(RColorBrewer)
library(chipmine)


rm(list = ls())

source(file = "E:/Chris_UM/Codes/GO_enrichment/topGO_functions.R")


path = "E:/Chris_UM/Analysis/21_CL2017_ChIPmix_ULAS_MIX/ULAS_AN"
setwd(path)

## This script is for exploratory analysis of the TF and polII expression data
mapFile = "E:/Chris_UM/Database/A_Nidulans/geneid2go.ANidulans.20171004.map"

##################################################################################
tfSample = "An_kdmB_48h_HA"
polIISample = "An_untagged_48h_polII"
name = "An_kdmB_48h_HA"


TF_dataPath = paste("TF_data/", tfSample, sep = "")

clusterFile = paste(TF_dataPath, "/", name, "_allGenes_clusters.tab", sep = "")
outDir = paste(TF_dataPath, "/distributionAnalysis", sep = "")

if(!dir.exists(outDir)){
  dir.create(outDir)
}

outPrefix_polIIDist = paste(outDir, "/", name, "_polIIDistribution", sep = "")

##################################################################################

clusterData = fread(file = clusterFile, sep = "\t",
                    header = T, stringsAsFactors = F, na.strings = "NA", data.table = F)


polIIIsExpCol = paste("is_expressed.", polIISample, sep = "")
hasPeakCol = paste("hasPeak.", tfSample, sep = "")
polIIExpCol = polIISample
tfExpCol = paste("peakExpr.", tfSample, sep = "")

##################################################################################
## plot polII expression distribution violin plot with respect to TF binding status
distViolin = polIIExp_vs_TFbinding_violin(df = clusterData,
                                          tf = tfSample,
                                          polII = polIISample,
                                          isExpCol = polIIIsExpCol,
                                          isBoundCol = hasPeakCol)

# draw Heatmap and add the annotation name decoration
png(filename = paste(outPrefix_polIIDist, ".png", sep = ""), width=4000, height=3500, res = 350)

print(distViolin)

dev.off()


##################################################################################
## venn diagram to show polII expression and TF binding overlap
geneList = list(
  polIIExp = unique(clusterData$gene[ which(clusterData[[polIIIsExpCol]]) ]),
  tfBound = unique(clusterData$gene[ which(clusterData[[hasPeakCol]]) ])
)


## plot Venn diagram for variants
vennOut = paste(outPrefix_polIIDist, "_venn.png", sep = "")

vennAll = venn.diagram(x = geneList,
                       category.names = c(polIIIsExpCol, hasPeakCol),
                       main = "polII expressed genes and TF bound genes",
                       imagetype = "png",
                       filename = vennOut,
                       height = 2000,
                       width = 4000,
                       resolution = 300,
                       euler.d = T,
                       fill = c("blue", "green"),
                       alpha = 0.5,
                       fontface = "bold",
                       main.cex = 1.5,
                       cex = 1.5,
                       cat.fontface = "bold",
                       cat.cex = 1,
                       margin = 0.2
)

## remove the vennDiagram log file
unlink(paste(vennOut, ".*.log", sep = ""))


##################################################################################
## GO enrichment for each set

polIIExpAndTfBound = intersect(x = geneList$polIIExp, geneList$tfBound)

## get GO enrichment table
goData = get_topGO_enrichment(goMapFile = mapFile, genesOfInterest = polIIExpAndTfBound)

fwrite(x = goData, file = paste(outPrefix_polIIDist, "_GO.tab", sep = ""),
       sep = "\t", col.names = T, quote = F)

goTitle = paste("GO enrichment for Genes expressed in ", polIISample, "\nAND bound by", tfSample, sep = " ")

topGoScatter = topGO_scatterPlot(df = goData, title = goTitle)

# draw Heatmap and add the annotation name decoration
ht = max(nrow(goData) * 80, 1500)
wd = (min(max(nchar(as.character(goData$Term))), 80) * 30) * 1.5
rs = max(min(wd, ht) / 12, 200)

png(filename = paste(outPrefix_polIIDist, "_GO.png", sep = ""), width = wd, height = ht, res = rs)

print(topGoScatter)

dev.off()


##################################################################################

##scatter plot for polII expression and TF binding level at promoter/peak region

scatterPlot = polII_tf_expression_scatter(df = clusterData,
                                          tf = tfSample,
                                          polII = polIISample,
                                          tfExp = tfExpCol,
                                          polIIExp = polIIExpCol,
                                          isExpCol = polIIIsExpCol,
                                          isBoundCol = hasPeakCol)


png(filename = paste(outPrefix_polIIDist, "_polIITf_scatter.png", sep = ""), width = 5000, height = 5000, res = 450)

print(scatterPlot)

dev.off()


##################################################################################








##################################################################################

exp = factor(clusterData$has_TSS_peak, levels = sort(unique(clusterData$has_TSS_peak), decreasing = T))

biserial.cor(x = clusterData$An_untagged_48h_polII, y = clusterData$has_TSS_peak)

biserial.cor(x = clusterData$An_untagged_48h_polII, y = exp)

polyserial(x = clusterData$An_untagged_48h_polII, y = clusterData$has_TSS_peak, std.err = T)

polyserial(x = clusterData$An_untagged_48h_polII, y = clusterData$has_TSS_peak, ML = T, std.err = T)


polyserial(x = clusterData$An_untagged_48h_polII, y = clusterData$`is_expressed(An_untagged_48h_polII)`,
           std.err = T)


##################################################################################









