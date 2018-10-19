library(chipmine)
library(org.Anidulans.eg.db)
library(foreach)
library(doParallel)

rm(list = ls())


path <- "E:/Chris_UM/Analysis/21_CL2017_ChIPmix_ULAS_MIX/ULAS_AN/data"
setwd(path)

cl <- makeCluster(4) #not to overload your computer
registerDoParallel(cl)

##################################################################################

exptInfoFile <- "referenceData/sampleInfo.txt"
TF_dataPath <- "TF_data"
polII_dataPath <- "polII_data"
histone_dataPath <- "histone_data"
otherDataPath <- "other_data"
geneCdsFile <- "E:/Chris_UM/Database/A_Nidulans/A_nidulans_FGSC_A4_version_s10-m04-r03_CDS_Unique.bed"
cdsUpRegionFile <- "E:/Chris_UM/Database/A_Nidulans/A_nidulans_FGSC_A4_version_s10-m04-r03_CDS_upstream500.bed"
file_genes <- "referenceData/AN_genesForPolII.bed"
orgDb <- org.Anidulans.eg.db


geneSet <- data.table::fread(file = file_genes, header = F,
                             col.names = c("chr", "start", "end", "gene", "score", "strand")) %>% 
  dplyr::select(-score) %>% 
  dplyr::mutate(length = end - start)

geneDesc <- AnnotationDbi::select(x = orgDb, keys = geneSet$gene, columns = "DESCRIPTION", keytype = "GID")

geneSet <- dplyr::left_join(x = geneSet, y = geneDesc, by = c("gene" = "GID"))

##################################################################################

# polIIsampleFile <- paste(polII_dataPath, "/", "polII_sample.list", sep = "")
# 
# polIISampleList <- fread(file = polIIsampleFile, sep = "\t", header = F,
#                          stringsAsFactors = F, col.names = c("id"), data.table = F)
# 
# # polIISampleList = data.frame(id = c("An_untagged_20h_polII_1", "An_untagged_48h_polII_1"), stringsAsFactors = F)
# 
# polII_info <- get_sample_information(exptInfoFile = exptInfoFile,
#                                      samples = polIISampleList$id,
#                                      dataPath = polII_dataPath,
#                                      matrixSource = "normalizedmatrix")
# 
# i <- 1
# 
# ## process all polII expression matrix
# foreach(i = 1:nrow(polII_info)) %dopar% {
#   
#   # polIIDf <- preProcess_polII_expression(expMat = polII_info$polIIExpMat[i],
#   #                                        title = polII_info$sampleId[i],
#   #                                        expFraction = 10,
#   #                                        polIIExpFile = polII_info$polIIExpFile[i])
#   
#   ## 2kb - 2kb - 1kb matrix
#   bwMat <- chipmine::bigwig_profile_matrix(bwFile = polII_info$bwFile[i],
#                                            bedFile = file_genes,
#                                            signalName = polII_info$sampleId[i],
#                                            genes = geneSet$gene,
#                                            readLocal = FALSE,
#                                            storeLocal = TRUE,
#                                            localPath = polII_info$matFile[i])
#   
#   mat5Kb <- gsub(pattern = ".tab.gz", replacement = "_5kb.tab.gz", x = polII_info$matFile[i])
#   
#   ## 5kb - 2kb - 1kb matrix
#   bwMat <- chipmine::bigwig_profile_matrix(bwFile = polII_info$bwFile[i],
#                                            bedFile = file_genes,
#                                            signalName = polII_info$sampleId[i],
#                                            genes = geneSet$gene,
#                                            readLocal = FALSE,
#                                            storeLocal = TRUE,
#                                            localPath = mat5Kb,
#                                            extend = c(5000, 1000),
#                                            target_ratio = 0.25)
#   
#   polII_info$bwFile[i]
# }


##################################################################################
## process all TF macs2 results to merge everything

tfSampleFile <- paste(TF_dataPath, "/", "tf_macs2_samples.list", sep = "")
# tfSampleFile <- paste(TF_dataPath, "/", "tf_samples.list", sep = "")

tfSampleList <- readr::read_tsv(file = tfSampleFile, col_names = c("id"),  comment = "#") %>% 
  as.data.frame()

# tfSampleList <- data.frame(id = c("An_untagged_48h_input_1", "An_untagged_20h_HA_1", "An_untagged_20h_HA_2", "An_untagged_48h_HA_1", "An_untagged_48h_HA_2", "An_untagged_20h_MYC_1", "An_untagged_20h_MYC_2", "An_untagged_48h_MYC_1", "An_untagged_48h_MYC_2"),
#                           stringsAsFactors = F)

# tfSampleList <- data.frame(id = c("An_kdmB_20h_HA_1", "An_kdmB_laeA_del_20h_HA_2", "An_sudA_20h_HA_1", "An_rpdA_20h_HA_1", "An_sntB_20h_HA_1", "An_ecoA_20h_HA_1", "An_kdmB_48h_HA_1", "An_rpdA_48h_HA_1", "An_sntB_48h_HA_1", "An_sudA_48h_HA_1", "An_ecoA_48h_HA_1"),
#                            stringsAsFactors = F)


tfInfo <- get_sample_information(exptInfoFile = exptInfoFile,
                                 samples = tfSampleList$id,
                                 dataPath = TF_dataPath,
                                 matrixSource = "normalizedmatrix")

i <- 1

foreach(i = 1:nrow(tfInfo)) %do% {
  
  tfDf <- chipmine::preProcess_macs2_results(title = tfInfo$sampleId[i],
                                             peakAnnoFile = tfInfo$narrowpeakAnno[i],
                                             cdsFile = geneCdsFile,
                                             peakFile = tfInfo$narrowpeakFile[i],
                                             bwFile = tfInfo$bwFile[i],
                                             outFile = tfInfo$tfPeakFile[i],
                                             bindingInGene = FALSE)
  
  # ## 2kb - 2kb - 1kb matrix
  # bwMat <- chipmine::bigwig_profile_matrix(bwFile = tfInfo$bwFile[i],
  #                                          bedFile = file_genes,
  #                                          signalName = tfInfo$sampleId[i],
  #                                          genes = geneSet$gene,
  #                                          readLocal = FALSE,
  #                                          storeLocal = TRUE,
  #                                          localPath = tfInfo$matFile[i])
  # 
  # mat5Kb <- gsub(pattern = ".tab.gz", replacement = "_5kb.tab.gz", x = tfInfo$matFile[i])
  # 
  # ## 5kb - 2kb - 1kb matrix
  # bwMat <- chipmine::bigwig_profile_matrix(bwFile = tfInfo$bwFile[i],
  #                                          bedFile = file_genes,
  #                                          signalName = tfInfo$sampleId[i],
  #                                          genes = geneSet$gene,
  #                                          readLocal = FALSE,
  #                                          storeLocal = TRUE,
  #                                          localPath = mat5Kb,
  #                                          extend = c(5000, 1000),
  #                                          target_ratio = 0.25)
  # 
  # matTss <- gsub(pattern = ".tab.gz", replacement = "_4kbTSS2kb.tab.gz", x = tfInfo$matFile[i])
  # ## -4kb - TSS - 2kb
  # bwMat <- chipmine::bigwig_profile_matrix(bwFile = tfInfo$bwFile[i],
  #                                          bedFile = file_genes,
  #                                          signalName = tfInfo$sampleId[i],
  #                                          genes = geneSet$gene,
  #                                          readLocal = FALSE,
  #                                          storeLocal = TRUE,
  #                                          localPath = matTss,
  #                                          extend = c(4000, 2000),
  #                                          target = "tss")
  # 
  # matTes <- gsub(pattern = ".tab.gz", replacement = "_2kbTES4kb.tab.gz", x = tfInfo$matFile[i])
  # ## -2kb - TES - 4kb
  # bwMat <- chipmine::bigwig_profile_matrix(bwFile = tfInfo$bwFile[i],
  #                                          bedFile = file_genes,
  #                                          signalName = tfInfo$sampleId[i],
  #                                          genes = geneSet$gene,
  #                                          readLocal = FALSE,
  #                                          storeLocal = TRUE,
  #                                          localPath = matTes,
  #                                          extend = c(2000, 4000),
  #                                          target = "tes")
  # 
  # matTss <- gsub(pattern = ".tab.gz", replacement = "_3kbTSS3kb.tab.gz", x = tfInfo$matFile[i])
  # ## -5kb - TSS - 1kb
  # bwMat <- chipmine::bigwig_profile_matrix(bwFile = tfInfo$bwFile[i],
  #                                          bedFile = file_genes,
  #                                          signalName = tfInfo$sampleId[i],
  #                                          genes = geneSet$gene,
  #                                          readLocal = FALSE,
  #                                          storeLocal = TRUE,
  #                                          localPath = matTss,
  #                                          extend = c(3000, 3000),
  #                                          target = "tss")
  # 
  # 
  # matTes <- gsub(pattern = ".tab.gz", replacement = "_3kbTES3kb.tab.gz", x = tfInfo$matFile[i])
  # ## -1kb - TES - 5kb
  # bwMat <- chipmine::bigwig_profile_matrix(bwFile = tfInfo$bwFile[i],
  #                                          bedFile = file_genes,
  #                                          signalName = tfInfo$sampleId[i],
  #                                          genes = geneSet$gene,
  #                                          readLocal = FALSE,
  #                                          storeLocal = TRUE,
  #                                          localPath = matTes,
  #                                          extend = c(3000, 3000),
  #                                          target = "tes")
  # 
  
  tfInfo$bwFile[i]
  
}

# 
# ## build peak region BED file for calculating peak region FPKM
# ## if gene has peak, use the peak region else use 1000bp upstream region
# for (i in 1:nrow(tfInfo)) {
#   genwise_peak_regions_bed(title = tfInfo$sampleId[i],
#                            tfPeakFile = tfInfo$tfPeakFile[i],
#                            cdsUpstreamFile = cdsUpRegionFile,
#                            outFile = tfInfo$tfRegionsBed[i])
# }
# 
# 
# ## generate polII FPKM matrix
# 
# 
# 
# ##################################################################################
# 
# 
# ## specific processing for samples where binding is seen in gene body
# # tfSampleList = c("An_laeA_20h_HA", "An_laeA_48h_HA", "An_kdmB_20h_HA", "An_kdmB_48h_HA")
# 
# samplesWithBindingInGene <- c("An_cclA_20h_HA_1", "An_cclA_20h_HA_2", "An_cclA_48h_HA_1", "An_cclA_48h_HA_2", "An_cclA_kdmA_del_20h_HA_1", "An_cclA_kdmA_del_20h_HA_2", "An_cclA_kdmA_del_48h_HA_1", "An_cclA_kdmA_del_48h_HA_2")
# 
# 
# tfInfo <- get_sample_information(exptInfoFile = exptInfoFile,
#                                  samples = samplesWithBindingInGene,
#                                  dataPath = TF_dataPath,
#                                  matrixSource = "normalizedmatrix")
# 
# 
# for (i in 1:nrow(tfInfo)) {
#   tfDf <- preProcess_macs2_results(title = tfInfo$sampleId[i],
#                                    peakAnnoFile = tfInfo$narrowpeakAnno[i],
#                                    cdsFile = geneCdsFile,
#                                    outFile = tfInfo$tfPeakFile[i],
#                                    bindingInGene = TRUE)
# }
# 
# 
# ## build peak region BED file for calculating peak region FPKM
# for (i in 1:nrow(tfInfo)) {
#   genwise_peak_regions_bed(title = tfInfo$sampleId[i],
#                            tfPeakFile = tfInfo$tfPeakFile[i],
#                            cdsUpstreamFile = cdsUpRegionFile,
#                            outFile = tfInfo$tfRegionsBed[i])
# }
# 

##################################################################################
## histone data



file_histData <- paste(histone_dataPath, "/", "histone_sample.list", sep = "")


histSamples <- readr::read_tsv(file = file_histData, col_names = c("id"),  comment = "#") %>% 
  as.data.frame()

# histSamples <- data.frame(id = c("An_H3_20h_HIST_1", "An_H3_48h_HIST_1", "An_H3_cclA_del_20h_HIST_1",
#                                  "An_H3_cclA_del_48h_HIST_1"),
#                           stringsAsFactors = F)

histInfo <- get_sample_information(exptInfoFile = exptInfoFile,
                                   samples = histSamples$id,
                                   dataPath = histone_dataPath,
                                   matrixSource = "normalizedmatrix")

i <- 1

## process all polII expression matrix
foreach(i = 1:nrow(histInfo)) %dopar% {
  
  
  ## 2kb - 2kb - 1kb matrix
  bwMat <- chipmine::bigwig_profile_matrix(bwFile = histInfo$bwFile[i],
                                           bedFile = file_genes,
                                           signalName = histInfo$sampleId[i],
                                           genes = geneSet$gene,
                                           readLocal = FALSE,
                                           storeLocal = TRUE,
                                           localPath = histInfo$matFile[i])

  mat5Kb <- gsub(pattern = ".tab.gz", replacement = "_5kb.tab.gz", x = histInfo$matFile[i])

  ## 5kb - 2kb - 1kb matrix
  bwMat <- chipmine::bigwig_profile_matrix(bwFile = histInfo$bwFile[i],
                                           bedFile = file_genes,
                                           signalName = histInfo$sampleId[i],
                                           genes = geneSet$gene,
                                           readLocal = FALSE,
                                           storeLocal = TRUE,
                                           localPath = mat5Kb,
                                           extend = c(5000, 1000),
                                           target_ratio = 0.25)
  
  matTss <- gsub(pattern = ".tab.gz", replacement = "_4kbTSS2kb.tab.gz", x = histInfo$matFile[i])
  ## -4kb - TSS - 2kb
  bwMat <- chipmine::bigwig_profile_matrix(bwFile = histInfo$bwFile[i],
                                           bedFile = file_genes,
                                           signalName = histInfo$sampleId[i],
                                           genes = geneSet$gene,
                                           readLocal = FALSE,
                                           storeLocal = TRUE,
                                           localPath = matTss,
                                           extend = c(4000, 2000),
                                           target = "tss")
  
  matTes <- gsub(pattern = ".tab.gz", replacement = "_2kbTES4kb.tab.gz", x = histInfo$matFile[i])
  ## -2kb - TES - 4kb
  bwMat <- chipmine::bigwig_profile_matrix(bwFile = histInfo$bwFile[i],
                                           bedFile = file_genes,
                                           signalName = histInfo$sampleId[i],
                                           genes = geneSet$gene,
                                           readLocal = FALSE,
                                           storeLocal = TRUE,
                                           localPath = matTes,
                                           extend = c(2000, 4000),
                                           target = "tes")
  
  
  matTss <- gsub(pattern = ".tab.gz", replacement = "_3kbTSS3kb.tab.gz", x = histInfo$matFile[i])
  ## -5kb - TSS - 1kb
  bwMat <- chipmine::bigwig_profile_matrix(bwFile = histInfo$bwFile[i],
                                           bedFile = file_genes,
                                           signalName = histInfo$sampleId[i],
                                           genes = geneSet$gene,
                                           readLocal = FALSE,
                                           storeLocal = TRUE,
                                           localPath = matTss,
                                           extend = c(3000, 3000),
                                           target = "tss")
  
  
  
  matTes <- gsub(pattern = ".tab.gz", replacement = "_3kbTES3kb.tab.gz", x = histInfo$matFile[i])
  ## -1kb - TES - 5kb
  bwMat <- chipmine::bigwig_profile_matrix(bwFile = histInfo$bwFile[i],
                                           bedFile = file_genes,
                                           signalName = histInfo$sampleId[i],
                                           genes = geneSet$gene,
                                           readLocal = FALSE,
                                           storeLocal = TRUE,
                                           localPath = matTes,
                                           extend = c(3000, 3000),
                                           target = "tes")
  
  histInfo$bwFile[i]
}


##################################################################################
## other data

# 
# file_otherData <- paste(otherDataPath, "/", "sample.list", sep = "")
# 
# otherSamples <- fread(file = file_otherData, sep = "\t", header = F,
#                       stringsAsFactors = F, col.names = c("id"), data.table = F)
# 
# i <- 1
# 
# ## generate profile matrix for all samples
# foreach(i = 1:nrow(otherSamples)) %dopar% {
#   
#   sampleDataPath <- paste(otherDataPath, "/", otherSamples$id[i], sep = "")
#   file_sampleBw <- paste(sampleDataPath, "/", otherSamples$id[i], "_normalized.bw", sep = "")
#   file_profileMat <- paste(sampleDataPath, "/", otherSamples$id[i], "_normalizedMatrix.tab.gz", sep = "")
#   
#   ## 2kb - 2kb - 1kb matrix
#   bwMat <- chipmine::bigwig_profile_matrix(bwFile = file_sampleBw,
#                                            bedFile = file_genes,
#                                            signalName = otherSamples$id[i],
#                                            genes = geneSet$gene,
#                                            readLocal = FALSE,
#                                            storeLocal = TRUE,
#                                            localPath = file_profileMat)
# 
#   mat5Kb <- gsub(pattern = ".tab.gz", replacement = "_5kb.tab.gz", x = file_profileMat)
# 
#   ## 5kb - 2kb - 1kb matrix
#   bwMat <- chipmine::bigwig_profile_matrix(bwFile = file_sampleBw,
#                                            bedFile = file_genes,
#                                            signalName = otherSamples$id[i],
#                                            genes = geneSet$gene,
#                                            readLocal = FALSE,
#                                            storeLocal = TRUE,
#                                            localPath = mat5Kb,
#                                            extend = c(5000, 1000),
#                                            target_ratio = 0.25)
# 
#   matTss <- gsub(pattern = ".tab.gz", replacement = "_4kbTSS2kb.tab.gz", x = file_profileMat)
#   ## -4kb - TSS - 2kb
#   bwMat <- chipmine::bigwig_profile_matrix(bwFile = file_sampleBw,
#                                            bedFile = file_genes,
#                                            signalName = otherSamples$id[i],
#                                            genes = geneSet$gene,
#                                            readLocal = FALSE,
#                                            storeLocal = TRUE,
#                                            localPath = matTss,
#                                            extend = c(4000, 2000),
#                                            target = "tss")
#   
#   matTes <- gsub(pattern = ".tab.gz", replacement = "_2kbTES4kb.tab.gz", x = file_profileMat)
#   ## -2kb - TES - 4kb
#   bwMat <- chipmine::bigwig_profile_matrix(bwFile = file_sampleBw,
#                                            bedFile = file_genes,
#                                            signalName = otherSamples$id[i],
#                                            genes = geneSet$gene,
#                                            readLocal = FALSE,
#                                            storeLocal = TRUE,
#                                            localPath = matTes,
#                                            extend = c(2000, 4000),
#                                            target = "tes")
#   
#   
#   matTss <- gsub(pattern = ".tab.gz", replacement = "_3kbTSS3kb.tab.gz", x = file_profileMat)
#   ## -5kb - TSS - 1kb
#   bwMat <- chipmine::bigwig_profile_matrix(bwFile = file_sampleBw,
#                                            bedFile = file_genes,
#                                            signalName = otherSamples$id[i],
#                                            genes = geneSet$gene,
#                                            readLocal = FALSE,
#                                            storeLocal = TRUE,
#                                            localPath = matTss,
#                                            extend = c(3000, 3000),
#                                            target = "tss")
#   
#   
#   
#   matTes <- gsub(pattern = ".tab.gz", replacement = "_3kbTES3kb.tab.gz", x = file_profileMat)
#   ## -1kb - TES - 5kb
#   bwMat <- chipmine::bigwig_profile_matrix(bwFile = file_sampleBw,
#                                            bedFile = file_genes,
#                                            signalName = otherSamples$id[i],
#                                            genes = geneSet$gene,
#                                            readLocal = FALSE,
#                                            storeLocal = TRUE,
#                                            localPath = matTes,
#                                            extend = c(3000, 3000),
#                                            target = "tes")
#   
#   otherSamples$id[i]
# }
# 
# 
# 



##################################################################################
stopCluster(cl)

