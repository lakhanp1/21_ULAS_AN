library(chipmine)
library(org.Anidulans.FGSCA4.eg.db)
library(TxDb.Anidulans.AspGD.GFF)
library(foreach)
library(doParallel)
library(here)

rm(list = ls())

# cl <- makeCluster(4) #not to overload your computer
# registerDoParallel(cl)  

##################################################################################

file_exptInfo <- here::here("data", "referenceData/sampleInfo.txt")
geneCdsFile <- "E:/Chris_UM/Database/A_Nidulans/A_nidulans_FGSC_A4_version_s10-m04-r03_CDS_Unique.bed"
file_genes <- here::here("data", "referenceData/AN_genesForPolII.bed")
orgDb <- org.Anidulans.FGSCA4.eg.db
txDb <- TxDb.Anidulans.AspGD.GFF

TF_dataPath <- here::here("data", "TF_data")
polII_dataPath <- here::here("data", "polII_data")
hist_dataPath <- here::here("data", "histone_data")
other_dataPath <- here::here("data", "other_data")

geneSet <- data.table::fread(file = file_genes, header = F,
                             col.names = c("chr", "start", "end", "geneId", "score", "strand")) %>%
  dplyr::select(geneId)

geneDesc <- AnnotationDbi::select(x = orgDb, keys = geneSet$geneId, columns = "DESCRIPTION", keytype = "GID")

geneSet <- dplyr::left_join(x = geneSet, y = geneDesc, by = c("geneId" = "GID"))

# ##################################################################################
# ## histone data
# 
# file_histData <- paste(hist_dataPath, "/", "sample_histone.list", sep = "")
# 
# 
# histSamples <- readr::read_tsv(file = file_histData, col_names = c("id"),  comment = "#") %>%
#   as.data.frame()
# 
# # histSamples <- data.frame(id = c("An_H3_20h_HIST_1", "An_H3_48h_HIST_1", "An_H3_cclA_del_20h_HIST_1",
# #                                  "An_H3_cclA_del_48h_HIST_1"),
# #                           stringsAsFactors = F)
# 
# histInfo <- get_sample_information(exptInfoFile = file_exptInfo,
#                                    samples = histSamples$id,
#                                    dataPath = hist_dataPath,
#                                    matrixSource = "normalizedmatrix")
# 
# i <- 1
# 
# ## process all polII expression matrix
# foreach(i = 1:nrow(histInfo),
#         .packages = c("chipmine"),
#         .verbose = TRUE) %dopar% {
# 
# 
#           # ## 2kb - 2kb - 1kb matrix
#           # bwMat <- chipmine::bigwig_profile_matrix(bwFile = histInfo$bwFile[i],
#           #                                          bedFile = file_genes,
#           #                                          signalName = histInfo$sampleId[i],
#           #                                          genes = geneSet$gene,
#           #                                          readLocal = FALSE,
#           #                                          storeLocal = TRUE,
#           #                                          localPath = histInfo$matFile[i])
# 
#           mat5Kb <- gsub(pattern = ".tab.gz", replacement = "_5kb.tab.gz", x = histInfo$matFile[i])
# 
#           ## 5kb - 2kb - 1kb matrix
#           bwMat <- chipmine::bigwig_profile_matrix(bwFile = histInfo$bwFile[i],
#                                                    bedFile = file_genes,
#                                                    signalName = histInfo$sampleId[i],
#                                                    genes = geneSet$gene,
#                                                    readLocal = FALSE,
#                                                    storeLocal = TRUE,
#                                                    localPath = mat5Kb,
#                                                    extend = c(5000, 1000),
#                                                    target_ratio = 0.25)
# 
# 
# 
#           # matTss <- gsub(pattern = ".tab.gz", replacement = "_3kbTSS3kb.tab.gz", x = histInfo$matFile[i])
#           # ## -3kb - TSS - 3kb
#           # bwMat <- chipmine::bigwig_profile_matrix(bwFile = histInfo$bwFile[i],
#           #                                          bedFile = file_genes,
#           #                                          signalName = histInfo$sampleId[i],
#           #                                          genes = geneSet$gene,
#           #                                          readLocal = FALSE,
#           #                                          storeLocal = TRUE,
#           #                                          localPath = matTss,
#           #                                          extend = c(3000, 3000),
#           #                                          target = "tss")
#           #
#           #
#           #
#           # matTes <- gsub(pattern = ".tab.gz", replacement = "_3kbTES3kb.tab.gz", x = histInfo$matFile[i])
#           # ## -3kb - TES - 3kb
#           # bwMat <- chipmine::bigwig_profile_matrix(bwFile = histInfo$bwFile[i],
#           #                                          bedFile = file_genes,
#           #                                          signalName = histInfo$sampleId[i],
#           #                                          genes = geneSet$gene,
#           #                                          readLocal = FALSE,
#           #                                          storeLocal = TRUE,
#           #                                          localPath = matTes,
#           #                                          extend = c(3000, 3000),
#           #                                          target = "tes")
# 
#           histInfo$bwFile[i]
#         }


# ##################################################################################
# ## other data
# 
# 
# file_otherData <- paste(other_dataPath, "/", "sample.list", sep = "")
# 
# otherSamples <- fread(file = file_otherData, sep = "\t", header = F,
#                       stringsAsFactors = F, col.names = c("id"), data.table = F)
# 
# i <- 1
# 
# ## generate profile matrix for all samples
# foreach(i = 1:nrow(otherSamples),
#         .packages = c("chipmine")) %dopar% {
#           
#           sampleDataPath <- paste(other_dataPath, "/", otherSamples$id[i], sep = "")
#           file_sampleBw <- paste(sampleDataPath, "/", otherSamples$id[i], "_normalized.bw", sep = "")
#           file_profileMat <- paste(sampleDataPath, "/", otherSamples$id[i], "_normalizedMatrix.tab.gz", sep = "")
#           
#           # ## 2kb - 2kb - 1kb matrix
#           # bwMat <- chipmine::bigwig_profile_matrix(bwFile = file_sampleBw,
#           #                                          bedFile = file_genes,
#           #                                          signalName = otherSamples$id[i],
#           #                                          genes = geneSet$gene,
#           #                                          readLocal = FALSE,
#           #                                          storeLocal = TRUE,
#           #                                          localPath = file_profileMat)
#           
#           # mat5Kb <- gsub(pattern = ".tab.gz", replacement = "_5kb.tab.gz", x = file_profileMat)
#           # 
#           # ## 5kb - 2kb - 1kb matrix
#           # bwMat <- chipmine::bigwig_profile_matrix(bwFile = file_sampleBw,
#           #                                          bedFile = file_genes,
#           #                                          signalName = otherSamples$id[i],
#           #                                          genes = geneSet$gene,
#           #                                          readLocal = FALSE,
#           #                                          storeLocal = TRUE,
#           #                                          localPath = mat5Kb,
#           #                                          extend = c(5000, 1000),
#           #                                          target_ratio = 0.25)
#           
#           
#           #
#           # matTss <- gsub(pattern = ".tab.gz", replacement = "_3kbTSS3kb.tab.gz", x = file_profileMat)
#           # ## -3kb - TSS - 3kb
#           # bwMat <- chipmine::bigwig_profile_matrix(bwFile = file_sampleBw,
#           #                                          bedFile = file_genes,
#           #                                          signalName = otherSamples$id[i],
#           #                                          genes = geneSet$gene,
#           #                                          readLocal = FALSE,
#           #                                          storeLocal = TRUE,
#           #                                          localPath = matTss,
#           #                                          extend = c(3000, 3000),
#           #                                          target = "tss")
#           #
#           #
#           #
#           # matTes <- gsub(pattern = ".tab.gz", replacement = "_3kbTES3kb.tab.gz", x = file_profileMat)
#           # ## -3kb - TES - 3kb
#           # bwMat <- chipmine::bigwig_profile_matrix(bwFile = file_sampleBw,
#           #                                          bedFile = file_genes,
#           #                                          signalName = otherSamples$id[i],
#           #                                          genes = geneSet$gene,
#           #                                          readLocal = FALSE,
#           #                                          storeLocal = TRUE,
#           #                                          localPath = matTes,
#           #                                          extend = c(3000, 3000),
#           #                                          target = "tes")
#           
#           otherSamples$id[i]
#         }
# 
# 




##################################################################################
# stopCluster(cl)

