library(chipmine)
library(foreach)
library(doParallel)

rm(list = ls())


path <- "E:/Chris_UM/Analysis/21_CL2017_ChIPmix_ULAS_MIX/ULAS_AN"
setwd(path)


##################################################################################


file_exptInfo <-"E:/Chris_UM/Analysis/21_CL2017_ChIPmix_ULAS_MIX/ULAS_AN/data/referenceData/sampleInfo.txt"
TF_dataPath <- "E:/Chris_UM/Analysis/21_CL2017_ChIPmix_ULAS_MIX/ULAS_AN/data/TF_data"
polII_dataPath <- "E:/Chris_UM/Analysis/21_CL2017_ChIPmix_ULAS_MIX/ULAS_AN/data/polII_data"
file_genes <- "E:/Chris_UM/Analysis/21_CL2017_ChIPmix_ULAS_MIX/ULAS_AN/data/referenceData/AN_genesForPolII.bed"

file_tfSamples <- paste(TF_dataPath, "/", "sample_tf_macs2.list", sep = "")



set.seed(20)
doClustering <- TRUE
clusters <- 7
tfYlim <- 0.996              ##0.999

geneFilter <- c("AN5245", "AN3245")

cl <- makeCluster(4) #not to overload your computer
registerDoParallel(cl)


##################################################################################

geneSet <- data.table::fread(file = file_genes, header = F,
                             col.names = c("chr", "start", "end", "name", "score", "strand")) %>% 
  dplyr::mutate(length = end - start) %>% 
  dplyr::filter(! name %in% geneFilter)


tfSampleList <- fread(file = file_tfSamples, sep = "\t", header = F,
                      stringsAsFactors = F, col.names = c("id"), data.table = F)


# tfSampleList <- data.frame(id = c("An_ecoA_20h_HA_1", "An_ecoA_48h_HA_1", "An_kdmB_20h_HA_1", "An_kdmB_48h_HA_1", "An_rpdA_20h_HA_1", "An_rpdA_48h_HA_1", "An_sntB_20h_HA_1", "An_sntB_48h_HA_1", "An_kdmB_20h_HA_2", "An_kdmB_48h_HA_2", "An_rpdA_20h_HA_2", "An_rpdA_48h_HA_2", "An_sntB_20h_HA_2", "An_sntB_48h_HA_2", "An_ecoA_kdmB_del_20h_HA_1", "An_ecoA_kdmB_del_48h_HA_1", "An_rpdA_kdmB_del_20h_HA_1", "An_rpdA_kdmB_del_48h_HA_1", "An_sntB_kdmB_del_20h_HA_1", "An_sntB_kdmB_del_48h_HA_1", "An_ecoA_20h_HA_2", "An_ecoA_48h_HA_2", "An_ecoA_kdmB_del_20h_HA_2", "An_ecoA_kdmB_del_48h_HA_2", "An_rpdA_kdmB_del_20h_HA_2", "An_rpdA_kdmB_del_48h_HA_2", "An_sntB_kdmB_del_20h_HA_2", "An_sntB_kdmB_del_48h_HA_2", "An_ecoA_sntB_del_20h_HA_2", "An_ecoA_sntB_del_48h_HA_2", "An_kdmB_laeA_del_20h_HA_1", "An_kdmB_laeA_del_48h_HA_1", "An_laeA_kdmB_del_20h_HA_1", "An_laeA_kdmB_del_48h_HA_1", "An_sudA_kdmB_del_20h_HA_1", "An_sudA_kdmB_del_48h_HA_1"))

tf_info <- get_sample_information(exptInfoFile = file_exptInfo,
                                  samples = tfSampleList$id,
                                  dataPath = TF_dataPath,
                                  matrixSource = "deeptools")

# i <- 20

foreach(i = 1:nrow(tf_info),
        .packages = c("chipmine")) %dopar% {
          ## read the profile matrix
          mat1 <- chipmine::import_profile_from_file(
            file = tf_info$matFile[i],
            source = "deeptools",
            signalName = tf_info$sampleId[i],
            selectGenes = geneSet$name)
          
          ## check the distribution in data
          quantile(mat1, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)
          # col_fun <- colorRamp2(quantile(mat1, c(0.50, 0.995), na.rm = T), c("white", "red"))
          
          km <- chipmine::profile_matrix_kmeans(
            mat = mat1,
            km = clusters,
            clustFile = tf_info$clusterFile[i],
            name = tf_info$sampleId[i])
          
          
          cat(as.character(Sys.time()), "Done...", tf_info$sampleId[i], "\n\n")
          
        }



stopCluster(cl)







