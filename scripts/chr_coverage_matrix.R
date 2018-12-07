library(chipmine)
library(rtracklayer)
library(GenomicRanges)
library(here)

## TF samples:
## 1) extract genome wide 1kb bin wise coverage for each TF sample
## 2) use background region bins i.e. bins which do not overlap with peak region
## for calculating summary statistics
##
## polII samples
## 1) extract genome wide 1kb bin wise coverage for each TF sample
## 2) use background region bins i.e. bins which do not overlap with gene region
## for calculating summary statistics

rm(list = ls())

##################################################################################
outPrefix <- here::here("generic_analysis", "mitochondrial_coverage", "chr_coverage_summary")


## genes to read
file_exptInfo <- here::here("data", "referenceData/sampleInfo.txt")
file_genes <- here::here("data", "referenceData/AN_genesForPolII.bed")
file_topGoMap <- "E:/Chris_UM/Database/A_Nidulans/ANidulans_OrgDb/geneid2go.ANidulans.topGO.map"
file_geneInfo <- "E:/Chris_UM/Database/A_Nidulans/A_nidulans_FGSC_A4_geneClasses.txt"
file_chrs <- here::here("data", "referenceData/A_Nidulans_chromosomes.bed")

TF_dataPath <- here::here("data", "TF_data")
polII_dataPath <- here::here("data", "polII_data")
hist_dataPath <- here::here("data", "histone_data")

chrGr <- import.bed(con = file_chrs)
seqlengths(chrGr) <- width(chrGr)

genomeBins <- tileGenome(seqlengths = seqlengths(chrGr),tilewidth = 1000, cut.last.tile.in.chrom = TRUE)
genomeBins$chr <- as.character(seqnames(genomeBins))
genomeBins$region <- paste(as.character(seqnames(genomeBins)), ":", start(genomeBins), "-", end(genomeBins), sep = "")

## import gene regions and find its overlap with genomewide 1kb bins
geneGr <- import.bed(con = file_genes)
binOpGene <- findOverlaps(query = genomeBins,
                          subject = geneGr)

## set background information to FALSE for the bins which overlap with gene regions
# mcols(genomeBins)[[ "intergenicBin" ]] <- TRUE
# mcols(genomeBins)[[ "intergenicBin" ]][ unique(queryHits(binOpGene)) ] <- FALSE
# intergenicBins <- genomeBins[genomeBins$intergenicBin]
# export.bed(object = intergenicBins, con = here::here("data", "referenceData", "intergenic_bins.bed"))

##################################################################################
tfSampleFile <- paste(TF_dataPath, "/", "sample_tfs.list", sep = "")

tfSampleList <- fread(file = tfSampleFile, sep = "\t", header = F,
                      stringsAsFactors = F, col.names = c("id"), data.table = F)


tfSampleList <- data.frame(id = c("An_kdmB_20h_HA_1", "An_sudA_20h_HA_1", "An_rpdA_20h_HA_1", "An_sntB_20h_HA_1", "An_ecoA_20h_HA_1", "An_kdmB_48h_HA_1", "An_rpdA_48h_HA_1", "An_sntB_48h_HA_1", "An_sudA_48h_HA_1", "An_ecoA_48h_HA_1"),
                           stringsAsFactors = F)


tfInfo <- get_sample_information(exptInfoFile = file_exptInfo,
                                 samples = tfSampleList$id,
                                 dataPath = TF_dataPath,
                                 matrixSource = "normalizedmatrix")

i <- 1

for (i in 1:nrow(tfInfo)) {
  ## import bigwig file
  bwGr <- rtracklayer::import.bw(con = tfInfo$bwFile[i], as = "RleList")
  
  sampleIdCol <- tfInfo$sampleId[i]
  backgroundCol <- paste("background.", sampleIdCol, sep = "")
  
  mcols(chrGr)[[ sampleIdCol ]] <- 0
  
  ## get coverage for each chromosome
  for (chr in seqinfo(chrGr)@seqnames) {
    chrScore <- sum(Views(subject = bwGr[[chr]],
                          ranges(chrGr[seqnames(chrGr) == chr]))) / width(chrGr[seqnames(chrGr) == chr])
    
    mcols(chrGr[seqnames(chrGr) == chr])[[ sampleIdCol ]]  <- chrScore
    
  }
  
  ## set the default values for sample's coverage and overlap with macs2 peak list
  mcols(genomeBins)[[ sampleIdCol ]] <- 0
  mcols(genomeBins)[[ backgroundCol ]] <- TRUE
  
  ## get binwise coverage
  for (chr in seqinfo(genomeBins)@seqnames) {
    binScore <- sum(
      Views(subject = bwGr[[chr]],
            ranges(genomeBins[seqnames(genomeBins) == chr]))
    ) / width(genomeBins[seqnames(genomeBins) == chr])
    
    mcols(genomeBins[seqnames(genomeBins) == chr])[[ sampleIdCol ]] <- binScore
  }
  
  
  if(file.exists(tfInfo$narrowpeakFile[i])){
    peaksGr <- rtracklayer::import(con = tfInfo$narrowpeakFile[i], format = "narrowPeak")
    
    ## find overlap of genome wide 1kb bins with macs2 peaks
    binOp <- findOverlaps(query = genomeBins,
                          subject = peaksGr)
    
    ## set background information to FALSE for the bins which overlap with peaks
    mcols(genomeBins)[[ backgroundCol ]][ unique(queryHits(binOp)) ] <- FALSE
  }

  

}


##################################################################################
## polII data

polIIsampleFile <- paste(polII_dataPath, "/", "sample_polII.list", sep = "")

polIISampleList <- fread(file = polIIsampleFile, sep = "\t", header = F,
                         stringsAsFactors = F, col.names = c("id"), data.table = F)

# polIISampleList = data.frame(id = c("An_untagged_20h_polII_1", "An_untagged_48h_polII_1"), stringsAsFactors = F)

polII_info <- get_sample_information(exptInfoFile = file_exptInfo,
                                     samples = polIISampleList$id,
                                     dataPath = polII_dataPath,
                                     matrixSource = "normalizedmatrix")


i <- 1

for (i in 1:nrow(polII_info)) {
  ## import bigwig file
  bwGr <- rtracklayer::import.bw(con = polII_info$bwFile[i], as = "RleList")
  
  sampleIdCol <- polII_info$sampleId[i]
  backgroundCol <- paste("background.", sampleIdCol, sep = "")
  
  mcols(chrGr)[[ sampleIdCol ]] <- 0
  
  ## get coverage for each chromosome
  for (chr in seqinfo(chrGr)@seqnames) {
    chrScore <- sum(Views(subject = bwGr[[chr]],
                          ranges(chrGr[seqnames(chrGr) == chr]))) / width(chrGr[seqnames(chrGr) == chr])
    
    mcols(chrGr[seqnames(chrGr) == chr])[[ polII_info$sampleId[i] ]]  <- chrScore
    
  }
  
  
  ## set the default values for sample's coverage and overlap with macs2 peak list
  mcols(genomeBins)[[ sampleIdCol ]] <- 0
  mcols(genomeBins)[[ backgroundCol ]] <- TRUE
  
  ## get binwise coverage
  for (chr in seqinfo(genomeBins)@seqnames) {
    binScore <- sum(
      Views(subject = bwGr[[chr]],
            ranges(genomeBins[seqnames(genomeBins) == chr]))
    ) / width(genomeBins[seqnames(genomeBins) == chr])
    
    mcols(genomeBins[seqnames(genomeBins) == chr])[[ sampleIdCol ]] <- binScore
  }
  
  ## set background information to FALSE for the bins which overlap with gene regions
  mcols(genomeBins)[[ backgroundCol ]][ unique(queryHits(binOpGene)) ] <- FALSE
  
}

##################################################################################
## store TF and polII summary
df <- as.data.frame(mcols(genomeBins))

summaryList <- sapply(X = c(tfInfo$sampleId, polII_info$sampleId), FUN = function(x){
  backgroundCol <- paste("background.", x, sep = "")
  quantile(x = df[[x]], seq(0, 1, 0.1))
  nBackground <- length(which(df[[backgroundCol]]))
  
  ## genome level summary
  genomicMedian <- median(df[[x]][ df$chr != "mito_A_nidulans_FGSC_A4" ])
  genomicMean <- mean(df[[x]][ df$chr != "mito_A_nidulans_FGSC_A4" ])
  
  ## genomic background (! macs2 peaks) region summary
  genomicBgMedian <- median(df[[x]][ df[[backgroundCol]] & df$chr != "mito_A_nidulans_FGSC_A4" ])
  genomicBgMean <- mean(df[[x]][ df[[backgroundCol]] & df$chr != "mito_A_nidulans_FGSC_A4" ])
  
  ## mitochondrial genome summary
  mtMedian <- median(df[[x]][ df$chr == "mito_A_nidulans_FGSC_A4" ])
  mtMean <- mean(df[[x]][ df$chr == "mito_A_nidulans_FGSC_A4" ])
  
  return(c(
    "nBackground" = nBackground,
    "genomicMedian" = genomicMedian,
    "genomicMean" = genomicMean,
    "genomicBgMedian" = genomicBgMedian,
    "genomicBgMean" = genomicBgMean,
    "mtMedian" = mtMedian,
    "mtMean" = mtMean))
  
}, USE.NAMES = TRUE, simplify = FALSE)


summaryDf <- as.data.frame(do.call("rbind", summaryList)) %>% 
  tibble::rownames_to_column(var = "sampleId") %>% 
  dplyr::mutate(type = if_else(condition = grepl(pattern = "_polII_", x = sampleId, fixed = TRUE),
                               true = "polII", false = "TF")) %>% 
  dplyr::select(sampleId, type, everything())

fwrite(x = summaryDf, file = paste(outPrefix, "_tf_polII.tab", sep = ""), sep = "\t", quote = FALSE, col.names = T)

##################################################################################
## histone data

file_histData <- paste(hist_dataPath, "/", "histone_sample.list", sep = "")

histSamples <- fread(file = file_histData, sep = "\t", header = F,
                     stringsAsFactors = F, col.names = c("id"), data.table = F)

histSamples <- data.frame(id = c("An_H3_20h_HIST_1", "An_H3_48h_HIST_1", "An_H3_cclA_del_20h_HIST_1",
                                 "An_H3_cclA_del_48h_HIST_1"),
                          stringsAsFactors = F)

histInfo <- get_sample_information(exptInfoFile = file_exptInfo,
                                   samples = histSamples$id,
                                   dataPath = hist_dataPath,
                                   matrixSource = "normalizedmatrix")


geneGr <- import.bed(con = file_genes)

tesGr <- get_TES(gr = geneGr, up = 100, down = 100)

i <- 1

for (i in 1:nrow(histInfo)) {
  ## import bigwig file
  bwGr <- rtracklayer::import.bw(con = histInfo$bwFile[i], as = "RleList")

  sampleIdCol <- histInfo$sampleId[i]

  mcols(chrGr)[[ sampleIdCol ]] <- 0
  
  for (chr in seqinfo(chrGr)@seqnames) {
    chrScore <- sum(Views(subject = bwGr[[chr]],
                          ranges(chrGr[seqnames(chrGr) == chr]))) / width(chrGr[seqnames(chrGr) == chr])
    
    mcols(chrGr[seqnames(chrGr) == chr])[[ histInfo$sampleId[i] ]]  <- chrScore
    
  }
  
  
  ## set the default values for sample's coverage and overlap with macs2 peak list
  mcols(tesGr)[[ sampleIdCol ]] <- 0

  ## get binwise coverage
  for (chr in seqinfo(tesGr)@seqnames) {
    binScore <- sum(
      Views(subject = bwGr[[chr]],
            ranges(tesGr[seqnames(tesGr) == chr]))
    ) / width(tesGr[seqnames(tesGr) == chr])
    
    mcols(tesGr[seqnames(tesGr) == chr])[[ sampleIdCol ]] <- binScore
  }
  

}



##################################################################################

df <- data.table::melt.data.table(data = as.data.table(mcols(chrGr)), 
                                  id.vars = "name",
                                  variable.name = "sampleId",
                                  value.name = "coverage") %>% 
  data.table::dcast.data.table(formula = sampleId ~ name,
                               value.var = "coverage")


fwrite(x = df, file = "chr_coverage.tab", sep = "\t", quote = FALSE, col.names = T)






melt.data.table()



