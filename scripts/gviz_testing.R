
library(rtracklayer)
library(Gviz)


rm(list = ls())


path <-  "E:/Chris_UM/Analysis/21_CL2017_ChIPmix_ULAS_MIX/ULAS_AN/data/TF_data/An_kdmB_20h_HA_1"
setwd(path)


file_bw <- "An_kdmB_20h_HA_1_normalized.bw"
file_peaks <- "An_kdmB_20h_HA_1_narrow_peaks.narrowPeak"
file_genes <- "E:/Chris_UM/Analysis/21_CL2017_ChIPmix_ULAS_MIX/ULAS_AN/data/referenceData/AN_genesForPolII.bed"
file_genes <- "E:/Chris_UM/Analysis/21_CL2017_ChIPmix_ULAS_MIX/ULAS_AN/data/referenceData/A_Nidulans_chromosomes.bed"



bwgr <- import.bw(con = file_bw, as = "RleList")
bedGr <- import.bed(con = file_genes)

bedGr$score <- 0

## calculate sum(read depth at each position for gene) / (gene length)
sum( Views(subject = bwgr$ChrIII_A_nidulans_FGSC_A4, ranges(head(bedGr, 1))) ) / width(head(bedGr, 1))

sampleScore <- numeric()

chr <- seqinfo(bedGr)@seqnames[1]

for (chr in seqinfo(bedGr)@seqnames) {
  chrScore <- sum(Views(subject = bwgr[[chr]],
                        ranges(bedGr[seqnames(bedGr) == chr]))) / width(bedGr[seqnames(bedGr) == chr])
  
  names(chrScore) <- bedGr[seqnames(bedGr) == chr]$name
  
  bedGr[seqnames(bedGr) == chr]$score <- chrScore
  # mcols(bedGr[seqnames(bedGr) == chr])[["score"]]  <- chrScore
  
  sampleScore <- append(sampleScore, chrScore)
}

cor(expressionData$An_untagged_48h_polII_1, sampleScore)
plot(expressionData$An_untagged_48h_polII_1, sampleScore)

sum( Views(subject = bwgr$ChrIII_A_nidulans_FGSC_A4, ranges(head(bedGr))) ) / width(head(bedGr))


# dat <- sin(seq(pi, 10*pi, len = 500))
# 
# dTrack.big <- DataTrack(
#   start = seq(1, 1e+05, len = 500),
#   width = 15,
#   chromosome = "chrX",
#   genome = "hg19",
#   name = "sinus",
#   data = sin(seq(pi, 5*pi, len = 500))*runif(500, 0.5, 1.5))
# 
# plotTracks(dTrack.big, type = "hist")
# 
# plotTracks(dTrack.big, type = "hist", window = "fixed", windowSize = 1000)


options(ucscChromosomeNames=FALSE)

dt <- DataTrack(range = bwgr, genome = "ANidulans", name = "An_cclA_20h_HA_1")

chromosome(dt)
seqinfo(dt)
seqlevels(dt)

plotTracks(trackList = dt,
           type = "heatmap",
           window = "fixed", windowSize = 5000,
           gradient = c("black", "yellow","green", "red"),
           aggregation = function(x){return(quantile(x, 0.99, names = FALSE, na.rm = TRUE))}
           )




##################################################################################

library(IRanges)

ir1 <- IRanges(start=1:10, width=10:1)
ir2 <- IRanges(start=1:10, end=11)
ir3 <- IRanges(end=11, width=10:1)
ir <- IRanges(c(1, 8, 14, 15, 19, 34, 40), width=c(12, 6, 6, 15, 6, 2, 7))

plotRanges <- function(x, xlim=x, main=deparse(substitute(x)),
                       col="black", sep=0.5, ...)
{
  height <- 1
  if (is(xlim, "IntegerRanges"))
    xlim <- c(min(start(xlim)), max(end(xlim)))
  bins <- disjointBins(IRanges(start(x), end(x) + 1))
  plot.new()
  plot.window(xlim, c(0, max(bins) * (height + sep)))
  ybottom <- bins *(sep + height) - height
  rect(start(x)-0.5, ybottom, end(x)+0.5, ybottom + height, col=col, ...)
  title(main)
  axis(1)
}

plotRanges(ir)

lambda <- c(rep(0.001, 4500), seq(0.001, 10, length=500), seq(10, 0.001, length=500))
xVector <- rpois(1e7, lambda)
yVector <- rpois(1e7, lambda[c(251:length(lambda), 1:250)])
xRle <- Rle(xVector)
yRle <- Rle(yVector)

cov <- coverage(ir)
cov <- as.vector(cov)
mat <- cbind(seq_along(cov)-0.5, cov)
d <- diff(cov) != 0
mat <- rbind(cbind(mat[d,1]+1, mat[d,2]), mat)
mat <- mat[order(mat[,1]),]
lines(mat, col="red", lwd=4)


plotRanges(disjoin(ir))

xViews <- Views(xRle, xRle >= 1)
xViews <- slice(xRle, 1)

head(viewSums(xViews), 10)


