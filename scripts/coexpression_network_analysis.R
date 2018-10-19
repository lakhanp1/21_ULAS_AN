library(dplyr)
library(data.table)
library(tibble)
library(purrr)


rm(list = ls())

path = "E:/Chris_UM/Analysis/23_ShuhuiMinorTasks/coexpression_net"
setwd(path)


file_mat = "correlation_plot_151_datasets_SM_46cluster.tab"

pccMat = fread(file = file_mat, sep = "\t", header = T, stringsAsFactors = F, data.table = F) %>% 
  column_to_rownames("V1")



## function to get the rank for each gene's PCC values as structured list
get_gene_rank = function(pccList, names){
  l1 = structure(pccList, names = names)
  r1 = rank(-pccList, ties.method= "min")
  
  n = map2(.x = l1, .y = r1, .f = function(x, y){list(val = x, rank = y)})
  
  return(n)
}



rankData = lapply(X = pccMat, FUN = get_gene_rank, names = rownames(pccMat))

g1 = "AN0029"
g2 = "AN0024"
rg12 = rankData[[g1]][[g2]]$rank
rg21 = rankData[[g2]][[g1]]$rank

mutualRank = sqrt(rg12 * rg21)


## calculate mutual rank score
mrMat = matrix(data = NA, nrow = nrow(pccMat), ncol = ncol(pccMat), dimnames = list(rownames(pccMat), colnames(pccMat)))

for (ic in colnames(mrMat)) {
  for (ir in rownames(mrMat)) {
    rg12 = rankData[[ic]][[ir]]$rank
    rg21 = rankData[[ir]][[ic]]$rank
    
    mrMat[ic, ir] = sprintf("%.3f", sqrt(rg12 * rg21))
  }
}


##
write.csv(x = mrMat, file = "MR_rank_scores.csv", row.names = T, col.names = T)





