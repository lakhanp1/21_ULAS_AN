
library(ggplot2)
library(dplyr)


## This script compare the expression value difference for old polII data and new polII data 

rm(list = ls())



path = "E:/Chris_UM/Analysis/21_CL2017_ChIPmix_ULAS_MIX/ULAS_AN"
setwd(path)

oldSample = "An_rstB_del_20h_polII"
newSample = "An_rstB_del_20h_polII_1"


file_oldExp = paste("oldBackup/polII_data/", oldSample, "/", oldSample, "_normalizedExpression.tab", sep = "")
oldExpCol = paste("is_expressed.", oldSample, sep = "")


file_newExp = paste("polII_data/", newSample, "/", newSample, "_normalizedExpression.tab", sep = "")
newExpCol = paste("is_expressed.", newSample, sep = "")


oldExp = fread(input = file_oldExp, sep = "\t", header = T, stringsAsFactors = T)
newExp = fread(input = file_newExp, sep = "\t", header = T, stringsAsFactors = T)

expData = left_join(x = oldExp, y = newExp, by = c("gene")) %>% 
  dplyr::mutate(ratio = log2(!! as.name(oldSample) / !! as.name(newSample)))



p = ggplot() +
  geom_freqpoly(data = expData, mapping = aes(x = ratio, color = "black"), breaks = seq(-1, 1, 0.04), size = 1) +
  geom_freqpoly(data = dplyr::filter(expData, !!as.name(oldExpCol) == TRUE), 
                mapping = aes(x = ratio, color = "red"), breaks = seq(-1, 1, 0.04), size = 1) +
  geom_freqpoly(data = dplyr::filter(expData, !!as.name(newExpCol) == TRUE),
                mapping = aes(x = ratio, color = "blue"), breaks = seq(-1, 1, 0.04), size = 1) +
  geom_vline(xintercept = 0, color = "green", linetype = "dashed", size = 1) + 
  scale_color_manual(
    name = newSample,
    values = c("black" = "black", "red" = "red", "blue" = "blue"),
    breaks = c("black", "red", "blue"),
    labels = c("black" = "all genes",
               "red" = "top 10% expressed in old seq",
               "blue" = "top 10% expressed in old+new seq")) + 
  theme_bw() +
  labs(title = "Difference between polII expression values between old and (old+new) sequencing data",
       x = "log2(old/(old+new))",
       y = "# genes") + 
  theme(legend.position = c(0.99, 0.99),
        legend.justification = c(1, 1),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20))


print(p)


outFileName = paste("old_vs_new_data/", newSample, ".png", sep = "")
png(filename = outFileName, width=1000, height=1000, res = 120)
print(p)
dev.off()



################################################################################
## with faceting
















