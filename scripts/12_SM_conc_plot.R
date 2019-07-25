library(data.table)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggforce)


## This script plots the SM concentration pie chart scatter plot for various strains and SM

rm(list = ls())

path <- "E:/Chris_UM/Analysis/21_CL2017_ChIPmix_ULAS_MIX/ULAS_AN/kdmA_analysis/SM_plot"

setwd(path)

smOrder <- c("Penicillin G", "Sterigmatocystin", "Norsolorinic acid", "Averufin", "Austinol", "Dehydroaustinol", "Emericellamide A", "Emericellamide C", "Citreorosein", "Emodin", "FK 9775 A", "FK 9775 B", "Orsellinic acid", "Versicolorin A", "Versicolorin C", "Averantin", "Nidurufin", "seco-Sterigmatocystin", "Emericellamide E", "Chrysophanol", "Dichlordiaportin", "Endocrocin", "Ilicicolin B", "Iso-Rhodoptilometrin")

strainOrder <- c("WT", "AN KdmA*Delta", "AN CclA*Delta", "AN EcmB*Delta", "AN RstB*Delta", "AN McmA-teton", "AN KdmA*Delta CclA*Delta", "AN KdmA*Delta EcmB*Delta", "AN KdmA*Delta RstB*Delta", "AN  RstB*Delta CclA*Delta", "AN RstB*Delta EcmB*Delta", "AN CclA*Delta EcmB*Delta")

##################################################################################

file_mat <- "SM_chromatograph_data.tab"

mat <- fread(file = file_mat, sep = "\t", header = T, stringsAsFactors = F) %>% 
  dplyr::mutate_if(.predicate = is.integer, .funs = as.numeric)

smOrder <- gsub(pattern = "\\s+", replacement = "~", x = smOrder, perl = T)

mat$Strain <- gsub(pattern = "\\s+", replacement = "~", x = mat$Strain, perl = T)
colnames(mat) <- gsub(pattern = "\\s+", replacement = "~", x = colnames(mat), perl = T)

## calculate the average SM concentration for each strain
mat2 <- melt.data.table(data = as.data.table(mat),
                        id.vars = c("Strain", "rowLabel", "time"),
                        variable.factor = FALSE, value.factor = FALSE,
                        variable.name = "sm",
                        value.name = "amount") %>% 
  dplyr::group_by(Strain, sm) %>% 
  dplyr::summarise(avgAmount = as.numeric(sprintf("%.2f", mean(amount)))) %>% 
  dplyr::ungroup()


if (all(levels(mat2$sm) %in% smOrder)) {
  mat2$sm <- factor(x = mat2$sm, levels = smOrder)
}

mat2$Strain <- factor(mat2$Strain, levels = unique(mat$Strain))



## find the maximum SM amount for each SM among different strains
mat3 <- dplyr::group_by(mat2, sm) %>% 
  dplyr::mutate(maxConc = max(avgAmount)) %>% 
  dplyr::mutate(blank = maxConc - avgAmount) %>% 
  dplyr::ungroup()


## convert data frame to long format by merging avgAmount and blank columns
mat4 <- melt.data.table(data = as.data.table(mat3),
                        id.vars = c("Strain", "sm", "maxConc"),
                        measure.vars = c("avgAmount", "blank"), variable.name = "grp", value.name = "SM_level") %>% 
  dplyr::group_by(sm, Strain) %>% 
  dplyr::arrange(sm, Strain, grp, .by_group = TRUE) %>% 
  dplyr::ungroup()


valueMat <- dplyr::filter(mat4, grp == "avgAmount")


p <- ggplot() +
  theme_no_axes() +
  coord_fixed() +
  geom_arc_bar(data = mat4,
               mapping = aes(x0 = 0, y0 = 0, r0 = 0, r = 1, amount = SM_level, fill = grp),
               stat = "pie", size = 1, linetype = "solid") +
  # geom_text(data = valueMat,
  #           mapping = aes(x = 0, y = -1.4, label = SM_level), color = "red", vjust = 0) +
  scale_fill_manual(values = c("avgAmount" = "#66c2a5", "blank" = "white")) +
  facet_grid(sm ~ Strain, labeller = label_parsed, switch = "y") +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 30, angle = 90, hjust = 0),
        strip.text.y = element_text(size = 30, angle = 180, hjust = 1),
        # panel.spacing = unit(1.5, "lines"),
        # rect = element_blank(),
        legend.position="none"
  )

png(filename = "SM_levels_2.png", width=6000, height=8000, res = 250)
p
dev.off()

pdf(file = "SM_levels_2.pdf", width = 20, height = 30)
p
dev.off()

svg(filename = "SM_levels.svg", width = 20, height = 30)
p
dev.off()




##################################################################################


file_mat <- "test.tab"


mat <- fread(file = file_mat, sep = "\t", header = T, stringsAsFactors = F)

mat$amount <- as.numeric(mat$amount)

mat2 <- dplyr::group_by(mat, sm) %>% 
  dplyr::mutate(grpMax = max(amount)) %>% 
  dplyr::mutate(blank = grpMax - amount) %>% 
  dplyr::ungroup()



mat3 <- melt.data.table(data = as.data.table(mat2),
                        id.vars = c("sm", "strain", "grpMax"), measure.vars = c("amount", "blank")) %>% 
  dplyr::group_by(sm, strain) %>% 
  dplyr::arrange(sm, strain, variable, .by_group = TRUE) %>%
  dplyr::ungroup()

mat3$variable <- factor(mat3$variable, levels = c("amount", "blank"))

mat4 <- dplyr::filter(mat3, variable == "amount")

# write.table(x = mat3, file = "clipboard", sep = "\t", col.names = T, row.names = F)

p <- ggplot() +
  theme_no_axes() +
  coord_fixed() +
  geom_arc_bar(data = mat3,
               mapping = aes(x0 = 0, y0 = 0, r0 = 0, r = 1, amount = value, fill = variable),
               stat = "pie", size = 2, linetype = "solid") +
  geom_text(data = mat4,
            mapping = aes(x = 0, y = -1.3, label = value)) +
  scale_fill_manual(values = c("amount" = "#fdb462", "blank" = "white")) +
  facet_grid(sm ~ strain, labeller = label_parsed, switch = "y") +
  theme(strip.text.x = element_text(size = 20, angle = 90, hjust = 0),
        strip.text.y = element_text(size = 20, angle = 180, hjust = 0),
        panel.spacing = unit(1, "lines"),
        # rect = element_blank(),
        strip.background = element_blank()
  )


png(filename = "test.png", width=4000, height=4000, res = 450)
p
dev.off()

##################################################################################



