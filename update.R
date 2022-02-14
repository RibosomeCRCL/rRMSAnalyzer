library(devtools)
library(roxygen2)

# load(file = "/home/hermes/Desktop/RMS/R/sysdata.rda")
# load(file = "/home/hermes/Desktop/RMS/data/mouse.methylated.rda")
# load(file = "/home/hermes/Desktop/RMS/data/mouse.suspected.rda")
# 
# #mouse.methylated <- read.csv("/home/hermes/Desktop/home/hermes/Bureau/CDD/Mouse_vs_Human_vs/Last_time/known_sites_mouse.csv")
# #mouse.suspected <- read.csv("/home/hermes/Desktop/home/hermes/Bureau/CDD/Mouse_vs_Human_vs/Last_time/suspected_sites_mouse.csv")
# names(mouse.methylated)[4] <- "NR_046235.Numbering"
# mouse.methylated$SNORD <- NA
# mouse.methylated$Mode.of.coding <- NA
# mouse.methylated$SNORD.host.gene <- NA
# #mouse.methylated <- known
# #mouse.suspected <- suspcted
# 
# all_z <- read.csv("/home/hermes/Desktop/home/hermes/Bureau/CDD/Mouse_vs_Human_vs/Other_Modification/all_zero.csv")
# all_z$Nomenclature <- paste(all_z$Sequence, all_z$NR_003279.1_28SrRNA_mouse_position, sep = "")
# 
# mouse.suspected <- dplyr::anti_join(mouse.suspected, all_z, by = "Nomenclature")
# use_data(mouse.suspected)
# 
# 
# # change name of human
# 
# names(human.methylated)[4] <- "NR_046235.Numbering"
# names(human.suspected)[4] <- "NR_046235.Numbering"
# 
# 
# use_data(wdw.probas.5s, wdw.probas.18s, wdw.probas.28s, wdw.probas.5.8s, human.suspected, human.methylated, 
#          mouse.wdw.probas.5.8s, mouse.wdw.probas.18s, mouse.wdw.probas.28s, mouse.suspected, mouse.methylated, internal=T, overwrite=T)
# 
# save(mouse.suspected, file = "/home/hermes/Desktop/RMS/data/mouse.suspected.rda")



# Update Human methylated and suspected sites


load("data/human_methylated_28S.Rdata")
load("data/human_methylated_18S.Rdata")
load("data/human_methylated_58S.Rdata")
load("data/human_suspected_18S.Rdata")
load("data/human_suspected_28S.Rdata")
load("data/human_suspected_58S.Rdata")

human_methylated <- rbind(df.methylated.58S, rbind(df.methylated.18S, df.methylated.28S))
colnames(df.suspected.18S)[9] <- "Ensembl"
human_suspected <- rbind(df.suspected.5.8S, rbind(df.suspected.18S, df.suspected.28S))

use_data(human_methylated, human_suspected, internal=F, overwrite=T)
