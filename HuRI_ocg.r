library(linkcomm)
library(openxlsx)
library(igraph)

# **2020.12.09 16:38** using largest connected component for community establishment
# **2020.12.07 17:11** remove those parameters
# original
huri <- read.csv("~/Documents/INET-work/references/HuRI_binaryPPI/HuRI_Tong_withSymbol.csv", header = T)
huri <- huri[, c(5:6)]
huri_graph <- graph_from_data_frame(huri, directed = FALSE)
huri_lc <- components(huri_graph, mode = "strong")
huri_sub <- induced_subgraph(huri_graph, names(huri_lc$membership[huri_lc[1]$membership == 1]))

huri_ocg <- getOCG.clusters(as_data_frame(huri_sub)); # save.image(file = "~/Documents/INET-work/virus_network/statistic_results/HuRI_ocg.RData");
# huri_lc <- getOCG.clusters(huri); # save.image(file = "~/Documents/INET-work/virus_network/statistic_results/HuRI_lc.RData");

hi_union <- read.table("~/Documents/INET-work/references/HuRI_binaryPPI/HI-union.txt", header = F, sep = "\t")
hi_ocg <- getOCG.clusters(hi_union); # save.image(file = "~/Documents/INET-work/virus_network/statistic_results/HI_union_ocg.RData")

# **2020.11.23 22:33** use individual assay, with LinkComm default parameters
huri_assay1 <- read.xlsx("~/Documents/INET-work/references/HuRI_binaryPPI/HuRI_individual_assay.xlsx", sheet = 1)
huri_a1_graph <- graph_from_data_frame(huri_assay1[, c(3, 4)], directed = FALSE)
huri_a1_lc <- components(huri_a1_graph)
huri_a1_sub <- induced_subgraph(huri_a1_graph, names(huri_a1_lc[1]$membership))

huri_ocg1 <- getOCG.clusters(as_data_frame(huri_a1_sub)); # save.image(file = "~/Documents/INET-work/virus_network/statistic_results/HuRI_ocg_assay1.RData");
#huri_lc1 <- getLinkCommunities(huri_assay1[, c(3, 4)]); # save.image(file = "~/Documents/INET-work/virus_network/statistic_results/HuRI_lc_assay1.RData");

huri_assay2 <- read.xlsx("~/Documents/INET-work/references/HuRI_binaryPPI/HuRI_individual_assay.xlsx", sheet = 2)
huri_a2_graph <- graph_from_data_frame(huri_assay2[, c(3, 4)], directed = FALSE)
huri_a2_lc <- components(huri_a2_graph)
huri_a2_sub <- induced_subgraph(huri_a2_graph, names(huri_a2_lc[1]$membership))

huri_ocg2 <- getOCG.clusters(as_data_frame(huri_a2_sub)); # save.image(file = "~/Documents/INET-work/virus_network/statistic_results/HuRI_ocg_assay2.RData");
# huri_lc2 <- getLinkCommunities(huri_assay1[, c(3, 4)]); # save.image(file = "~/Documents/INET-work/virus_network/statistic_results/HuRI_lc_assay2.RData");

huri_assay3 <- read.xlsx("~/Documents/INET-work/references/HuRI_binaryPPI/HuRI_individual_assay.xlsx", sheet = 3)
huri_a3_graph <- graph_from_data_frame(huri_assay3[, c(3, 4)], directed = FALSE)
huri_a3_lc <- components(huri_a3_graph)
huri_a3_sub <- induced_subgraph(huri_a3_graph, names(huri_a3_lc[1]$membership))

huri_ocg3 <- getOCG.clusters(as_data_frame(huri_a3_sub)); # save.image(file = "~/Documents/INET-work/virus_network/statistic_results/HuRI_ocg_assay3.RData");
# huri_lc3 <- getLinkCommunities(huri_assay1[, c(3, 4)]); # save.image(file = "~/Documents/INET-work/virus_network/statistic_results/HuRI_lc_assay3.RData");

# it's because the above scripts are manual simultaneous running, now manual re-load
# setwd("~/Documents/INET-work/virus_network/statistic_results/")
# load("HuRI_ocg_assay1.RData"); load("HuRI_ocg_assay2.RData"); load("HuRI_ocg_assay3.RData");

wb <- createWorkbook()
addWorksheet(wb, "assay1")
addWorksheet(wb, "assay2")
addWorksheet(wb, "assay3")
addWorksheet(wb, "HuRI")
addWorksheet(wb, "HI-union")
writeData(wb, sheet = "assay1", huri_ocg1[4])
writeData(wb, sheet = "assay2", huri_ocg2[4])
writeData(wb, sheet = "assay3", huri_ocg3[4])
writeData(wb, sheet = "HuRI", huri_ocg[4])
writeData(wb, sheet = "HI-union", hi_ocg[4])
saveWorkbook(wb, "~/Documents/INET-work/virus_network/statistic_results/HuRI_assays_OCG.xlsx", overwrite = TRUE)
