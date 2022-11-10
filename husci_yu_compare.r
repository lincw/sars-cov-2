# comparision between HuSCI and SARS-CoV-2 contactome from Haiyuan Yu
# 11.08.2022 **17:11**
# Lin Chung-wen

library(igraph)
library(openxlsx)
library(RCy3)

husci <- read.csv('~/Documents/INET-work/virus_network/final_list/HuSCI_PPI.csv')
husci_source <- read.xlsx("~/Documents/INET-work/virus_network/manuscript/Nature_Biotechnology/66004_0_data_set_567160_r7rvsn.xlsx", sheet = 3, start = 4)
husci <- husci[!apply(husci == "", 1, all), ]
hai <- read.xlsx('~/Documents/INET-work/virus_network/references/PPIs/Zhou_Y2HandAPMS/TableS1.xlsx', na.string = '')
hai$Bait <- toupper(hai$Bait)
hai$'Prey.Symbol' <- toupper(hai$'Prey.Symbol')

husci_g <- graph_from_data_frame(husci, directed = F)
E(husci_g)$source <- ifelse(husci_source$HuSCIHIS3 == 1, "HIS3", "GFP")

hai_g_main <- graph_from_data_frame(hai[, c(1, 5)], directed = F)
hai_g_y2h <- graph_from_data_frame(hai[!is.na(hai$Source_Y2H), c(1, 5)], directed = F)
hai_g_apms <- graph_from_data_frame(hai[!is.na(hai$Source_APMS), c(1, 5)], directed = F)

h_y2h_inter <- husci_g %s% hai_g_y2h
to_delete <- which(degree(h_y2h_inter) == 0)
h_y2h_inter_fin <- delete.vertices(h_y2h_inter, to_delete)
createNetworkFromIgraph(h_y2h_inter_fin, title = "intersection_y2h", collection = 'Merged')

h_main_inter <- husci_g %s% hai_g_main
to_delete <- which(degree(h_main_inter) == 0)
h_main_inter_fin <- delete.vertices(h_main_inter, to_delete)
createNetworkFromIgraph(h_main_inter_fin, title = "intersection_main", collection = 'Merged')
