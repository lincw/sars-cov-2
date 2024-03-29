# COVID-19 GWAS + 1 community detection
# Lin Chung-wen

######
# loading packages
library(linkcomm)

######
# load dataset
gwas <- read.csv("~/Documents/INET-work/virus_network/statistic_results/GWAS/HuSCI_and_HuRI_+GWAS_interaction_v2.csv", header = T)
gwas <- gwas[, c(9, 10)]
gwas_graph <- graph_from_data_frame(gwas, directed = FALSE)
gwas_graph <- simplify(gwas_graph, remove.loops = FALSE)
gwas_strong <- components(gwas_graph, mode = "strong")
gwas_graph_sub <- induced_subgraph(gwas_graph, names(gwas_strong$membership[gwas_strong[1]$membership == 1])) # only the largest connected component was used for the OCG community detection

######
# community detection
# 1. edge betweenness
gwas_eb_whole <- cluster_edge_betweenness(gwas_graph, directed = FALSE)
gwas_eb <- cluster_edge_betweenness(gwas_graph_sub, directed = FALSE)
# 2. edge similarity
gwas_strong_df <- as_data_frame(gwas_graph_sub)
gwas_ocg <- getOCG.clusters(gwas_strong_df, init.class.sys = 1, cent.class.sys = 0) # with Maximal cliques and final class system
gwas_ocg2 <- getOCG.clusters(gwas_strong_df, cent.class.sys = 0) # with Centered cliques

######
# save R data
save.image("~/Documents/INET-work/virus_network/statistic_results/GWAS/GWAS_community_v2.RData")
