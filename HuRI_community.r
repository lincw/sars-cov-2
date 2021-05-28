# network community detection for HuRI
# Lin Chung-wen
# 18.05.2021 (adapted from "coronavirus_subnetwork.rmd", 02.11.2020)

library(linkcomm)
# library(igraph)

huri <- read.csv("~/Documents/INET-work/references/HuRI_binaryPPI/HuRI_Tong_withSymbol.csv", header = T)
huri <- huri[, c(5:6)]
huri_g <- graph_from_data_frame(huri, directed = FALSE)
huri_g <- simplify(huri_g, remove.loops = FALSE)
huri_g_strong <- components(huri_g, mode = "strong")
huri_g_sub <- induced_subgraph(huri_g, names(huri_g_strong$membership[huri_g_strong[1]$membership == 1]))

# walk trap
cw <- cluster_walktrap(huri_g)
cw_sub <- cluster_walktrap(huri_g_sub)

# Community structure via greedy optimization of modularity
cfg <- cluster_fast_greedy(huri_g)
cfg_sub <- cluster_fast_greedy(huri_g_sub)

# propagating labels
clp <- cluster_label_prop(huri_g)
clp_sub <- cluster_f(huri_g_sub)

# Community structure detecting based on the leading eigenvector of the community matrix
cle <- cluster_leading_eigen(huri_g)

# Finding community structure by multi-level optimization of modularity
cl <- cluster_louvain(huri_g)

# edge betweenness
ceb <- cluster_edge_betweenness(huri_g, directed = FALSE)

# !!!Optimal community structure
co <- cluster_optimal(huri_g)
#! At optimal_modularity.c:85 : GLPK is not available, Unimplemented function call

# spin-glass
cls <- cluster_spinglass(huri_g, spins = 2)

# OCG
ocg <- getOCG.clusters(as_data_frame(huri_g_sub))
ocg_max_Cliques <- getOCG.clusters(as_data_frame(huri_g_sub), init.class.sys = 1, cent.class.sys = 0)
ocg_edge_Cliques <- getOCG.clusters(as_data_frame(huri_g_sub), init.class.sys = 2, cent.class.sys = 0, min.class = 100)
ocg_edge_Cliques20 <- getOCG.clusters(as_data_frame(huri_g_sub), init.class.sys = 2, cent.class.sys = 0, min.class = 20)

ocg_min_cluster10 <- getOCG.clusters(as_data_frame(huri_g_sub), min.class = 10, cent.class.sys = 0) # not good, 1 huge community (# 6349)

ocg_min_cluster6 <- getOCG.clusters(as_data_frame(huri_g_sub), min.class = 6, cent.class.sys = 0) # not good, 1 huge community ()

huri_community <- list(
    walktrap = data.frame(gene = cw$names, membership = cw$membership),
    greedy = data.frame(gene = cfg$names, membership = cfg$membership),
    propagating = data.frame(gene = clp$names, membership = clp$membership),
    leading_eigenvector = data.frame(gene = cle$names, membership = cle$membership),
    multi_level = data.frame(gene = cl$names, membership = cl$membership))
write.xlsx(huri_community, file = "/tmp/HuRI_communities.xlsx")
