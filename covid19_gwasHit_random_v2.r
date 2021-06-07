# COVID-19 GWAS hit analysis
# Lin Chung-win
## v2: this is very wired! One typo from Pascal makes me unable to replicate the analysis?!

######
# load package
library(openxlsx)
library(igraph)
library(rethinking)
library(gprofiler2)
library(linkcomm)
source("~/Documents/INET-work/virus_network/src/plotOCGGraph.r")

go2check <- function(x, organism) {
    goquery <- gost(query = x, organism = organism, correction_method = "bonferroni", evcodes = T)
    goquery$result$inCommunity <- paste0( goquery$meta$query_metadata$queries$query_1, collapse = ",")
    goquery$result$annotatedInCommunity <- paste0(goquery$meta$genes_metadata$query$query_1$ensgs, collapse = ",")
    return(goquery$result)
}

######
# load dataset
huri <- read.csv("~/Documents/INET-work/references/HuRI_binaryPPI/HuRI_Tong_withSymbol.csv", header = T)
gwas <- read.xlsx("~/Documents/INET-work/virus_network/statistic_results/GWAS/COVID_GWAS hits_v2.xlsx")
husci <- read.csv("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/binary_node_1126.csv", header = TRUE)
gwas_2 <- read.xlsx("~/Documents/INET-work/virus_network/statistic_results/GWAS/COVID_GWAS_hit_inHUSCI_v2.xlsx") # to filter hits with high degree
    # the gwas_2 was generated after 1st analysis
gwas_file <- "~/Documents/INET-work/virus_network/statistic_results/GWAS/COVID_GWAS hits_source_v2.xlsx"
gwas_candidate <- read.xlsx(gwas_file, sheet = "all")
gwas_ctcl <- read.xlsx(gwas_file, sheet = "critical")
gwas_hosp <- read.xlsx(gwas_file, sheet = "hospitalization")
gwas_infct <- read.xlsx(gwas_file, sheet = "infection")
######
# HuRI graph generation
huri_symbol <- huri[, c(5:6)]
huri_g <- graph_from_data_frame(huri_symbol, directed = FALSE) # V:8274, E:52573
huri_g <- simplify(huri_g, remove.loops = FALSE) # V:8274, E:52558
######
# protein list filter
husci_sym <- husci[husci$group == "human", "node"]
husci_huri <- V(huri_g)$name[V(huri_g)$name %in% husci_sym] # HuSCI in HuRI whole
gwas_huri <- gwas$All.LD[gwas$All.LD %in% V(huri_g)$name] # GWAS hit in HuRI
######
# functions
check <- function(x) {
    value <- table(V(x)$name %in% husci_sym)[2]
    value <- ifelse(is.na(value), "0", value)
    return(as.numeric(value))
}
checkUpdate <- function(x) {
    V(x)$name[V(x)$name %in% husci_sym]
}
combineNetwork <- function(network, node) {
    gwas_random_g <- make_ego_graph(network, nodes = node, order = 1, mode = "all")
    gwas_random_list_df <- lapply(gwas_random_g, as_data_frame)
    gwas_random_df <- do.call(rbind, gwas_random_list_df)
    gwas_random_g_merge <- graph_from_data_frame(gwas_random_df, directed = FALSE)
    gwas_random_final <- simplify(induced_subgraph(huri_g, names(V(gwas_random_g_merge))), remove.loops = FALSE)
    return(gwas_random_final)
}
######
# interactor of GWAS hit
gwas_hit_1st <- make_ego_graph(huri_g, nodes = sort(gwas$All.LD[gwas$All.LD %in% V(huri_g)$name]), order = 1, mode = "all") #17 of 42 in HuRI
gwas_check <- length(unique(unlist(lapply
    (gwas_hit_1st, checkUpdate)))) #11 interactors
gwas_all_candidate <- expand.grid(strsplit(gwas_candidate[, 5], split = ","))
gwas_all_candidate_uniq <- unname(unique(as.matrix(unlist(gwas_all_candidate)))[, 1])
gwas_all_candidate_huri <- gwas_all_candidate_uniq[gwas_all_candidate_uniq %in% V(huri_g)$name]
gwas_ctcl_candidate <- expand.grid(strsplit(gwas_ctcl[, 5], split = ","))
gwas_ctcl_candidate_uniq <- unname(unique(as.matrix(unlist(gwas_ctcl_candidate)))[, 1])
gwas_ctcl_candidate_huri <- gwas_ctcl_candidate_uniq[gwas_ctcl_candidate_uniq %in% V(huri_g)$name]
gwas_hosp_candidate <- expand.grid(strsplit(gwas_hosp[, 5], split = ","))
gwas_hosp_candidate_uniq <- unname(unique(as.matrix(unlist(gwas_hosp_candidate)))[, 1])
gwas_hosp_candidate_huri <- gwas_hosp_candidate_uniq[gwas_hosp_candidate_uniq %in% V(huri_g)$name]
gwas_infct_candidate <- expand.grid(strsplit(gwas_infct[, 5], split = ","))
gwas_infct_candidate_uniq <- unname(unique(as.matrix(unlist(gwas_infct_candidate)))[, 1])
gwas_infct_candidate_huri <- gwas_infct_candidate_uniq[gwas_infct_candidate_uniq %in% V(huri_g)$name]
######
# GO enrichment analysis
cluster_list <- lapply(gwas_hit_1st, function(x) {
    df <- as_data_frame(x)
    unique(c(df$from, df$to))
})
names(cluster_list) <- sort(gwas$All.LD[gwas$All.LD %in% V(huri_g)$name])

bp_cust <- upload_GMT_file(gmtfile = "~/workplace/database/hsapien_HuRI_GOBP_EXP.gmt")
mf_cust <- upload_GMT_file(gmtfile = "~/workplace/database/hsapien_HuRI_GOMF_EXP.gmt")
cc_cust <- upload_GMT_file(gmtfile = "~/workplace/database/hsapien_HuRI_GOCC_EXP.gmt")
## individual subnetwork
gobp_check <- lapply(cluster_list, function(x) {go2check(x, bp_cust)})
gomf_check <- lapply(cluster_list, function(x) {go2check(x, mf_cust)})
gocc_check <- lapply(cluster_list, function(x) {go2check(x, cc_cust)})

bp <- do.call(rbind.data.frame, gobp_check[c(1:12, 14:16)])
mf <- do.call(rbind.data.frame, gomf_check[c(2, 4:12, 14, 16)])
cc <- do.call(rbind.data.frame, gocc_check[c(3, 5:7, 9, 12, 14)])
write.csv(rbind(bp, mf, cc)[, c(1:13, 15:18)], file = "/tmp/GWAS_GO_individual.csv")
## whole GWAS hits set
whole_gwas_bp <- go2check(unique(as.character(unlist(lapply(gwas_hit_1st, as_data_frame)))), bp_cust)
whole_gwas_mf <- go2check(unique(as.character(unlist(lapply(gwas_hit_1st, as_data_frame)))), mf_cust)
whole_gwas_cc <- go2check(unique(as.character(unlist(lapply(gwas_hit_1st, as_data_frame)))), cc_cust)
write.csv(rbind(whole_gwas_mf)[, c(1:13, 15:18)], file = "~/Documents/INET-work/virus_network/statistic_results/GWAS/whole_GWAS_GO_v2.csv")
######
# rewiring analysis of HuRI, to see if the HuSCI viral target is significant.
# load gwas loci info, with 3 phenotype (critical illness, hospitalization and infection)
gwas_all <- unique(gwas$All.LD)
gwas_all <- gwas_all[gwas_all != ""]

######
# subnetwork of GWAS hit from HuRI
## inherite from above code
gwas_all_list_df <- lapply(gwas_hit_1st, as_data_frame)
gwas_all_df <- do.call(rbind, gwas_all_list_df)
gwas_all_g_merge <- graph_from_data_frame(gwas_all_df, directed = FALSE)
## to have interaction between 1st interactors
gwas_all_final <- simplify(induced_subgraph(huri_g, names(V(gwas_all_g_merge))), remove.loops = FALSE) # V:118, E:370
# subgraph for critical illness, hospitalization and infection
gwas_ctcl_graph <- make_ego_graph(huri_g, nodes = gwas_ctcl_candidate_huri, mode = "all")
gwas_hosp_graph <- make_ego_graph(huri_g, nodes = gwas_hosp_candidate_huri, mode = "all")
gwas_infct_graph <- make_ego_graph(huri_g, nodes = gwas_infct_candidate_huri, mode = "all")
### HuSCI in 1st interactor
gwas_1stN_inHuSCI <- table(names(V(gwas_all_final)) %in% husci_sym)['TRUE']

gwas_ctcl_inHuSCI <- table( unique( unlist( lapply(gwas_ctcl_graph, function(x) names(V(x))) ) ) %in% husci_sym)["TRUE"]
gwas_hosp_inHuSCI <- table( unique( unlist( lapply(gwas_hosp_graph, function(x) names(V(x))) ) ) %in% husci_sym)["TRUE"]
gwas_infct_inHuSCI <- table( unique( unlist( lapply(gwas_infct_graph, function(x) names(V(x))) ) ) %in% husci_sym)["TRUE"]
######
# degree preserving randomized HuRI network and count final node list appeared in HuSCI
randomGwas <- function(network, node) {
    gwas_random_final <- combineNetwork(network, node)
    gwas_1stN_inHuSCI <- table(names(V(gwas_random_final)) %in% husci_sym)['TRUE']
    return(gwas_1stN_inHuSCI)
}
random_gwas_all <- c()
random_gwas_all <- c(random_gwas_all, mcreplicate(10000, randomGwas(huri_g, gwas_all_candidate_huri)), mc.cores = detectCores())
random_gwas_all_final <- as.numeric(random_gwas_all)
random_gwas_all_final[is.na(random_gwas_all_final)] <- 0

random_gwas_ctcl <- c()
random_gwas_ctcl <- c(random_gwas_ctcl, mcreplicate(10000, randomGwas(huri_g, ctcl), mc.cores = detectCores()))
random_gwas_ctcl_final <- as.numeric(random_gwas_ctcl)
random_gwas_ctcl_final[is.na(random_gwas_ctcl_final)] <- 0

random_gwas_hosp <- c()
random_gwas_hosp <- c(random_gwas_hosp, mcreplicate(10000, randomGwas(huri_g, hosp), mc.cores = detectCores()))
random_gwas_hosp_final <- as.numeric(random_gwas_hosp)
random_gwas_hosp_final[is.na(random_gwas_hosp_final)] <- 0

random_gwas_infct <- c()
random_gwas_infct <- c(random_gwas_infct, mcreplicate(10000, randomGwas(huri_g, infct), mc.cores = detectCores()))
random_gwas_infct_final <- as.numeric(random_gwas_infct)
random_gwas_infct_final[is.na(random_gwas_infct_final)] <- 0

######
# plot GWAS hits of 3 phenotypes permutation
source("gwas_phenotype_plot.r")

######
# 3 phenotypes GO enrichment analysis
ctcl_gobp <- go2check(unique(as.character(unlist(lapply(gwas_ctcl_graph, as_data_frame)))), bp_cust)
ctcl_gomf <- go2check(unique(as.character(unlist(lapply(gwas_ctcl_graph, as_data_frame)))), mf_cust)
ctcl_gocc <- go2check(unique(as.character(unlist(lapply(gwas_ctcl_graph, as_data_frame)))), cc_cust)
write.csv(rbind(ctcl_gomf)[, c(1:13, 15:18)], file = "/tmp/ctcl_GWAS_go.csv")

hosp_gobp <- go2check(unique(as.character(unlist(lapply(gwas_hosp_graph, as_data_frame)))), bp_cust)
hosp_gomf <- go2check(unique(as.character(unlist(lapply(gwas_hosp_graph, as_data_frame)))), mf_cust)
hosp_gocc <- go2check(unique(as.character(unlist(lapply(gwas_hosp_graph, as_data_frame)))), cc_cust)
write.csv(rbind(hosp_gomf)[, c(1:13, 15:18)], file = "/tmp/hosp_GWAS_go.csv")

infct_gobp <- go2check(unique(as.character(unlist(lapply(gwas_infct_graph, as_data_frame)))), bp_cust)
infct_gomf <- go2check(unique(as.character(unlist(lapply(gwas_infct_graph, as_data_frame)))), mf_cust)
infct_gocc <- go2check(unique(as.character(unlist(lapply(gwas_infct_graph, as_data_frame)))), cc_cust)
write.csv(rbind(infct_gomf, infct_gocc)[, c(1:13, 15:18)], file = "/tmp/infct_GWAS_go.csv")
######
# community detected based on 3 different phenotype
ctcl_network <- combineNetwork(huri_g, gwas_ctcl_candidate_huri)
ctcl_strong <- components(ctcl_network, mode = "strong")
ctcl_network_sub <- induced_subgraph(ctcl_network, names(ctcl_strong$membership[ctcl_strong[1]$membership == 1]))
hosp_network <- combineNetwork(huri_g, hosp)
hosp_strong <- components(hosp_network, mode = "strong")
hosp_network_sub <- induced_subgraph(hosp_network, names(hosp_strong$membership[hosp_strong[1]$membership == 1]))
infct_network <- combineNetwork(huri_g, infct)
infct_strong <- components(infct_network, mode = "strong")
infct_network_sub <- induced_subgraph(infct_network, names(infct_strong$membership[infct_strong[1]$membership == 1]))

ctcl_ocg <- getOCG.clusters(as_data_frame(ctcl_network_sub), init.class.sys = 1, cent.class.sys = 0)
hosp_ocg <- getOCG.clusters(as_data_frame(hosp_network_sub), init.class.sys = 1, cent.class.sys = 0)
infct_ocg <- getOCG.clusters(as_data_frame(infct_network_sub), init.class.sys = 1, cent.class.sys = 0)
######
# plotting graph
plotOCGraph()
plotOCGraph(ctcl_ocg)
plotOCGraph(hosp_ocg)
plotOCGraph(infct_ocg)
plot(gwas_all_final, vertex.size = 3, vertex.label = NA)
plot(ctcl_network, vertex.size = 3, vertex.label = NA)
plot(hosp_network, vertex.size = 3, vertex.label = NA)
plot(infct_network, vertex.size = 3, vertex.label = NA)
