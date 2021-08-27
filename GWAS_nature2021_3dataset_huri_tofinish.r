# COVID-19 GWAS hit analysis
# Lin Chung-win
# 11.08.2021
# result saved as "~/INET/GWAS_3data_HuRI.RData"

######
# load package
library(openxlsx)
library(igraph)
library(rethinking)
library(gprofiler2)
library(linkcomm)
library(gplots)
# functions

subnetwork <- function(network, node) {
    gwas_hit_1st <- make_ego_graph(network, nodes = node, order = 1, mode = "all")
    gwas_all_list_df <- lapply(gwas_hit_1st, as_data_frame)
    gwas_all_df <- do.call(rbind, gwas_all_list_df)
    gwas_all_g_merge <- graph_from_data_frame(gwas_all_df, directed = FALSE)
    gwas_all_final <- simplify(induced_subgraph(network, names(V(gwas_all_g_merge))), remove.loops = F)
    return(gwas_all_final)
}

rewireMulti <- function(network, gwas, ctcl, hosp, infct, husci, gordon, stukalov, v_from, v_to) {
    df <- c()
    re <- simplify(rewire(network, keeping_degseq(niter = gsize(network) * 10)), remove.loops = FALSE)
    merged <- subnetwork(re, gwas) # get GWAS+1 subnetwork
    df <- c(df, table(V(merged)$name %in% husci)["TRUE"]) # get viral targets from HuSCI in GWAS+1 subnetwork
    df <- c(df, table(V(merged)$name %in% gordon)["TRUE"]) # get viral targets from Gordon in GWAS+1 subnetwork
    df <- c(df, table(V(merged)$name %in% stukalov)["TRUE"]) # get viral targets from Stukalov in GWAS+1 subnetwork
    df <- c(df, gsize(merged)) # get network size of GWAS+1 subnetwork
    # average shortest path of GWAS subnetwork from the rewired HuRI
    df <- c(df, mean(distances(merged, v = v_from, to = v_to, mode = "all")))

    merged <- subnetwork(re, ctcl) # get GWAS+1 subnetwork only from critical illness candidates
    df <- c(df, table(V(merged)$name %in% husci)["TRUE"])
    df <- c(df, table(V(merged)$name %in% gordon)["TRUE"])
    df <- c(df, table(V(merged)$name %in% stukalov)["TRUE"])
    df <- c(df, gsize(merged))
    df <- c(df, mean(distances(merged, v = v_from, to = v_to, mode = "all")))

    merged <- subnetwork(re, hosp) # get GWAS+1 subnetwork only from hospitalized candidates
    df <- c(df, table(V(merged)$name %in% husci)["TRUE"])
    df <- c(df, table(V(merged)$name %in% gordon)["TRUE"])
    df <- c(df, table(V(merged)$name %in% stukalov)["TRUE"])
    df <- c(df, gsize(merged))
    df <- c(df, mean(distances(merged, v = v_from, to = v_to, mode = "all")))

    merged <- subnetwork(re, infct) # get GWAS+1 subnetwork only from reported infection candidates
    df <- c(df, table(V(merged)$name %in% husci)["TRUE"])
    df <- c(df, table(V(merged)$name %in% gordon)["TRUE"])
    df <- c(df, table(V(merged)$name %in% stukalov)["TRUE"])
    df <- c(df, gsize(merged))
    df <- c(df, mean(distances(merged, v = v_from, to = v_to, mode = "all")))

    return(df)
}
plotHist <- function(value, title, length, y1, y2) {
    med <- median(value)
    dens_gwas <- hist(value, breaks = c(0:(max(value) + 1)), plot = FALSE, right = FALSE)
    plot(dens_gwas, col = rgb(0.75, 0.75, 0.75, 1/2), freq = FALSE, border = NA, las = 1, xlim = c(0, 20), xaxt = "n", xlab = "Number of viral targets", ylab = "Frequency", main = "", cex.sub = 0.5)
    mytitle <- paste0("COVID19 GWAS subnetwork\nviral targets in ", title)
    mtext(side = 3, line = 1, cex = 1, mytitle)
    mtext(side = 3, line = 0.2, cex = 0.8, "subnetwork extracted from HuRI")
    axis(side = 1, at = seq(0, 20, by = 5) + 0.5, labels = seq(0, 20, by = 5))
    arrows(length + 0.5, y1, length + 0.5, 0, col = "#922687", lwd = 2, length = 0.1)
    text(med + 2, max(dens_gwas$counts / 10000), paste0("median = ", med), col = "grey", cex = 0.5)
    text(length - 2, y2, paste0("observed = ", length, "\np = ", table(value >= length)["TRUE"] / 10000), cex = 0.4, pos = 4)
}

plotInteraction <- function(value, ymax, observe, phenotype) {
    dens_gwas <- hist(value, breaks = 20, plot = FALSE, right = FALSE)
    plot(dens_gwas, col = rgb(0.75, 0.75, 0.75, 1/2), border = NA, las = 1, yaxt = "n", xlab = "Number of interactions", main = "", cex.sub = 0.5)
    mtext(side = 3, line = 1, cex = 1, paste0("COVID19 GWAS subnetwork: ", phenotype))
    mtext(side = 3, line = 0.2, cex = 0.8, "subnetwork extracted from HuRI")
    axis(side = 2, at = seq(0, ymax, by = 500), labels = seq(0, ymax/10000, by = 0.05), las = 1)
    arrows(observe, 200, observe, 0, col = "#922687", lwd = 2, length = 0.1)
    text(median(value), max(dens_gwas$counts), paste0("median = ", median(value)), col = "grey", cex = 0.5)
    text(observe - (observe / 10), 350, paste0("observed = ", observe, "\np = ", table(value >= observe)["TRUE"]/10000), cex = 0.4, pos = 4)    
}

plotDistance <- function(value, ymax, observe, phenotype) {
    dens_gwas <- hist(value, breaks = 20, plot = FALSE, right = FALSE)
    plot(dens_gwas, col = rgb(0.75, 0.75, 0.75, 1/2), border = NA, las = 1, xlim = c(1, round(max(value), 0)), yaxt = "n", xlab = "Average shortest path", main = "", cex.sub = 0.5)
    mtext(side = 3, line = 1, cex = 1, paste0("COVID19 GWAS subnetwork: ", phenotype))
    mtext(side = 3, line = 0.2, cex = 0.8, "subnetwork extracted from HuRI")
    axis(side = 2, at = seq(0, ymax, by = 500), labels = seq(0, ymax/10000, by = 0.05), las = 1)
    arrows(observe, 300, observe, 0, col = "#922687", lwd = 2, length = 0.1)
    text(round(median(value), 2), max(dens_gwas$counts), paste0("median = ", round(median(value), 2)), col = "grey", cex = 0.5)
    text(round(observe, 1) - 0.3, 400, paste0("observed = ", round(observe, 2), "\np = ", table(value >= round(observe, 2))["TRUE"]/10000), cex = 0.4, pos = 4)        
}
######
# load dataset
huri <- read.csv("~/Documents/INET-work/references/HuRI_binaryPPI/HuRI_Tong_withSymbol.csv", header = T)
gwas <- read.xlsx("~/Documents/INET-work/virus_network/statistic_results/GWAS/COVID_GWAS hits_v2.xlsx") # ref: Mapping the human genetic architecture of COVID-19 (url: https://doi.org/10.1038/s41586-021-03767-x)
husci <- read.csv("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/binary_node_1126.csv", header = TRUE)

# add Gordon and Stukalov data
gordon <- read.xlsx("/Volumes/GoogleDrive/My Drive/VirHostome_CW/GitHub/data/extended_table/Extended_Table_2_PPIs.xlsx", sheet = "Gordon")
stukalov <- read.xlsx("/Volumes/GoogleDrive/My Drive/VirHostome_CW/GitHub/data/extended_table/Extended_Table_2_PPIs.xlsx", sheet = "Stukalov")

######
# 1. HuRI graph generation
huri_symbol <- huri[, c(5:6)]
huri_g_ori <- graph_from_data_frame(huri_symbol, directed = FALSE) # V:8274, E:52573
huri_g <- simplify(huri_g_ori, remove.loops = FALSE) # V:8274, E:52558

gwas_huri <- gwas$All.LD[gwas$All.LD %in% V(huri_g)$name] # GWAS hit in HuRI
gwas_huri_paralogs <- gwas_huri[c(1, 3, 5:9, 11:17)] # omit OAS2, ICAM3 and ICAM4

# HuSCI, Gordon et al and Stukalov at al in HuRI
husci_sym <- husci[husci$group == "human", "node"]
husci_huri <- V(huri_g)$name[V(huri_g)$name %in% husci_sym] 

gordon_sym <- unique(gordon$PreyGene)
gordon_huri <- V(huri_g)$name[V(huri_g)$name %in% gordon_sym]

stukalov_sym <- unique(stukalov$human)
stukalov_huri <- V(huri_g)$name[V(huri_g)$name %in% stukalov_sym]

# individual phenotype
ctcl <- gwas[, 2][gwas[, 5] == 1]
ctcl <- unique(ctcl[!is.na(ctcl)])
ctcl_huri <- ctcl[ctcl %in% V(huri_g)$name]
ctcl_huri_paralogs <- ctcl_huri[c(1, 3, 5:7, 9, 10)] # only OAS1, ICAM3

hosp <- gwas[, 2][gwas[, 6] == 1]
hosp <- unique(hosp[!is.na(hosp)])
hosp_huri <- hosp[hosp %in% V(huri_g)$name]
hosp_huri_paralogs <- hosp_huri[c(1, 3, 5:7, 9:13)] # only OAS1, ICAM3

infct <- gwas[, 2][gwas[, 7] == 1]
infct <- unique(infct[!is.na(infct)])
infct_huri <- infct[infct %in% V(huri_g)$name]
infct_huri_paralogs <- infct_huri[c(1:3, 5, 6)] # only OAS1

######
# 2. observation
# subnetwork establishment
observation_all <- subnetwork(huri_g, gwas_huri)
print("\n"); print("GWAS subnetwork")
observation_all
observation_paralogs <- subnetwork(huri_g, gwas_huri_paralogs)
print("\n"); print("GWAS subnetwork, without paralogs")
observation_paralogs

observation_ctcl <- subnetwork(huri_g, ctcl_huri)
print("\n"); print("Critical illness GWAS subnetwork")
observation_ctcl
observation_ctcl_paralogs <- subnetwork(huri_g, ctcl_huri_paralogs)
print("\n"); print("Critical illness GWAS subnetwork, without paralogs")
observation_ctcl_paralogs

observation_hosp <- subnetwork(huri_g, hosp_huri)
print("\n"); print("Hospitalization GWAS subnetwork")
observation_hosp
observation_hosp_paralogs <- subnetwork(huri_g, hosp_huri_paralogs)
print("\n"); print("Hospitalization GWAS subnetwork, without paralogs")
observation_hosp_paralogs

observation_infct <- subnetwork(huri_g, infct_huri)
print("\n"); print("Reported infection GWAS subnetwork")
observation_infct
observation_infct_paralogs <- subnetwork(huri_g, infct_huri_paralogs)
print("\n"); print("Reported infection GWAS subnetwork, without paralogs")
observation_infct_paralogs

# a. viral targets, all GWAS proteins
husci_viral_targets_all <- V(observation_all)$name[V(observation_all)$name %in% husci_sym]
print("\n"); print("Virtal targets from HuSCI in GWAS subnetwork")
husci_viral_targets_all
gordon_viral_targets_all <- V(observation_all)$name[V(observation_all)$name %in% gordon_sym]
print("\n"); print("Virtal targets from Gordon at al in GWAS subnetwork")
gordon_viral_targets_all
stukalov_viral_targets_all <- V(observation_all)$name[V(observation_all)$name %in% stukalov_sym]
print("\n"); print("Virtal targets from Stukalov et al in GWAS subnetwork")
stukalov_viral_targets_all

# viral targets, paralogs
husci_viral_targets_paralogs <- V(observation_paralogs)$name[V(observation_paralogs)$name %in% husci_sym]
print("\n"); print("Virtal targets from HuSCI in GWAS subnetwork, without paralogs")
husci_viral_targets_paralogs
gordon_viral_targets_paralogs <- V(observation_paralogs)$name[V(observation_paralogs)$name %in% gordon_sym]
print("\n"); print("Virtal targets from Gordon at al in GWAS subnetwork, without paralogs")
gordon_viral_targets_paralogs
stukalov_viral_targets_paralogs <- V(observation_paralogs)$name[V(observation_paralogs)$name %in% stukalov_sym]
print("\n"); print("Virtal targets from Gordon at al in GWAS subnetwork, without paralogs")
stukalov_viral_targets_paralogs

# viral targets, critical illness GWAS proteins
husci_viral_targets_ctcl <- V(observation_ctcl)$name[V(observation_ctcl)$name %in% husci_sym]
print("\n"); print("Virtal targets from HuSCI in critical illness GWAS subnetwork")
husci_viral_targets_ctcl
gordon_viral_targets_ctcl <- V(observation_ctcl)$name[V(observation_ctcl)$name %in% gordon_sym]
print("\n"); print("Virtal targets from Gordon at al in critical illness GWAS subnetwork")
gordon_viral_targets_ctcl
stukalov_viral_targets_ctcl <- V(observation_ctcl)$name[V(observation_ctcl)$name %in% stukalov_sym]
print("\n"); print("Virtal targets from Gordon at al in critical illness GWAS subnetwork")
stukalov_viral_targets_ctcl

# viral targets, critical illness paralogs
husci_viral_targets_ctcl_paralogs <- V(observation_ctcl_paralogs)$name[V(observation_ctcl_paralogs)$name %in% husci_sym]
print("\n"); print("Virtal targets from HuSCI in critical illness GWAS subnetwork, without paralogs")
husci_viral_targets_ctcl_paralogs
gordon_viral_targets_ctcl_paralogs <- V(observation_ctcl_paralogs)$name[V(observation_ctcl_paralogs)$name %in% gordon_sym]
print("\n"); print("Virtal targets from Gordon at al in critical illness GWAS subnetwork, without paralogs")
gordon_viral_targets_ctcl_paralogs
stukalov_viral_targets_ctcl_paralogs <- V(observation_ctcl_paralogs)$name[V(observation_ctcl_paralogs)$name %in% stukalov_sym]
print("\n"); print("Virtal targets from Gordon at al in critical illness GWAS subnetwork, without paralogs")
stukalov_viral_targets_ctcl_paralogs

# viral targets, hospitalization GWAS proteins
husci_viral_targets_hosp <- V(observation_hosp)$name[V(observation_hosp)$name %in% husci_sym]
print("\n"); print("Virtal targets from HuSCI in hospitalization GWAS subnetwork")
husci_viral_targets_hosp
gordon_viral_targets_hosp <- V(observation_hosp)$name[V(observation_hosp)$name %in% gordon_sym]
print("\n"); print("Virtal targets from Gordon at al in hospitalization GWAS subnetwork")
gordon_viral_targets_hosp
stukalov_viral_targets_hosp <- V(observation_hosp)$name[V(observation_hosp)$name %in% stukalov_sym]
print("\n"); print("Virtal targets from Gordon at al in hospitalization GWAS subnetwork")
stukalov_viral_targets_hosp

# viral targets, hospiptalization paralogs
husci_viral_targets_hosp_paralogs <- V(observation_hosp_paralogs)$name[V(observation_hosp_paralogs)$name %in% husci_sym]
print("\n"); print("Virtal targets from HuSCI in hospitalization GWAS subnetwork")
husci_viral_targets_hosp_paralogs
gordon_viral_targets_hosp_paralogs <- V(observation_hosp_paralogs)$name[V(observation_hosp_paralogs)$name %in% gordon_sym]
print("\n"); print("Virtal targets from Gordon at al in hospitalization GWAS subnetwork")
gordon_viral_targets_hosp_paralogs
stukalov_viral_targets_hosp_paralogs <- V(observation_hosp_paralogs)$name[V(observation_hosp_paralogs)$name %in% stukalov_sym]
print("\n"); print("Virtal targets from Gordon at al in hospitalization GWAS subnetwork")
stukalov_viral_targets_hosp_paralogs

# viral targets, reported infection GWAS proteins
husci_viral_targets_infct <- V(observation_infct)$name[V(observation_infct)$name %in% husci_sym]
print("\n"); print("Virtal targets from HuSCI in reported infection GWAS subnetwork")
husci_viral_targets_infct
gordon_viral_targets_infct <- V(observation_infct)$name[V(observation_infct)$name %in% gordon_sym]
print("\n"); print("Virtal targets from Gordon at al in reported infection GWAS subnetwork")
gordon_viral_targets_infct
stukalov_viral_targets_infct <- V(observation_infct)$name[V(observation_infct)$name %in% stukalov_sym]
print("\n"); print("Virtal targets from Gordon at al in reported infection GWAS subnetwork")
stukalov_viral_targets_infct

# viral targets, reported infection paralogs
husci_viral_targets_infct_paralogs <- V(observation_infct_paralogs)$name[V(observation_infct_paralogs)$name %in% husci_sym]
print("\n"); print("Virtal targets from HuSCI in reported infection GWAS subnetwork")
husci_viral_targets_infct_paralogs
gordon_viral_targets_infct_paralogs <- V(observation_infct_paralogs)$name[V(observation_infct_paralogs)$name %in% gordon_sym]
print("\n"); print("Virtal targets from Gordon at al in reported infection GWAS subnetwork")
gordon_viral_targets_infct_paralogs
stukalov_viral_targets_infct_paralogs <- V(observation_infct_paralogs)$name[V(observation_infct_paralogs)$name %in% stukalov_sym]
print("\n"); print("Virtal targets from Gordon at al in reported infection GWAS subnetwork")
stukalov_viral_targets_infct_paralogs

# b. interaction
interaction_all <- gsize(observation_all)
print("\n"); print("Interaction of GWAS subnetwork:" )
interaction_all

interaction_ctcl <- gsize(observation_ctcl)
print("\n"); print("Interaction of critical illness GWAS subnetwork:" )
interaction_ctcl

interaction_hosp <- gsize(observation_hosp)
print("\n"); print("Interaction of hospitalization GWAS subnetwork:" )
interaction_hosp

interaction_infct <- gsize(observation_infct)
print("\n"); print("Interaction of reported infection GWAS subnetwork:" )
interaction_infct

interaction_paralogs <- gsize(observation_paralogs)
print("\n"); print("Interaction of GWAS subnetwork, without paralogs:" )
interaction_paralogs

interaction_ctcl_paralogs <- gsize(observation_ctcl_paralogs)
print("\n"); print("Interaction of critical illness GWAS subnetwork, without paralogs:" )
interaction_ctcl_paralogs

interaction_hosp_paralogs <- gsize(observation_hosp_paralogs)
print("\n"); print("Interaction of hospitalization GWAS subnetwork, without paralogs:" )
interaction_hosp_paralogs

interaction_infct_paralogs <- gsize(observation_infct_paralogs)
print("\n"); print("Interaction of reported infection GWAS subnetwork, without paralogs:" )
interaction_infct_paralogs

# c. average shortest path between GWAS proteins
gwas_protein_shortest_path_all <- mean(distances(observation_all, v = gwas_huri, to = gwas_huri, mode = "all"))
gwas_protein_shortest_path_paralogs <- mean(distances(observation_all, v = gwas_huri_paralogs, to = gwas_huri_paralogs, mode = "all"))

gwas_protein_shortest_path_ctcl <- mean(distances(observation_ctcl, v = ctcl_huri, to = ctcl_huri, mode = "all"))
gwas_protein_shortest_path_ctcl_paralogs <- mean(distances(observation_ctcl_paralogs, v = ctcl_huri_paralogs, to = ctcl_huri_paralogs, mode = "all"))

gwas_protein_shortest_path_hosp <- mean(distances(observation_hosp, v = hosp_huri, to = hosp_huri, mode = "all"))
gwas_protein_shortest_path_hosp_paralogs <- mean(distances(observation_hosp, v = hosp_huri_paralogs, to = hosp_huri_paralogs, mode = "all"))

# 3.1. do statistical analysis for all 17 GWAS hit candidate genes # time consuming
# huri_list <- list()
# for (i in 1:10000) {
#     huri_list[[i]] <- huriRewire()
# }

# HuSCI, Gordon and Stukalov
all_re <- c()
all_re <- c(all_re, mcreplicate(10000, huriRewireMulti(gwas_huri, ctcl_huri, hosp_huri, infct_huri, husci_sym, gordon_sym, stukalov_sym), mc.cores = detectCores()))
all_re[is.na(all_re)] <- 0

all_re2 <- c()
all_re2 <- c(all_re2, mcreplicate(10000, huriRewireMulti(gwas_huri_paralogs, ctcl_huri2, hosp_huri2, infct_huri2, husci_sym, gordon_sym, stukalov_sym), mc.cores = detectCores()))
all_re2[is.na(all_re2)] <- 0

# choose 1 of 2 randomization results
{
    # randomization <- all_re
    randomization <- all_re2
}

all_re_df <- data.frame(matrix(randomization, ncol = 20, byrow = T))
names(all_re_df) <- c(
    "allGWAS_viral_target_inHuSCI",
    "allGWAS_viral_target_inGordon",
    "allGWAS_viral_target_inStukalov",
    "allGWAS_subnetworkSize",
    "allGWAS_shortestpath",
    "ctclGWAS_viral_target_inHuSCI",
    "ctclGWAS_viral_target_inGordon",
    "ctclGWAS_viral_target_inStukalov",
    "ctclGWAS_subnetworkSize",
    "ctclGWAS_shortestpath",
    "hospGWAS_viral_target_inHuSCI",
    "hospGWAS_viral_target_inGordon",
    "hospGWAS_viral_target_inStukalov",
    "hospGWAS_subnetworkSize",
    "hospGWAS_shortestpath",
    "infctGWAS_viral_target_inHuSCI",
    "infctGWAS_viral_target_inGordon",
    "infctGWAS_viral_target_inStukalov",
    "infctGWAS_subnetworkSize",
    "infctGWAS_shortestpath"
)
write.xlsx(all_re_df, file = "~/Documents/INET-work/virus_network/statistic_results/GWAS/Nature2021b_3dataset_paralog.xlsx", overwrite = T)
all_re_df_plot <- all_re_df[, c(1:3, 6:8, 11:13, 16:18)]
all_length <- c(gwas_all_husci_length,
    gwas_all_gordon_length,
    gwas_all_stukalov_length,
    gwas_ctcl_husci_length,
    gwas_ctcl_gordon_length,
    gwas_ctcl_stukalov_length,
    gwas_hosp_husci_length,
    gwas_hosp_gordon_length,
    gwas_hosp_stukalov_length,
    gwas_infct_husci_length,
    gwas_infct_gordon_length,
    gwas_infct_stukalov_length
)
title <- rep(c("HuSCI", "Gordon et al", "Stukalov et al"), 4)
phenotype <- rep(c(
    "all 17 genes and 1st interactors",
    "critical illness, 10 genes and 1st interactors",
    "hospitalization, 13 genes and 1st interactors",
    "reported infection, 6 genes and 1st interactors"
    ), each = 3)

phenotype2 <- rep(c(
    "all 14 genes",
    "critical illness, 7 genes",
    "hospitalization, 10 genes",
    "reported infection, 5 genes"
    ), each = 3)

y1 <- c(0.03, 0.03, 0.03, 0.03, 0.04, 0.03, 0.03, 0.03, 0.03, 0.04, 0.05, 0.05)
y2 <- c(0.05, 0.05, 0.05, 0.05, 0.06, 0.05, 0.05, 0.05, 0.05, 0.06, 0.07, 0.07)
# plotting
pdf("~/Documents/INET-work/virus_network/figure_results/GWAS/Nature2021b_3dataset_HuRI_paralog.pdf", width = 3, height = 3)
par(mgp = c(2, 0.7, 0), ps = 8)
for (i in 1:12) {
    plotHist(
        all_re_df_plot[, i],
        title[i],
        # phenotype2[i],
        all_length[i],
        y1[i], y2[i]
        )
}

# interaction, all GWAS
plotInteraction(all_re_df[, 4], 1000, gsize(gwas_all_final), "all GWAS")
# interaction, all Critical illness
plotInteraction(all_re_df[, 9], 2000, gsize(ctcl_1st), "critical")
# interaction, all hospitalization
plotInteraction(all_re_df[, 14], 1000, gsize(hosp_1st), "hospitalization")
# interaction, all reported infection
plotInteraction(all_re_df[, 19], 1000, gsize(infct_1st), "infection")

# average shortest path, all GWAS
plotDistance(all_re_df[, 5], 1200, gwas_all_mean_dist, "all GWAS")
# average shortest path, all critical illness
plotDistance(all_re_df[, 10], 2000, gwas_ctcl_mean_dist, "critical illness")
# average shortest path, all hospitalization
plotDistance(all_re_df[, 15], 2000, gwas_hosp_mean_dist, "hospitalization")
# average shortest path, all reported infection
plotDistance(all_re_df[, 20], 1500, gwas_infct_mean_dist, "reported infection")
dev.off()

boxplot(inHuSCI_summary, las = 1)
pdf("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/random_GWAS_viral_target_v2_left_closed.pdf", width = 4, height = 4)
for (i in seq(1, length(inHuSCI_summary), by = 2)) {
    dens_gwas <- hist(inHuSCI_summary[[i]], breaks = seq(0, max(inHuSCI_summary[[i]]), by = 1), plot = FALSE, right = TRUE)
    plot(dens_gwas, col = rgb(0.75, 0.75, 0.75, 1/2), freq = FALSE, border = NA, las = 1, xlim = c(0, max(inHuSCI_summary[[i]])), xaxt = "n", xlab = "Number of viral targets", ylab = "Frequency", main = "", cex.sub = 0.5)
    lines(dens_gwas$mids, dens_gwas$density, col = "darkgrey", lwd = 3)
    mytitle <- paste0("COVID19 GWAS loci candidate genes\n", names(inHuSCI_summary[i]))
    mtext(side = 3, line = 1, cex = 1, mytitle)
    axis(side = 1, at = seq(min(dens_gwas$mids), max(dens_gwas$mids) , 2), labels = seq(min(dens_gwas$mids) - 0.5, max(dens_gwas$mids) - 0.5, 2))
    arrows(inHuSCI_length[[i]] + 0.5, 0.04, inHuSCI_length[[i]] + 0.5, 0.01, col = "#922687", lwd = 2, length = 0.1)
    text(inHuSCI_length[[i]], 0.05, paste0("observed = ", inHuSCI_length[[i]], "\n p = ", table(inHuSCI_summary[[i]] >= inHuSCI_length[[i]])["TRUE"]/10000), cex = 0.4, pos = 4)
}
dev.off()

pdf("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/random_GWAS_viral_target_v2_right_closed.pdf", width = 4, height = 4)
for (i in seq(1, length(inHuSCI_summary), by = 2)) {
    dens_gwas <- hist(inHuSCI_summary[[i]], breaks = seq(0, max(inHuSCI_summary[[i]]), by = 1), plot = FALSE, right = FALSE)
    plot(dens_gwas, col = rgb(0.75, 0.75, 0.75, 1/2), freq = FALSE, border = NA, las = 1, xlim = c(0, max(inHuSCI_summary[[i]])), xaxt = "n", xlab = "Number of viral targets", ylab = "Frequency", main = "", cex.sub = 0.5)
    lines(dens_gwas$mids, dens_gwas$density, col = "darkgrey", lwd = 3)
    mytitle <- paste0("COVID19 GWAS loci candidate genes\n", names(inHuSCI_summary[i]))
    mtext(side = 3, line = 1, cex = 1, mytitle)
    axis(side = 1, at = seq(min(dens_gwas$mids), max(dens_gwas$mids) , 2), labels = seq(min(dens_gwas$mids) - 0.5, max(dens_gwas$mids) - 0.5, 2))
    arrows(inHuSCI_length[[i]] + 0.5, 0.04, inHuSCI_length[[i]] + 0.5, 0.01, col = "#922687", lwd = 2, length = 0.1)
    text(inHuSCI_length[[i]], 0.05, paste0("observed = ", inHuSCI_length[[i]], "\n p = ", table(inHuSCI_summary[[i]] >= inHuSCI_length[[i]])["TRUE"]/10000), cex = 0.4, pos = 4)
}
dev.off()

#######
# 4. keep one ortholog at a time
ortholog_count <- list()
for (i in 1:dim(all_candidate)[1]) {
    ortholog_gra <- combineNetwork(huri_g, c(as.matrix(all_candidate[i, ])))
    V(ortholog_gra)$color <- ifelse(V(ortholog_gra)$name %in% husci_sym, "red", "grey")
    name <- paste0(c(as.matrix(all_candidate[i, ])), collapse = "-")
    ortholog_count[[name]] <- data.frame(Vcount = vcount(ortholog_gra), Ecount = ecount(ortholog_gra))
    pdf(paste0("/tmp/", name, ".pdf"))
    plot(ortholog_gra, vertex.size = 4, vertex.label.cex = .5, vertex.label.dist = 1, vertex.label.color = "black")
    title(name, cex.main = .5)
    dev.off()
}

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
# 5. community detected
gwas_all_strong <- components(gwas_all_final, mode = "strong")
gwas_all_strong1 <- induced_subgraph(gwas_all_final, names(gwas_all_strong$membership[gwas_all_strong$membership == 1]))
gwas_all_strong1 <- simplify(gwas_all_strong1) # looks like doesn't hurt
gwas_all_ocg <- getOCG.clusters(as_data_frame(gwas_all_strong1), init.class.sys = 3, cent.class.sys = 0)
## visualized with whole network
pdf('/tmp/gwas_ocg.pdf')
plotOCGraph(gwas_all_ocg)
dev.off()

gwas_all_ocg_dif <- sort(gwas_all_ocg$clustsizes, decreasing = T)[1] / sort(gwas_all_ocg$clustsizes, decreasing = T)[2]
## randomized network, gwas all
all_ocg_rand_count <- list()
for (i in 1:1350) {
    tryCatch({
        all_ocg_rand_count[[i]] <- ocgRewireRatio(gwas_huri)
    }, error = function(e){
        print("Wired error!")
    })
}
all_ocg_sum <- as.numeric(unlist(lapply(all_ocg_rand_count, function(x) x[1])))[c(1:1000)]
all_ocg_sum3 <- as.numeric(unlist(lapply(all_ocg_rand_count, function(x) x[3])))[c(1:1000)]
## visualized with boxplot
boxplot(all_ocg_sum, horizontal = TRUE, cex.axis = 2)
points(gwas_all_ocg_dif, 1, pch = 19, col = "red", cex = 1.5)
## visualized with histogram
pdf("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/random_GWAS_commuRatio_v2.pdf", width = 3, height = 3)
dens <- hist(all_ocg_sum, breaks = 15, plot = FALSE)
par(mgp = c(2, 0.7, 0), ps = 8)
plot(dens, col = rgb(0.75, 0.75, 0.75, 1/2), freq = FALSE, border = NA, las = 1, xlim = c(0, 20), xlab = "Ratio of communities", ylab = "Frequency", main = "", cex.sub = 0.5)
arrows(gwas_all_ocg_dif + 0.5, 0.1, gwas_all_ocg_dif + 0.5, 0.02, col = "#922687", lwd = 2, length = 0.1)
text(gwas_all_ocg_dif - 2, 0.12, paste0("observed = ", gwas_all_ocg_dif, "\np = ", table(all_ocg_sum >= gwas_all_ocg_dif)["TRUE"]/1000, "/", table(all_ocg_sum >= gwas_all_ocg_dif)["FALSE"]/1000), cex = 0.4, pos = 4)
dev.off()
## randomized network to reveal viral target significance


# 5.1.  based on 3 different phenotype
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

######
# save result
save.image("~/Documents/INET-work/virus_network/statistic_results/GWAS/Nature2021b_3data_HuRI.RData")
