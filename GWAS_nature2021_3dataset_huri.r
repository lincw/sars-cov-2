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
source("~/Documents/INET-work/virus_network/src/plotOCGGraph.r")

check <- function(x) {
    value <- table(V(x)$name %in% husci_sym)[2]
    value <- ifelse(is.na(value), "0", value)
    return(as.numeric(value))
}

checkUpdate <- function(x) {
    V(x)$name[V(x)$name %in% husci_sym]
}

source("~/Documents/INET-work/virus_network/src/combineNetwork.r")

huriRewire <- function(remove.loops = FALSE, ...) {
    huri_re <- rewire(huri_g, keeping_degseq(niter = gsize(huri_g) * 10))
    huri_sim <- simplify(huri_re, remove.loops = remove.loops)
    return(huri_sim)
}

huriRewireMulti <- function(gwas, ctcl, hosp, infct, husci, gordon, stukalov, remove.loops = FALSE) {
    df <- c()
    re <- huriRewire(remove.loops) # rewire huri
    merged <- combineNetwork(re, gwas) # get GWAS+1 subnetwork
    df <- c(df, table(V(merged)$name %in% husci)["TRUE"]) # get viral targets from HuSCI in GWAS+1 subnetwork
    df <- c(df, table(V(merged)$name %in% gordon)["TRUE"]) # get viral targets from Gordon in GWAS+1 subnetwork
    df <- c(df, table(V(merged)$name %in% stukalov)["TRUE"]) # get viral targets from Stukalov in GWAS+1 subnetwork
    df <- c(df, gsize(merged)) # get network size of GWAS+1 subnetwork
    # average shortest path of GWAS subnetwork from the rewired HuRI
    df <- c(df, mean_distance(merged))

    merged <- combineNetwork(re, ctcl) # get GWAS+1 subnetwork only from critical illness candidates
    df <- c(df, table(V(merged)$name %in% husci)["TRUE"])
    df <- c(df, table(V(merged)$name %in% gordon)["TRUE"])
    df <- c(df, table(V(merged)$name %in% stukalov)["TRUE"])
    df <- c(df, gsize(merged))
    df <- c(df, mean_distance(merged))

    merged <- combineNetwork(re, hosp) # get GWAS+1 subnetwork only from hospitalized candidates
    df <- c(df, table(V(merged)$name %in% husci)["TRUE"])
    df <- c(df, table(V(merged)$name %in% gordon)["TRUE"])
    df <- c(df, table(V(merged)$name %in% stukalov)["TRUE"])
    df <- c(df, gsize(merged))
    df <- c(df, mean_distance(merged))

    merged <- combineNetwork(re, infct) # get GWAS+1 subnetwork only from reported infection candidates
    df <- c(df, table(V(merged)$name %in% husci)["TRUE"])
    df <- c(df, table(V(merged)$name %in% gordon)["TRUE"])
    df <- c(df, table(V(merged)$name %in% stukalov)["TRUE"])
    df <- c(df, gsize(merged))
    df <- c(df, mean_distance(merged))

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
# edit plot parameters
trace("plot.igraph", edit = T)
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
huri_g_noloop <- simplify(huri_g_ori) # V:8274 E:52078 without self-loops
# protein list filter
husci_sym <- husci[husci$group == "human", "node"]
husci_huri <- V(huri_g)$name[V(huri_g)$name %in% husci_sym] # HuSCI in HuRI whole
gwas_huri <- gwas$All.LD[gwas$All.LD %in% V(huri_g)$name] # GWAS hit in HuRI

gwas_huri2 <- gwas_huri[c(1, 2, 5:9, 11:17)] # omit OAS2, ICAM3 and ICAM4

# Gordon and Stukalov in HuRI
gordon_sym <- unique(gordon$PreyGene)
gordon_huri <- V(huri_g)$name[V(huri_g)$name %in% gordon_sym]

stukalov_sym <- unique(stukalov$human)
stukalov_huri <- V(huri_g)$name[V(huri_g)$name %in% stukalov_sym]

######
# 2. interactor of GWAS hit
# gwas_hit_1st <- make_ego_graph(huri_g, nodes = sort(gwas$All.LD[gwas$All.LD %in% V(huri_g)$name]), order = 1, mode = "all") #17 of 42 in HuRI
gwas_hit_1st <- make_ego_graph(huri_g, nodes = sort(gwas_huri2), order = 1, mode = "all") #14 of 42 in HuRI

# 2.1 interactors of GWAS hits with critical illness phenotypes
ctcl <- gwas[, 2][gwas[, 5] == 1]
ctcl <- unique(ctcl[!is.na(ctcl)])
ctcl_huri <- ctcl[ctcl %in% V(huri_g)$name]

ctcl_huri2 <- ctcl_huri[c(1, 2, 5:7, 9, 10)] # only OAS1, ICAM1

ctcl_1st <- combineNetwork(huri_g, ctcl_huri2)
gwas_ctcl_husci <- V(ctcl_1st)$name[V(ctcl_1st)$name %in% husci_sym]
gwas_ctcl_husci_length <- length(gwas_ctcl_husci)
# average shortest path
gwas_ctcl_mean_dist <- mean_distance(ctcl_1st)

gwas_ctcl_gordon <- V(ctcl_1st)$name[V(ctcl_1st)$name %in% gordon_sym]
gwas_ctcl_gordon_length <- length(gwas_ctcl_gordon)
gwas_ctcl_stukalov <- V(ctcl_1st)$name[V(ctcl_1st)$name %in% stukalov_sym]
gwas_ctcl_stukalov_length <- length(gwas_ctcl_stukalov)

hosp <- gwas[, 2][gwas[, 6] == 1]
hosp <- unique(hosp[!is.na(hosp)])
hosp_huri <- hosp[hosp %in% V(huri_g)$name]

hosp_huri2 <- hosp_huri[c(1, 2, 5:7, 9:13)] # only OAS1, ICAM1

hosp_1st <- combineNetwork(huri_g, hosp_huri2)
gwas_hosp_husci <- V(hosp_1st)$name[V(hosp_1st)$name %in% husci_sym]
gwas_hosp_husci_length <- length(gwas_hosp_husci)
# average shortest path
gwas_hosp_mean_dist <- mean_distance(hosp_1st)

gwas_hosp_gordon <- V(hosp_1st)$name[V(hosp_1st)$name %in% gordon_sym]
gwas_hosp_gordon_length <- length(gwas_hosp_gordon)
gwas_hosp_stukalov <- V(hosp_1st)$name[V(hosp_1st)$name %in% stukalov_sym]
gwas_hosp_stukalov_length <- length(gwas_hosp_stukalov)

infct <- gwas[, 2][gwas[, 7] == 1]
infct <- unique(infct[!is.na(infct)])
infct_huri <- infct[infct %in% V(huri_g)$name]

infct_huri2 <- infct_huri[c(1:3, 5, 6)] # only OAS1

infct_1st <- combineNetwork(huri_g, infct_huri2)
gwas_infct_husci <- V(infct_1st)$name[V(infct_1st)$name %in% husci_sym]
gwas_infct_husci_length <- length(gwas_infct_husci)
# average shortest path
gwas_infct_mean_dist <- mean_distance(infct_1st)

gwas_infct_gordon <- V(infct_1st)$name[V(infct_1st)$name %in% gordon_sym]
gwas_infct_gordon_length <- length(gwas_infct_gordon)
gwas_infct_stukalov <- V(infct_1st)$name[V(infct_1st)$name %in% stukalov_sym]
gwas_infct_stukalov_length <- length(gwas_infct_stukalov)

overlap_gwas <- list(
    ctcl = ctcl,
    hosp = hosp,
    infct = infct
)
venn(overlap_gwas); title("Overlap between all LD")
overlap_huri <- list(
    ctcl = ctcl_huri,
    hosp = hosp_huri,
    infct = infct_huri
)
venn(overlap_huri); title("Overlap between all LD in HuRI")
overlap_husci <- list(
    ctcl = V(ctcl_1st)$name[V(ctcl_1st)$name %in% husci_sym],
    hosp = V(hosp_1st)$name[V(hosp_1st)$name %in% husci_sym],
    infct = V(infct_1st)$name[V(infct_1st)$name %in% husci_sym]
)
venn(overlap_husci); title("Overlap between 1st node in HuSCI")
######
# 3. **rewiring analysis of HuRI**, to see if the HuSCI viral target is significant.
# load gwas loci info, with 3 phenotype (critical illness, hospitalization and infection)
# subnetwork of GWAS hit from HuRI
# inherite from above code
gwas_all_list_df <- lapply(gwas_hit_1st, as_data_frame)
gwas_all_df <- do.call(rbind, gwas_all_list_df)
gwas_all_g_merge <- graph_from_data_frame(gwas_all_df, directed = FALSE)
# to have interaction between 1st interactors
gwas_all_final <- simplify(induced_subgraph(huri_g, names(V(gwas_all_g_merge))), remove.loops = F) 
gwas_all_mean_dist <- mean_distance(gwas_all_final)

gwas_all_husci <- V(gwas_all_final)$name[V(gwas_all_final)$name %in% husci_sym] # 11 viral targets
gwas_all_husci_length <- length(gwas_all_husci)

gwas_all_gordon <- V(gwas_all_final)$name[V(gwas_all_final)$name %in% gordon_sym] # 4 viral targets
gwas_all_gordon_length <- length(gwas_all_gordon)
gwas_all_stukalov <- V(gwas_all_final)$name[V(gwas_all_final)$name %in% stukalov_sym] # 9 viral targets
gwas_all_stukalov_length <- length(gwas_all_stukalov)

# visualization of gwas subnetwork
bd <- ifelse(V(gwas_all_final)$name %in% husci_sym, "orange", NA)
label <- ifelse(V(gwas_all_final)$name %in% c(husci_sym, gwas_huri), V(gwas_all_final)$name, NA)
color <- ifelse(V(gwas_all_final)$name %in% gwas_huri, "red", "grey")
pdf("/tmp/GWAS_subnetwork.pdf")
plot(gwas_all_final, vertex.frame.color = bd, vertex.size = 3, vertex.label.dist = 1, vertex.label.color = "black", vertex.label = label, vertex.color = color, vertex.frame.width = 2)
write_graph(gwas_all_final, "/tmp/GWAS_subnetwork.gml", format = "gml")
dev.off()
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
all_re2 <- c(all_re2, mcreplicate(10000, huriRewireMulti(gwas_huri2, ctcl_huri2, hosp_huri2, infct_huri2, husci_sym, gordon_sym, stukalov_sym), mc.cores = detectCores()))
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
