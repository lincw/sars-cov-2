# COVID-19 GWAS locus analysis with HuRI degree preserving randomized network
# Lin Chung-wen
######
# load package
library(openxlsx)
library(linkcomm)
library(rethinking)
library(plotrix)
source("~/Documents/INET-work/virus_network/src/plotOCGGraph.r")
# trace("plot.igraph",edit=TRUE)
# as....
# 179             cp <- matrix(c(x0, y0, x0 + 0.08, y0 + 0.08, x0 + 0.08,
# 180                 y0 - 0.08, x0, y0), ncol = 2, byrow = TRUE)
# end...
######
# functions
# 1. degree preserving randomized GWAS hits network from HuRI
randomHuri <- function(...) {
    huri_random <- simplify(rewire(huri_g, with = keeping_degseq(niter = gsize(huri_g) * 10, loops = FALSE)), remove.loops = FALSE)
    return(huri_random)
}
randomGwas <- function(node) {
    huri_random <- randomHuri()
    gwas_random_g <- make_ego_graph(huri_random, nodes = node, order = 1, mode = "all")
    gwas_random_list_df <- lapply(gwas_random_g, as_data_frame)
    gwas_random_df <- do.call(rbind, gwas_random_list_df)
    gwas_random_g_merge <- graph_from_data_frame(gwas_random_df, directed = FALSE)
    gwas_random_final <- simplify(induced_subgraph(huri_g, names(V(gwas_random_g_merge))), remove.loops = FALSE) # V: 209, E: 1058
    return(gwas_random_final)
}
# 2. average shortest path with permutated HuRI
avgShortInHuRI <- function(node) {
    huri_random <- randomHuri()
    ct <- mean(distances(huri_random, v = node, to = node, mode = "all"))
    return(ct)
}
# 3. degree preserving randomized HuRI network
source("~/Documents/INET-work/virus_network/src/countRandomize.r")
######
# load dataset
huri <- read.csv("~/Documents/INET-work/references/HuRI_binaryPPI/HuRI_Tong_withSymbol.csv", header = T)
husci <- read.csv("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/binary_node_1126.csv", header = TRUE)
gwas <- read.csv("~/Documents/INET-work/virus_network/statistic_results/GWAS/HuSCI_and_HuRI_+GWAS_v02_pfb.csv", header = T)
gwas_loci <- read.csv("~/Documents/INET-work/virus_network/statistic_results/GWAS/COVID_GWAS+1 node_pfb.csv", header = T)
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
gwas_all <- unique(gwas_loci$All.LD)
gwas_all <- gwas_all[gwas_all != ""]
######
# GWAS graph generation
gwas <- gwas[, c(9, 10)]
gwas_graph <- graph_from_data_frame(gwas, directed = FALSE)
gwas_graph <- simplify(gwas_graph, remove.loops = FALSE)
gwas_strong <- components(gwas_graph, mode = "strong")
gwas_graph_sub <- induced_subgraph(gwas_graph, names(gwas_strong$membership[gwas_strong[1]$membership == 1])) # V:207, E:1057
ecount <- gsize(gwas_graph)
vcount <- vcount(gwas_graph)
str_ecount <- gsize(gwas_graph_sub)
str_vcount <- vcount(gwas_graph_sub)
######
# average shortest path
huri_average_SP <- mean_distance(huri_g, directed = FALSE)
gwas_average_SP <- mean_distance(gwas_graph, directed = FALSE)
gwas_16g_aspl <- avgShortInHuRI(gwas_all)
######
# interactor of GWAS hit
gwas_ori_eb <- cluster_edge_betweenness(gwas_graph_sub, directed = FALSE)
eb_size <- length(sizes(gwas_ori_eb))
gwas_ori_ocg <- getOCG.clusters(as_data_frame(gwas_graph_sub), init.class.sys = 1, cent.class.sys = 0)
ocg_size <- length(gwas_ori_ocg$clustsizes)
# community size distribution (1st - 2nd)
eb_size_1_2 <- sort(sizes(gwas_ori_eb), decreasing = TRUE)[c(1, 2)]
eb_size_df <- eb_size_1_2[1] - eb_size_1_2[2]
ocg_size_1_2 <- sort(gwas_ori_ocg$clustsizes, decreasing = TRUE)[c(1, 2)]
ocg_size_df <- ocg_size_1_2[1] - ocg_size_1_2[2]

######
# community profiles calculation
# 1. edge and node counts in entire network and largely connected component, number of communities via edge-betweenness and OCG algorithms
all <- list()
for (i in 1:100) {
    # all[[i]] <- as.data.frame(countRandomize(gwas_all))
    all[[i]] <- countRandomize2(gwas_all)
    print(i)
    # one possible reason `mcreplicate` is not useable here is because the OCG algorithm will generate a temp interaction matrix
}
all_df <- as.data.frame(t(do.call(cbind, all)))
names(all_df) <- c("ecount", "vcount", "str_ecount", "str_vcount", "eb_size", "ocg_size")
# 2. average shortest path length
aspl_16g_huri <- c()
aspl_16g_huri <- c(aspl_16g_huri, mcreplicate(11000, avgShortInHuRI(gwas_all), mc.cores = detectCores())) # with unknown reason there are some infinity values
all_df$aspl <- aspl_16g_huri[aspl_16g_huri != Inf][c(1:10000)]
# plotting rewiring result/graph
pdf(file = "GWAS_loci_OCG_original.pdf")
plotOCGraph(gwas_ori_ocg)
df <- as.data.frame(gwas_ori_ocg$clustsizes)
names(df) <- "cluster size"
mtext(side = 3, "GWAS Observed", line = 1)
addtable2plot(-1, 0.7, df, display.rownames = TRUE, hlines = TRUE, cex = 0.6, bty = "o")
dev.off()
######
# plotting results of community profiles after rewiring
pdf(file = "~/Documents/INET-work/virus_network/statistic_results/GWAS/GWAS_fraction.pdf", width = 3, height = 3)
par(mgp = c(2, 0.7, 0), ps = 8)
# ecount
dens_gwas <- hist(all_df$ecount, plot = FALSE)
plot(dens_gwas, xlim = c(350, 1200), col = rgb(0.75, 0.75, 0.75, 1/2), border = NA, las = 1, yaxt = "n", xlab = "Number of edges", ylab = "Frequency density", main = "", cex.sub = 0.5)
mytitle <- "Subnetwork of COVID19 GWAS analysis"
mysubtitle <- "Number of edges in subnetwork of 16 identified locus"
mtext(side = 3, line = 1, cex = 1, mytitle)
mtext(side = 3, line = 0.5, cex = .7, mysubtitle)
axis(side = 2, at = seq(0, 2500, 500), labels = seq(0, 2500, 500)/10000, las = 1)
arrows(ecount, 300, ecount, 0, col = "#922687", lwd = 2, length = 0.1)
text(ecount - 100, 400, "p < 1e-4", cex = 0.4, pos = 4)
# vcount
dens_gwas <- hist(all_df$vcount, plot = FALSE)
plot(dens_gwas, col = rgb(0.75, 0.75, 0.75, 1/2), border = NA, las = 1, yaxt = "n", xlab = "Number of nodes", ylab = "Frequency density", main = "", cex.sub = 0.5)
mytitle <- "Subnetwork of COVID19 GWAS analysis"
mysubtitle <- "Number of nodes in subnetwork of 16 identified locus"
mtext(side = 3, line = 1, cex = 1, mytitle)
mtext(side = 3, line = 0.5, cex = .7, mysubtitle)
axis(side = 2, at = seq(0, 2500, 500), labels = seq(0, 2500, 500)/10000, las = 1)
arrows(vcount, 300, vcount, 0, col = "#922687", lwd = 2, length = 0.1)
text(vcount - 3, 400, "p = 0.8591", cex = 0.4, pos = 4)
# strongest ecount
dens_gwas <- hist(all_df$str_ecount, plot = FALSE)
plot(dens_gwas, xlim = c(350, 1100), col = rgb(0.75, 0.75, 0.75, 1/2), border = NA, las = 1, yaxt = "n", xlab = "Number of edges in HuRI", ylab = "Frequency density", main = "", cex.sub = 0.5)
mytitle <- "Subnetwork of COVID19 GWAS analysis"
mysubtitle <- "Number of edges in subnetwork of 16 identified locus"
mtext(side = 3, line = 1, cex = 1, mytitle)
mtext(side = 3, line = 0.5, cex = .7, mysubtitle)
axis(side = 2, at = seq(0, 2500, 500), labels = seq(0, 2500, 500)/10000, las = 1)
arrows(str_ecount, 300, str_ecount, 0, col = "#922687", lwd = 2, length = 0.1)
text(str_ecount - 50, 400, "p < 1e-4", cex = 0.4, pos = 4)
# strongest vcount
dens_gwas <- hist(all_df$str_vcount, plot = FALSE)
plot(dens_gwas, xlim = c(120, 220), col = rgb(0.75, 0.75, 0.75, 1/2), border = NA, las = 1, yaxt = "n", xlab = "Number of nodes in HuRI", ylab = "Frequency density", main = "", cex.sub = 0.5)
mytitle <- "Subnetwork of COVID19 GWAS analysis"
mysubtitle <- "Number of nodes in subnetwork of 16 identified locus"
mtext(side = 3, line = 1, cex = 1, mytitle)
mtext(side = 3, line = 0.5, cex = .7, mysubtitle)
axis(side = 2, at = seq(0, 2500, 500), labels = seq(0, 2500, 500)/10000, las = 1)
arrows(str_vcount, 300, str_vcount, 0, col = "#922687", lwd = 2, length = 0.1)
text(str_vcount - 5, 400, "p < 1e-04", cex = 0.4, pos = 4)
# number of clusters via edge-betweenness
dens_gwas <- hist(all_df$eb_size, breaks = 100, plot = FALSE)
plot(dens_gwas, xlim = c(0, 100), col = rgb(0.75, 0.75, 0.75, 1/2), border = NA, las = 1, yaxt = "n", xlab = "Number of clusters via EB", ylab = "Frequency density", main = "", cex.sub = 0.5)
mytitle <- "Subnetwork of COVID19 GWAS analysis"
mysubtitle <- "Number of clusters in subnetwork of 16 identified locus"
mtext(side = 3, line = 1, cex = 1, mytitle)
mtext(side = 3, line = 0.5, cex = .7, mysubtitle)
axis(side = 2, at = seq(0, 350, 50), labels = seq(0, 350, 50)/10000, las = 1)
arrows(eb_size, 50, eb_size, 0, col = "#922687", lwd = 2, length = 0.1)
text(eb_size - 5, 80, "p = 0.9999", cex = 0.4, pos = 4)
# number of clusters via OCG
# dens_gwas <- hist(all_df$ocg_size, plot = FALSE)
barplot(prop.table(table(all_df$ocg_size)),
    xlim = c(0, 14), col = rgb(0.75, 0.75, 0.75, 1/2), border = NA, las = 1, xaxt = "n", space = 0, xlab = "Number of clusters via EB", ylab = "Frequency density", main = "", cex.sub = 0.5)
mytitle <- "Subnetwork of COVID19 GWAS analysis"
mysubtitle <- "Number of clusters in subnetwork of 16 identified locus"
mtext(side = 3, line = 1, cex = 1, mytitle)
mtext(side = 3, line = 0.5, cex = .7, mysubtitle)
axis(side = 1, at = seq(0, 14, 2) + 0.5, labels = seq(2, 16, 2), line = 0.3)
arrows(0.5, 0.02, 0.5, 0, col = "#922687", lwd = 2, length = 0.1)
text(0, 0.04, "p = 1", cex = 0.4, pos = 4)
dev.off()
# boxplot of community profiles after rewiring
pdf(file = "~/Documents/INET-work/virus_network/statistic_results/GWAS/GWAS_fraction_boxplot.pdf", width = 7, height = 3)
par(mar = c(8, 4, 2, 2), mgp = c(3, 0.7, 0), ps = 10, mfrow = c(1, 4))
names <- c("entire network", "largely connected")
boxplot(all_df[, c(1, 3)], las = 2, ylab = "", xlab = "", main = "edge\ndistribution", names = names, frame = F)
title(ylab = "count", line = 2)
boxplot(all_df[, c(2, 4)], las = 2, ylab = "", xlab = "", main = "node\ndistribution", names = names, frame = F)
title(ylab = "count", line = 2)
boxplot(all_df[, c(5, 6)], las = 2, ylab = "", xlab = "", main = "community\nnumber", names = c("edge-betweenness", "OCG"), frame = F)
title(ylab = "count", line = 2)
boxplot(all_df[, 7], ylab = "", xlab = "", main = "average shortest\n path length", names = "HuRI", frame = F)
title(ylab = "average shortest path length", line = 2)
title(xlab = "randomized HuRI", line = 1)
dev.off()
# visualization of community
for (i in 1:20) {
    gwas_rand1 <- random_gwas(gwas_all)
    gwas_rand1_strong <- components(gwas_rand1, mode = "strong")
    biggest_cluster_id <- which.max(gwas_rand1_strong$csize)
    vert_ids <- V(gwas_rand1)[gwas_rand1_strong$membership == biggest_cluster_id]
    gwas_rand1_strong_graph <- induced_subgraph(gwas_rand1, vert_ids)
    gwas_rand1_eb <- cluster_edge_betweenness(gwas_rand1_strong_graph, directed = FALSE)
    gwas_rand1_ocg <- getOCG.clusters(as_data_frame(gwas_rand1), init.class.sys = 1, cent.class.sys = 0)
    pdf(file = paste0("/tmp/GAWS_loci_OCG_x1000_noLoop_", i, ".pdf"))
    plotOCGraph(gwas_rand1_ocg)
    df <- as.data.frame(gwas_rand1_ocg$clustsizes)
    names(df) <- "cluster size"
    mtext(side = 3, paste0("randomization analysis #", i), line = 1)
    addtable2plot(-1, 0.7, df, display.rownames = TRUE, hlines = TRUE, cex = 0.6)
    dev.off()
}
######
# degree preserving HuRI and extract network with 3 phenotype associated GWAS hits
# 1. network profile
ctcl_all <- list()
for (i in 1:1000) {
    # ctcl_all[[i]] <- as.data.frame(countRandomize(ctcl))
    ctcl_all[[i]] <- countRandomize2(ctcl)
}
ctcl_all_df <- t(as.data.frame(do.call(cbind, ctcl_all)))
hosp_all <- list()
for (i in 1:10000) {
    hosp_all[[i]] <- as.data.frame(countRandomize(hosp))
}
hosp_all_df <- t(as.data.frame(do.call(cbind, hosp_all)))
infct_all <- list()
for (i in 1:10000) {
    infct_all[[i]] <- as.data.frame(countRandomize(infct))
}
infct_all_df <- t(as.data.frame(do.call(cbind, infct_all)))
# plotting
# boxplot of community profiles after rewiring
pdf(file = "~/Documents/INET-work/virus_network/statistic_results/GWAS/GWAS_phenotypes_permutation.pdf", width = 3, height = 3)
par(mgp = c(2, 0.7, 0), ps = 8)
# critical illness, ecount
dens_gwas <- hist(ctcl_all_df[, 1], plot = FALSE)
plot(dens_gwas, xlim = c(0, 120), col = rgb(0.75, 0.75, 0.75, 1/2), border = NA, las = 1, yaxt = "n", xlab = "Number of edges", ylab = "Frequency density", main = "", cex.sub = 0.5)
mytitle <- "Subnetwork of COVID19 GWAS analysis"
mysubtitle <- "Number of edges, critical illness"
mtext(side = 3, line = 1, cex = 1, mytitle)
mtext(side = 3, line = 0.5, cex = .7, mysubtitle)
axis(side = 2, at = seq(0, 3000, 500), labels = seq(0, 3000, 500)/10000, las = 1)
arrows(110, 500, 110, 0, col = "#922687", lwd = 2, length = 0.1)
text(110 - 10, 600, "p = 0.0006", cex = 0.4, pos = 4)
# critical illness, largely connect ecount
dens_gwas <- hist(ctcl_all_df[, 3], plot = FALSE)
plot(dens_gwas, xlim = c(0, 120), col = rgb(0.75, 0.75, 0.75, 1/2), border = NA, las = 1, yaxt = "n", xlab = "Number of edges", ylab = "Frequency density", main = "", cex.sub = 0.5)
mytitle <- "Subnetwork of COVID19 GWAS analysis"
mysubtitle <- "Number of edges, largely connect, critical illness"
mtext(side = 3, line = 1, cex = 1, mytitle)
mtext(side = 3, line = 0.5, cex = .7, mysubtitle)
axis(side = 2, at = seq(0, 3000, 500), labels = seq(0, 3000, 500)/10000, las = 1)
arrows(109, 500, 109, 0, col = "#922687", lwd = 2, length = 0.1)
text(109 - 10, 600, "p = 0.0003", cex = 0.4, pos = 4)
# hospitalization, ecount
dens_gwas <- hist(hosp_all_df[, 1], plot = FALSE)
plot(dens_gwas, xlim = c(0, 1100), col = rgb(0.75, 0.75, 0.75, 1/2), border = NA, las = 1, yaxt = "n", xlab = "Number of edges", ylab = "Frequency density", main = "", cex.sub = 0.5)
mytitle <- "Subnetwork of COVID19 GWAS analysis"
mysubtitle <- "Number of edges, hospitalization"
mtext(side = 3, line = 1, cex = 1, mytitle)
mtext(side = 3, line = 0.5, cex = .7, mysubtitle)
axis(side = 2, at = seq(0, 2500, 500), labels = seq(0, 2500, 500)/10000, las = 1)
arrows(1011, 500, 1011, 0, col = "#922687", lwd = 2, length = 0.1)
text(1011 - 200, 600, "p < 1e-4", cex = 0.4, pos = 4)
# hospitalization, largely connect, ecount
dens_gwas <- hist(hosp_all_df[, 3], plot = FALSE)
plot(dens_gwas, xlim = c(0, 1100), col = rgb(0.75, 0.75, 0.75, 1/2), border = NA, las = 1, yaxt = "n", xlab = "Number of edges", ylab = "Frequency density", main = "", cex.sub = 0.5)
mytitle <- "Subnetwork of COVID19 GWAS analysis"
mysubtitle <- "Number of edges, hospitalization"
mtext(side = 3, line = 1, cex = 1, mytitle)
mtext(side = 3, line = 0.5, cex = .7, mysubtitle)
axis(side = 2, at = seq(0, 2500, 500), labels = seq(0, 2500, 500)/10000, las = 1)
arrows(1010, 500, 1010, 0, col = "#922687", lwd = 2, length = 0.1)
text(1010 - 200, 600, "p < 1e-4", cex = 0.4, pos = 4)
# infection, ecount
dens_gwas <- hist(infct_all_df[, 1], plot = FALSE)
plot(dens_gwas, xlim = c(0, 50), col = rgb(0.75, 0.75, 0.75, 1/2), border = NA, las = 1, yaxt = "n", xlab = "Number of edges", ylab = "Frequency density", main = "", cex.sub = 0.5)
mytitle <- "Subnetwork of COVID19 GWAS analysis"
mysubtitle <- "Number of edges, infection"
mtext(side = 3, line = 1, cex = 1, mytitle)
mtext(side = 3, line = 0.5, cex = .7, mysubtitle)
axis(side = 2, at = seq(0, 1500, 250), labels = seq(0, 1500, 250)/10000, las = 1)
arrows(44, 250, 44, 0, col = "#922687", lwd = 2, length = 0.1)
text(44 - 5, 350, "p < 1e-4", cex = 0.4, pos = 4)
# infection, ecount
dens_gwas <- hist(infct_all_df[, 3], plot = FALSE)
plot(dens_gwas, xlim = c(0, 50), col = rgb(0.75, 0.75, 0.75, 1/2), border = NA, las = 1, yaxt = "n", xlab = "Number of edges", ylab = "Frequency density", main = "", cex.sub = 0.5)
mytitle <- "Subnetwork of COVID19 GWAS analysis"
mysubtitle <- "Number of edges, infection"
mtext(side = 3, line = 1, cex = 1, mytitle)
mtext(side = 3, line = 0.5, cex = .7, mysubtitle)
axis(side = 2, at = seq(0, 1500, 250), labels = seq(0, 1500, 250)/10000, las = 1)
arrows(39, 250, 39, 0, col = "#922687", lwd = 2, length = 0.1)
text(39 - 5, 350, "p = 0.0002", cex = 0.4, pos = 4)

dev.off()
# 2. average shortest path length
aspl_ctcl_huri <- c()
aspl_ctcl_huri <- c(aspl_ctcl_huri, mcreplicate(11000, avgShortInHuRI(ctcl), mc.cores = detectCores())) # with unknown reason there are some infinity values
