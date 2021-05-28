# COVID-19 GWAS locus analysis with HuRI degree preserving randomized network
# Lin Chung-wen
######
# load package
library(openxlsx)
library(linkcomm)
library(rethinking)
library(plotrix)
source("plotOCGGraph.r")
# trace("plot.igraph",edit=TRUE)
# as....
# 179             cp <- matrix(c(x0, y0, x0 + 0.08, y0 + 0.08, x0 + 0.08,
# 180                 y0 - 0.08, x0, y0), ncol = 2, byrow = TRUE)
# end...

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
# interactor of GWAS hit
gwas_ori_eb <- cluster_edge_betweenness(gwas_graph_sub, directed = FALSE)
eb_size <- length(sizes(gwas_ori_eb))
gwas_ori_ocg <- getOCG.clusters(as_data_frame(gwas_graph_sub), init.class.sys = 1, cent.class.sys = 0)
ocg_size <- length(gwas_ori_ocg$clustsizes)
######
# degree preserving randomized HuRI network
random_gwas <- function(node) {
    huri_random <- simplify(rewire(huri_g, keeping_degseq(niter = gsize(huri_g) * 10)), remove.loops = FALSE)
    gwas_random_g <- make_ego_graph(huri_random, nodes = node, order = 1, mode = "all")
    gwas_random_list_df <- lapply(gwas_random_g, as_data_frame)
    gwas_random_df <- do.call(rbind, gwas_random_list_df)
    gwas_random_g_merge <- graph_from_data_frame(gwas_random_df, directed = FALSE)
    gwas_random_final <- simplify(induced_subgraph(huri_g, names(V(gwas_random_g_merge))), remove.loops = FALSE) # V: 209, E: 1058
    return(gwas_random_final)
}
count_randomize <- function(...) {
    count <- c()
    gwas_rand1 <- random_gwas(gwas_all)
    count[1] <- gsize(gwas_rand1) # ecount
    count[2] <- vcount(gwas_rand1) # vcount
    gwas_rand1_strong <- components(gwas_rand1, mode = "strong")
    biggest_cluster_id <- which.max(gwas_rand1_strong$csize)
    vert_ids <- V(gwas_rand1)[gwas_rand1_strong$membership == biggest_cluster_id]
    gwas_rand1_strong_graph <- induced_subgraph(gwas_rand1, vert_ids)
    count[3] <- gsize(gwas_rand1_strong_graph) # str_ecount
    count[4] <- vcount(gwas_rand1_strong_graph) # str_vcount
    gwas_rand1_eb <- cluster_edge_betweenness(gwas_rand1_strong_graph, directed = FALSE)
    count[5] <- length(table(gwas_rand1_eb$membership)) # eb_size
    gwas_rand1_ocg <- getOCG.clusters(as_data_frame(gwas_rand1), init.class.sys = 1, cent.class.sys = 0, verbose = FALSE)
    count[6] <- length(gwas_rand1_ocg$clustsizes) # ocg_size
    return(count)
}
all <- list()
for (i in 1001:10000) {
    all[[i]] <- as.data.frame(count_randomize())
}
all_df <- as.data.frame(t(do.call(cbind, all)))
names(all_df) <- c("ecount", "vcount", "str_ecount", "str_vcount", "eb_size", "ocg_size")
# plotting rewiring result/graph
pdf(file = "/tmp/GWAS_loci_OCG_original.pdf")
plotOCGraph(gwas_ori_ocg)
df <- as.data.frame(gwas_ori_ocg$clustsizes)
names(df) <- "cluster size"
mtext(side = 3, "GWAS Observed", line = 1)
addtable2plot(-1, 0.7, df, display.rownames = TRUE, hlines = TRUE, cex = 0.6, bty = "o")
dev.off()

note <- "toCount"
pdf(file = "/tmp/GWAS_fraction.pdf", width = 3, height = 3)
par(mgp = c(2, 0.7, 0), ps = 8)
# ecount
dens_gwas <- hist(all_df$ecount, breaks = 15, plot = FALSE)
plot(dens_gwas, xlim = c(350, 1100), col = rgb(0.75, 0.75, 0.75, 1/2), freq = FALSE, border = NA, las = 1, xlab = "Number of edges in HuRI", ylab = "Frequency density", main = "", cex.sub = 0.5)
mytitle <- "Subnetwork of COVID19 GWAS analysis"
mysubtitle <- "Number of edges in subnetwork of 16 identified locus"
mtext(side = 3, line = 1, cex = 1, mytitle)
mtext(side = 3, line = 0.5, cex = .7, mysubtitle)
arrows(ecount, 0.0005, ecount, 0, col = "#922687", lwd = 2, length = 0.1)
text(ecount - 2, 0.001, paste0("p < ", note), cex = 0.4, pos = 4)
# vcount
dens_gwas <- hist(all_df$vcount, breaks = 15, plot = FALSE)
plot(dens_gwas, xlim = c(200, 220), col = rgb(0.75, 0.75, 0.75, 1/2), freq = FALSE, border = NA, las = 1, xlab = "Number of nodes in HuRI", ylab = "Frequency density", main = "", cex.sub = 0.5)
mytitle <- "Subnetwork of COVID19 GWAS analysis"
mysubtitle <- "Number of nodes in subnetwork of 16 identified locus"
mtext(side = 3, line = 1, cex = 1, mytitle)
mtext(side = 3, line = 0.5, cex = .7, mysubtitle)
arrows(vcount, 0.01, vcount, 0, col = "#922687", lwd = 2, length = 0.1)
text(vcount - 1, 0.02, paste0("p < ", note), cex = 0.4, pos = 4)
# strongest ecount
dens_gwas <- hist(all_df$str_ecount, breaks = 15, plot = FALSE)
plot(dens_gwas, xlim = c(350, 1100), col = rgb(0.75, 0.75, 0.75, 1/2), freq = FALSE, border = NA, las = 1, xlab = "Number of edges in HuRI", ylab = "Frequency density", main = "", cex.sub = 0.5)
mytitle <- "Subnetwork of COVID19 GWAS analysis"
mysubtitle <- "Number of edges in subnetwork of 16 identified locus"
mtext(side = 3, line = 1, cex = 1, mytitle)
mtext(side = 3, line = 0.5, cex = .7, mysubtitle)
arrows(str_ecount, 0.0005, str_ecount, 0, col = "#922687", lwd = 2, length = 0.1)
text(str_ecount - 2, 0.001, paste0("p < ", note), cex = 0.4, pos = 4)
# strongest vcount
dens_gwas <- hist(all_df$str_vcount, breaks = 15, plot = FALSE)
plot(dens_gwas, xlim = c(120, 220), col = rgb(0.75, 0.75, 0.75, 1/2), freq = FALSE, border = NA, las = 1, xlab = "Number of nodes in HuRI", ylab = "Frequency density", main = "", cex.sub = 0.5)
mytitle <- "Subnetwork of COVID19 GWAS analysis"
mysubtitle <- "Number of nodes in subnetwork of 16 identified locus"
mtext(side = 3, line = 1, cex = 1, mytitle)
mtext(side = 3, line = 0.5, cex = .7, mysubtitle)
arrows(str_vcount, 0.005, str_vcount, 0, col = "#922687", lwd = 2, length = 0.1)
text(str_vcount - 1, 0.01, paste0("p < ", note), cex = 0.4, pos = 4)
# number of clusters via edge-betweenness
dens_gwas <- hist(all_df$eb_size, breaks = 15, plot = FALSE)
plot(dens_gwas, xlim = c(0, 100), col = rgb(0.75, 0.75, 0.75, 1/2), freq = FALSE, border = NA, las = 1, xlab = "Number of clusters via EB", ylab = "Frequency density", main = "", cex.sub = 0.5)
mytitle <- "Subnetwork of COVID19 GWAS analysis"
mysubtitle <- "Number of clusters in subnetwork of 16 identified locus"
mtext(side = 3, line = 1, cex = 1, mytitle)
mtext(side = 3, line = 0.5, cex = .7, mysubtitle)
arrows(eb_size, 0.002, eb_size, 0, col = "#922687", lwd = 2, length = 0.1)
text(eb_size - 1, 0.004, paste0("p < ", note), cex = 0.4, pos = 4)
# number of clusters via OCG
dens_gwas <- hist(all_df$ocg_size, breaks = 15, plot = FALSE)
plot(dens_gwas, xlim = c(0, 15), col = rgb(0.75, 0.75, 0.75, 1/2), freq = FALSE, border = NA, las = 1, xlab = "Number of clusters via EB", ylab = "Frequency density", main = "", cex.sub = 0.5)
mytitle <- "Subnetwork of COVID19 GWAS analysis"
mysubtitle <- "Number of clusters in subnetwork of 16 identified locus"
mtext(side = 3, line = 1, cex = 1, mytitle)
mtext(side = 3, line = 0.5, cex = .7, mysubtitle)
arrows(ocg_size, 0.02, ocg_size, 0, col = "#922687", lwd = 2, length = 0.1)
text(ocg_size - 1, 0.04, paste0("p < ", note), cex = 0.4, pos = 4)
dev.off()

for (i in 1:20) {
    gwas_rand1 <- random_gwas(gwas_all)
    gwas_rand1_strong <- components(gwas_rand1, mode = "strong")
    biggest_cluster_id <- which.max(gwas_rand1_strong$csize)
    vert_ids <- V(gwas_rand1)[gwas_rand1_strong$membership == biggest_cluster_id]
    gwas_rand1_strong_graph <- induced_subgraph(gwas_rand1, vert_ids)
    gwas_rand1_eb <- cluster_edge_betweenness(gwas_rand1_strong_graph, directed = FALSE)
    gwas_rand1_ocg <- getOCG.clusters(as_data_frame(gwas_rand1), init.class.sys = 1, cent.class.sys = 0)
    pdf(file = paste0("/tmp/GAWS_loci_OCG_", i, ".pdf"))
    plotOCGraph(gwas_rand1_ocg)
    df <- as.data.frame(gwas_rand1_ocg$clustsizes)
    names(df) <- "cluster size"
    mtext(side = 3, paste0("randomization analysis #", i), line = 1)
    addtable2plot(-1, 0.7, df, display.rownames = TRUE, hlines = TRUE, cex = 0.6)
    dev.off()
}
