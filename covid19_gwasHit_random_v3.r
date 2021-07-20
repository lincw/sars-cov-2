# COVID-19 GWAS hit analysis
# Lin Chung-win
## v3: based on the list from https://doi.org/10.1038/s41586-020-03065-y

######
# load package
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

huriRewireHusci <- function(node, remove.loops) {
    huri_re <- huriRewire(remove.loops)
    merged <- combineNetwork(huri_re, node)
    merged_inHuSCI <- as.numeric(table(V(merged)$name %in% husci_sym)["TRUE"])
    output <- c(merged_inHuSCI, gsize(merged))
    return(output)
}

# edit plot parameters
trace("plot.igraph", edit = T)
######
# load dataset
huri <- read.csv("~/Documents/INET-work/references/HuRI_binaryPPI/HuRI_Tong_withSymbol.csv", header = T)
husci <- read.csv("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/binary_node_1126.csv", header = TRUE)
gwas <- read.csv("~/Documents/INET-work/virus_network/statistic_results/GWAS/fromNature_node.csv", header = T)
######
# 1. HuRI graph generation
huri_symbol <- huri[, c(5:6)]
huri_g_ori <- graph_from_data_frame(huri_symbol, directed = FALSE) # V:8274, E:52573
huri_g <- simplify(huri_g_ori, remove.loops = FALSE) # V:8274, E:52558
# protein list filter
husci_sym <- husci[husci$group == "human", "node"]
husci_huri <- V(huri_g)$name[V(huri_g)$name %in% husci_sym] # HuSCI in HuRI whole
gwas_huri <- gwas$name[!is.na(gwas$ctl == 1)] # GWAS hit in HuRI

######
# 2. interactor of GWAS hit
gwas_hit_1st <- make_ego_graph(huri_g, nodes = gwas_huri, order = 1, mode = "all")

######
# 3. **rewiring analysis of HuRI**, to see if the HuSCI viral target is significant.
# load gwas loci info, with 3 phenotype (critical illness, hospitalization and infection)
# subnetwork of GWAS hit from HuRI
# inherite from above code
gwas_all_list_df <- lapply(gwas_hit_1st, as_data_frame)
gwas_all_df <- do.call(rbind, gwas_all_list_df)
gwas_all_g_merge <- graph_from_data_frame(gwas_all_df, directed = FALSE)
# to have interaction between 1st interactors
gwas_all_final <- simplify(induced_subgraph(huri_g, names(V(gwas_all_g_merge))), remove.loops = F) # V:169, E:1108
gwas_all_husci <- V(gwas_all_final)$name[V(gwas_all_final)$name %in% husci_sym] # 20 viral targets
gwas_all_husci_length <- length(gwas_all_husci)
# visualization of gwas subnetwork
bd <- ifelse(V(gwas_all_final)$name %in% husci_sym, "orange", NA)
label <- ifelse(V(gwas_all_final)$name %in% c(husci_sym, gwas_huri), V(gwas_all_final)$name, NA)
color <- ifelse(V(gwas_all_final)$name %in% gwas_huri, "red", "grey")
pdf("/tmp/GWAS_subnetwork.pdf")
plot(gwas_all_final, vertex.frame.color = bd, vertex.size = 3, vertex.label.dist = 1, vertex.label.color = "black", vertex.label = label, vertex.color = color, vertex.frame.width = 2)
dev.off()
# 3.1. do statistical analysis for all 5 GWAS hit candidate genes
gwas_rand_r2 <- c()
gwas_rand_r2 <- c(gwas_rand_r2, mcreplicate(10000, huriRewireHusci(gwas_huri, FALSE), mc.cores = detectCores()))
gwas_rand_df_r2 <- data.frame(matrix(gwas_rand_r2, ncol = 2, byrow = T))
names(gwas_rand_df_r2) <- c("viral_target", "interactions")
write.csv(gwas_rand_df_r2, file = "~/Documents/INET-work/virus_network/statistic_results/GWAS/10000_randomHuRI_v3_r2_0720.csv", row.names = FALSE)

# plotting
pdf("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/random_GWAS_viral_target_v3_r2_0720.pdf", width = 3, height = 3)
par(mgp = c(2, 0.7, 0), ps = 8)
dens_gwas <- hist(gwas_rand_df_r2[, "viral_target"], breaks = 10, plot = FALSE, right = F)
plot(dens_gwas, col = rgb(0.75, 0.75, 0.75, 1/2), xlim = c(0, 25), border = NA, las = 1, yaxt = "n", xlab = "Number of viral targets", ylab = "Frequency density", main = "", cex.sub = 0.5)
mytitle <- "COVID19 GWAS loci candidate genes"
mtext(side = 3, line = 1, cex = 1, mytitle)
axis(side = 2, at = seq(0, 2500, by = 500), labels = seq(0, 0.25, by = 0.05), las = 2)
arrows(gwas_all_husci_length, 500, gwas_all_husci_length, 10, col = "#922687", lwd = 2, length = 0.1)
text(gwas_all_husci_length - 2, 700, paste0("observed = 20 \np = ", table(gwas_rand_df_r2[, "viral_target"] >= 20)["TRUE"]/10000), cex = 0.4, pos = 4)

dens_gwas <- hist(gwas_rand_df_r2[, "interactions"], breaks = 35, plot = FALSE)
plot(dens_gwas, col = rgb(0.75, 0.75, 0.75, 1/2), border = NA, las = 1, xlim = c(200, 1200), yaxt = "n", xlab = "Number of interactions", ylab = "Frequency density", main = "", cex.sub = 0.5)
mytitle <- "COVID19 GWAS loci candidate genes"
mtext(side = 3, line = 1, cex = 1, mytitle)
axis(side = 2, at = seq(0, 1400, by = 200), labels = seq(0, 0.14, by = 0.02), las = 2)
arrows(gsize(gwas_all_final) + 0.5, 200, gsize(gwas_all_final) + 0.5, 20, col = "#922687", lwd = 2, length = 0.1)
text(gsize(gwas_all_final) - 200, 400, paste0("observed = ", gsize(gwas_all_final), "\np < 0.0001"), cex = 0.4, pos = 4)

dev.off()

######
# according to the raw data on the supplementary code (folder 05)
all_0720 <- read.xlsx("/Volumes/GoogleDrive/My Drive/Paper_VirHostome_CoV2/05_Source Data/Source Data permutation analyses.xlsx")
# plotting
pdf("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/random_GWAS_viral_target_v3_r2_fig4.pdf", width = 3, height = 3)
par(mgp = c(2, 0.7, 0), ps = 8)
dens_gwas <- hist(as.numeric(all_0720$X9[c(3:10002)]), breaks = 10, plot = FALSE, right = F)
plot(dens_gwas, col = rgb(0.75, 0.75, 0.75, 1 / 2), xlim = c(0, 25), border = NA, las = 1, yaxt = "n", xlab = "Number of viral targets", ylab = "Frequency density", main = "", cex.sub = 0.5)
mytitle <- "COVID19 GWAS loci candidate genes"
mtext(side = 3, line = 1, cex = 1, mytitle)
axis(side = 2, at = seq(0, 2500, by = 500), labels = seq(0, 0.25, by = 0.05), las = 2)
arrows(gwas_all_husci_length, 500, gwas_all_husci_length, 10, col = "#922687", lwd = 2, length = 0.1)
text(gwas_all_husci_length - 2, 700, paste0("observed = 20 \np = ", table(gwas_rand_df_r2[, "viral_target"] >= 20)["TRUE"] / 10000), cex = 0.4, pos = 4)

dens_gwas <- hist(as.numeric(all_0720$X8[c(3:10002)]), breaks = 25, plot = FALSE, right = F)
plot(dens_gwas, col = rgb(0.75, 0.75, 0.75, 1 / 2), border = NA, las = 1, xlim = c(200, 1200), yaxt = "n", xlab = "Number of interactions", ylab = "Frequency density", main = "", cex.sub = 0.5)
mytitle <- "COVID19 GWAS loci candidate genes"
mtext(side = 3, line = 1, cex = 1, mytitle)
axis(side = 2, at = seq(0, 1400, by = 200), labels = seq(0, 0.14, by = 0.02), las = 2)
arrows(gsize(gwas_all_final) + 0.5, 200, gsize(gwas_all_final) + 0.5, 20, col = "#922687", lwd = 2, length = 0.1)
text(gsize(gwas_all_final) - 200, 400, paste0("observed = ", gsize(gwas_all_final), "\np < 0.0001"), cex = 0.4, pos = 4)
dev.off()

