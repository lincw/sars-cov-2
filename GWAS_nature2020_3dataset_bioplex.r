# COVID-19 GWAS hit analysis
# plus Gordon and Stukalov dataset
# Lin Chung-wen
# Date: 28.07.2021 **23:49**

setwd("/Volumes/GoogleDrive/My Drive/VirHostome_CW/GitHub/src")
######
# load package
library(igraph)
library(rethinking)
library(gprofiler2)
library(gplots)
library(openxlsx)
library(plotrix) # add table to plot

source("combineNetwork.r")

bioplexRewire <- function(remove.loops = FALSE, ...) {
    bioplex_re <- rewire(bioplex_g, keeping_degseq(niter = gsize(bioplex_g) * 10))
    bioplex_sim <- simplify(bioplex_re, remove.loops = remove.loops)
    return(bioplex_sim)
}
bioplexRewireDataset <- function(node, remove.loops = FALSE) {
    count <- c()
    re <- bioplexRewire(remove.loops)
    merged <- combineNetwork(re, node)
    # merged_inHuSCI
    count <- c(count, as.numeric(table(V(merged)$name %in% husci_sym)["TRUE"]))
    # merged_inGordon
    count <- c(count, as.numeric(table(V(merged)$name %in% gordon_sym)["TRUE"]))
    # merged_inStukalov
    count <- c(count, as.numeric(table(V(merged)$name %in% stukalov_sym)["TRUE"]))
    count <- c(count, gsize(merged))
    return(count)
}
plotHist <- function(value, title, phenotype, length, xmax, y1, y2, density = TRUE) {
    if (density == TRUE) {
        dens_gwas <- hist(value, breaks = c(0:(max(value) + 1)), plot = FALSE, right = F)
        plot(dens_gwas, xlim = c(0, xmax), col = rgb(0.75, 0.75, 0.75, 1/2), border = NA, las = 1, xaxt = "n", freq = FALSE, xlab = "Number of viral targets", ylab = "Frequency", main = "", cex.sub = 0.5)
    } else {
        dens_gwas <- hist(value, plot = FALSE, right = F)
        plot(dens_gwas, xlim = c(140, xmax), col = rgb(0.75, 0.75, 0.75, 1/2), border = NA, las = 1, xaxt = "n", yaxt = "n", xlab = "Number of viral targets", ylab = "Frequency", main = "", cex.sub = 0.5)
    }
    mytitle <- paste0("COVID19 GWAS subnetwork\n(", phenotype, ")\nviral targets in ", title)
    mtext(side = 3, line = 1, cex = 1, mytitle)
    mtext(side = 3, line = 0.2, cex = .8, "subnetwork extracted from BioPlex3.0")
    if (density == TRUE) {
        axis(side = 1, at = seq(0, xmax, by = 5) + 0.5, labels = seq(0, xmax, by = 5))
    } else {
        axis(side = 1, at = seq(140, 300, by = 20), labels = seq(140, 300, by = 20))
        axis(side = 2, at = seq(0, 2500, by = 500), labels = seq(0, 0.25, by = 0.05), las = 2)
    }
    arrows(length + 0.5, y1, length + 0.5, 0, col = "#922687", lwd = 2, length = 0.1)
    text(median(value) + 4, ifelse(density == TRUE, max(dens_gwas$counts / 10000), max(dens_gwas$counts)), paste0("median = ", median(value)), col = "grey", cex = 0.5)
    text(length - 2, y2, paste0("observed = ", length, "\np = ", table(value >= length)["TRUE"]/10000), cex = 0.4, pos = 4)
}
######
# load dataset
bioplex <- read.delim("../data/extended_table/BioPlex.3.0_edge.tsv", header = T)
husci <- read.csv("../data/HuSCI_node.csv", header = TRUE)
gwas <- read.csv("~/Documents/INET-work/virus_network/references/GWAS/Genetic\ mechanisms\ of\ critical\ illness\ in\ COVID-19/table1.csv", header = T)
gordon <- read.xlsx("../data/extended_table/Extended_Table_2_PPIs.xlsx", sheet = "Gordon")
stukalov <- read.xlsx("../data/extended_table/Extended_Table_2_PPIs.xlsx", sheet = "Stukalov")

gwas <- c(gwas$Locus[c(1:4, 6:8)], "OAS1", "OAS2", "OAS3")
######
# 1. BioPlex graph generation
bioplex_symbol <- bioplex[, c(5:6)]
bioplex_g_ori <- graph_from_data_frame(bioplex_symbol, directed = FALSE) # V:13957, E:118162
bioplex_g <- simplify(bioplex_g_ori, remove.loops = FALSE) # V:13957, E:118162
# protein list filter
husci_sym <- husci$node
husci_bioplex <- V(bioplex_g)$name[V(bioplex_g)$name %in% husci_sym] # HuSCI in BioPlex whole, V:132

# GWAS hit in BioPlex
gwas_bioplex <- gwas[gwas %in% V(bioplex_g)$name] # V:7
# Gordon and Stukalov in BioPlex
gordon_sym <- unique(gordon$PreyGene)
gordon_bioplex <- V(bioplex_g)$name[V(bioplex_g)$name %in% gordon_sym] # V:346

stukalov_sym <- unique(stukalov$human)
stukalov_bioplex <- V(bioplex_g)$name[V(bioplex_g)$name %in% stukalov_sym] # V:723

######
# 2. interactor of GWAS hit
gwas_hit_1st <- make_ego_graph(bioplex_g, nodes = gwas_bioplex, order = 1, mode = "all")

######
# 3. **rewiring analysis of HuRI**, to see if the HuSCI viral target is significant.
# subnetwork of GWAS hit from HuRI
# inherit from above code
gwas_all_list_df <- lapply(gwas_hit_1st, as_data_frame)
gwas_all_df <- do.call(rbind, gwas_all_list_df)
gwas_all_g_merge <- graph_from_data_frame(gwas_all_df, directed = FALSE)
# to have interaction between 1st interactors
gwas_all_final <- simplify(induced_subgraph(bioplex_g, names(V(gwas_all_g_merge))), remove.loops = F)

# GWAS hit in HuSCI
gwas_all_husci <- V(gwas_all_final)$name[V(gwas_all_final)$name %in% husci_sym]
gwas_all_husci_length <- length(gwas_all_husci)

# GWAS hit in Gordon
gwas_all_gordon <- V(gwas_all_final)$name[V(gwas_all_final)$name %in% gordon_sym]
gwas_all_gordon_length <- length(gwas_all_gordon)

# GWAS hit in Stukalov
gwas_all_stukalov <- V(gwas_all_final)$name[V(gwas_all_final)$name %in% stukalov_sym]
gwas_all_stukalov_length <- length(gwas_all_stukalov)

######
# permutation analysis
gwas_rand_r2 <- c()
gwas_rand_r2 <- c(gwas_rand_r2, mcreplicate(10000, bioplexRewireDataset(gwas_bioplex, FALSE), mc.cores = detectCores()))
gwas_rand_r2[is.na(gwas_rand_r2)] <- 0
gwas_rand_df_r2 <- data.frame(matrix(gwas_rand_r2, ncol = 4, byrow = T))
names(gwas_rand_df_r2) <- c("HuSCI_viral_target", "Gordon_viral_target", "Stukalov_viral_target", "interactions")

######
# plot
pdf(file = "~/Documents/INET-work/virus_network/figure_results/GWAS/Nature2020_3dataset_bioPlex.pdf", width = 3, height = 3)
par(mgp = c(2, 0.7, 0), ps = 8)
# HuSCI viral target in GWAS subnetwork
plotHist(gwas_rand_df_r2$HuSCI_viral_target, "HuSCI", "criticall illness, 7 genes and 1st interactors", gwas_all_husci_length, 20, 0.05, 0.07)

# Gordon viral target in GWAS subnetwork
plotHist(gwas_rand_df_r2$Gordon_viral_target, "Gordon et al", "criticall illness, 7 genes and 1st interactors",gwas_all_gordon_length, 20, 0.03, 0.05)

# Stukalov viral target in GWAS subnetwork
plotHist(gwas_rand_df_r2$Stukalov_viral_target, "Stukalov et al", "criticall illness, 7 genes and 1st interactors",gwas_all_stukalov_length, 20, 0.03, 0.05)

# Interconnectivity result of GWAS+1 subnetwork in BioPlex3
plotHist(gwas_rand_df_r2$interactions, "Interaction", "GWAS+1 subnetwork in BioPlex", length(E(gwas_all_final)), 300, 300, 500, density = FALSE)
dev.off()

######
# display degree of viral targets in HuRI
husci_deg <- data.frame(degree(bioplex_g, v = husci_bioplex))
names(husci_deg) <- "degree"
gordon_deg <- data.frame(degree(bioplex_g, v = gordon_bioplex))
names(gordon_deg) <- "degree"
stukalov_deg <- data.frame(degree(bioplex_g, v = stukalov_bioplex))
names(stukalov_deg) <- "degree"

hist(husci_deg$degree, xlab = "degree", main = "Degree of HuSCI proteins in BioPlex")
addtable2plot(100, 50, summary(husci_deg), vlines = TRUE, bty = "l", cex = 2)

hist(gordon_deg$degree, xlab = "degree", main = "Degree of Gordon proteins in BioPlex")
addtable2plot(150, 50, summary(gordon_deg), vlines = TRUE, bty = "l", cex = 2)

hist(stukalov_deg$degree, xlab = "degree", main = "Degree of Stukalov proteins in BioPlex")
addtable2plot(70, 100, summary(stukalov_deg), vlines = TRUE, bty = "l", cex = 2)

wb <- createWorkbook()
addWorksheet(wb, "HuSCI")
writeData(wb, "HuSCI", husci_deg, rowNames = TRUE)

addWorksheet(wb, "Gordon et al")
writeData(wb, "Gordon et al", gordon_deg, rowNames = TRUE)

addWorksheet(wb, "Stukalov et al")
writeData(wb, "Stukalov et al", stukalov_deg, rowNames = TRUE)

saveWorkbook(wb, "~/Documents/INET-work/virus_network/statistic_results/GWAS/3dataset_degree_BioPlex.xlsx", overwrite = TRUE)

######
# save workarea data
save.image("~/Documents/INET-work/virus_network/statistic_results/GWAS/3dataset_degree_BioPlex.RData")
