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

huriRewire <- function(remove.loops = FALSE, ...) {
    huri_re <- rewire(huri_g, keeping_degseq(niter = gsize(huri_g) * 10))
    huri_sim <- simplify(huri_re, remove.loops = remove.loops)
    return(huri_sim)
}

huriRewireDataset <- function(node, remove.loops) {
    count <- c()
    huri_re <- huriRewire(remove.loops)
    merged <- combineNetwork(huri_re, node)
    # merged_inHuSCI
    count <- c(count, as.numeric(table(V(merged)$name %in% husci_sym)["TRUE"]))
    # merged_inGordon
    count <- c(count, as.numeric(table(V(merged)$name %in% gordon_sym)["TRUE"]))
    # merged_inStukalov
    count <- c(count, as.numeric(table(V(merged)$name %in% stukalov_sym)["TRUE"]))
    count <- c(count, gsize(merged))
    return(count)
}

######
# load dataset
huri <- read.xlsx("../data/extended_table/Extended_Table_2_PPIs.xlsx", sheet = "HuRI")
husci <- read.csv("../data/HuSCI_node.csv", header = TRUE)
gwas <- read.csv("../data/GWAS_hits.csv", header = T)
gordon <- read.xlsx("../data/extended_table/Extended_Table_2_PPIs.xlsx", sheet = "Gordon")
stukalov <- read.xlsx("../data/extended_table/Extended_Table_2_PPIs.xlsx", sheet = "Stukalov")

######
# 1. HuRI graph generation
huri_symbol <- huri[, c(5:6)]
huri_g_ori <- graph_from_data_frame(huri_symbol, directed = FALSE) # V:8274, E:52573
huri_g <- simplify(huri_g_ori, remove.loops = FALSE) # V:8274, E:52558
# protein list filter
husci_sym <- husci$node
husci_huri <- V(huri_g)$name[V(huri_g)$name %in% husci_sym] # HuSCI in HuRI whole

# GWAS hit in HuRI
gwas_huri <- gwas$name[!is.na(gwas$ctl == 1)]
# Gordon and Stukalov in HuRI
gordon_sym <- unique(gordon$PreyGene)
gordon_huri <- V(huri_g)$name[V(huri_g)$name %in% gordon_sym]

stukalov_sym <- unique(stukalov$human)
stukalov_huri <- V(huri_g)$name[V(huri_g)$name %in% stukalov_sym]

######
# 2. interactor of GWAS hit
gwas_hit_1st <- make_ego_graph(huri_g, nodes = gwas_huri, order = 1, mode = "all")

######
# 3. **rewiring analysis of HuRI**, to see if the HuSCI viral target is significant.
# subnetwork of GWAS hit from HuRI
# inherit from above code
gwas_all_list_df <- lapply(gwas_hit_1st, as_data_frame)
gwas_all_df <- do.call(rbind, gwas_all_list_df)
gwas_all_g_merge <- graph_from_data_frame(gwas_all_df, directed = FALSE)
# to have interaction between 1st interactors
gwas_all_final <- simplify(induced_subgraph(huri_g, names(V(gwas_all_g_merge))), remove.loops = F)

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
gwas_rand_r2 <- c(gwas_rand_r2, mcreplicate(10000, huriRewireDataset(gwas_huri, FALSE), mc.cores = detectCores()))
gwas_rand_r2[is.na(gwas_rand_r2)] <- 0
gwas_rand_df_r2 <- data.frame(matrix(gwas_rand_r2, ncol = 4, byrow = T))
names(gwas_rand_df_r2) <- c("HuSCI_viral_target", "Gordon_viral_target", "Stukalov_viral_target", "interactions")

sign_husci <-
######
# plot
pdf(file = "~/Documents/INET-work/virus_network/figure_results/GWAS/GWAS_3dataset.pdf", width = 3, height = 3)
par(mgp = c(2, 0.7, 0), ps = 8)
# HuSCI viral target in GWAS subnetwork
dens_gwas <- hist(gwas_rand_df_r2$HuSCI_viral_target, plot = FALSE, right = F)
plot(dens_gwas, xlim = c(0, 25), col = rgb(0.75, 0.75, 0.75, 1/2), border = NA, las = 1, yaxt = ifelse(sum(dens_gwas$density) != 1, "n", "s"), freq = ifelse(sum(dens_gwas$density) != 1, TRUE, FALSE), xlab = "Number of viral targets", ylab = "Frequency density", main = "", cex.sub = 0.5)
mytitle <- "Subnetwork of COVID19 GWAS"
mysubtitle <- "# HuSCI viral targets"
mtext(side = 3, line = 1, cex = 1, mytitle)
mtext(side = 3, line = 0.5, cex = .7, mysubtitle)
if(sum(dens_gwas$density) != 1) {
    axis(side = 2, at = seq(0, 2500, 500), labels = seq(0, 2500, 500)/10000, las = 1)
} else {
    next
}
arrows(gwas_all_husci_length, 300, gwas_all_husci_length, 0, col = "#922687", lwd = 2, length = 0.1)
text(gwas_all_husci_length, 400, paste0("observed = 20 \np = ", table(gwas_rand_df_r2[, "HuSCI_viral_target"] >= 20)["TRUE"]/10000), cex = 0.4, pos = 4)
# Gordon viral target in GWAS subnetwork
dens_gwas <- hist(gwas_rand_df_r2$Gordon_viral_target, plot = FALSE, right = F)
plot(dens_gwas, xlim = c(0, 25), col = rgb(0.75, 0.75, 0.75, 1/2), border = NA, las = 1, yaxt = ifelse(sum(dens_gwas$density) != 1, "n", "s"), freq = ifelse(sum(dens_gwas$density) != 1, TRUE, FALSE), xlab = "Number of viral targets", ylab = "Frequency density", main = "", cex.sub = 0.5)
mytitle <- "Subnetwork of COVID19 GWAS"
mysubtitle <- "# Gordon et al. viral targets"
mtext(side = 3, line = 1, cex = 1, mytitle)
mtext(side = 3, line = 0.5, cex = .7, mysubtitle)
if(sum(dens_gwas$density) != 1) {
    axis(side = 2, at = seq(0, 2500, 500), labels = seq(0, 2500, 500)/10000, las = 1)
}
arrows(gwas_all_gordon_length, 0.03, gwas_all_gordon_length, 0, col = "#922687", lwd = 2, length = 0.1)
text(gwas_all_gordon_length, 0.04, paste0("observed = 7 \np = ", table(gwas_rand_df_r2[, "Gordon_viral_target"] >= 7)["TRUE"]/10000), cex = 0.4, pos = 4)
# Stukalov viral target in GWAS subnetwork
dens_gwas <- hist(gwas_rand_df_r2$Stukalov_viral_target, plot = FALSE, right = F)
plot(dens_gwas, xlim = c(0, 25), col = rgb(0.75, 0.75, 0.75, 1/2), border = NA, las = 1, yaxt = ifelse(sum(dens_gwas$density) != 1, "n", "s"), freq = ifelse(sum(dens_gwas$density) != 1, TRUE, FALSE), xlab = "Number of viral targets", ylab = "Frequency density", main = "", cex.sub = 0.5)
mytitle <- "Subnetwork of COVID19 GWAS"
mysubtitle <- "# Stukalov et al. viral targets"
mtext(side = 3, line = 1, cex = 1, mytitle)
mtext(side = 3, line = 0.5, cex = .7, mysubtitle)
if(sum(dens_gwas$density) != 1) {
    axis(side = 2, at = seq(0, 2500, 500), labels = seq(0, 2500, 500)/10000, las = 1)
}
arrows(gwas_all_stukalov_length, 0.03, gwas_all_stukalov_length, 0, col = "#922687", lwd = 2, length = 0.1)
text(gwas_all_stukalov_length, 0.04, paste0("observed = 2 \np = ", table(gwas_rand_df_r2[, "Stukalov_viral_target"] >= 2)["TRUE"]/10000), cex = 0.4, pos = 4)
dev.off()

df <- data.frame(HuRI = V(huri_g)$name, inGWASsubnetwork = V(huri_g)$name %in% gwas$name, inHuSCI = V(huri_g)$name %in% husci_sym, inGordon = V(huri_g)$name %in% gordon_sym, inGordoninGWAS = V(huri_g)$name %in% gordon_sym[gordon_sym %in% gwas$name], inStukalov = V(huri_g)$name %in% stukalov_sym, inStukalovinGWAS = V(huri_g)$name %in% stukalov_sym[stukalov_sym %in% gwas$name])
df_raw <- list(HuSCI = husci_sym, Gordon = gordon_sym, Stukalov = stukalov_sym)
attributes(df_raw) <- list(names = names(df_raw), row.names = 1:max(length(husci_sym), length(gordon_sym), length(stukalov_sym)), class = 'data.frame') # ref: https://stackoverflow.com/questions/7196450/create-a-data-frame-of-unequal-lengths
write.xlsx(df, file = "~/Documents/INET-work/virus_network/statistic_results/GWAS/3dataset_df.xlsx", quote = TRUE, overwrite = TRUE)
write.table(husci_sym, file = "/tmp/husci_sym.tsv", sep = "\t", row.names = F, quote = F)
write.table(gordon_sym, file = "/tmp/gordon_sym.tsv", sep = "\t", row.names = F, quote = F)
write.table(stukalov_sym, file = "/tmp/stukalov_sym.tsv", sep = "\t", row.names = F, quote = F)

######
# display degree of viral targets in HuRI
husci_deg <- data.frame(degree(huri_g, v = husci_huri))
names(husci_deg) <- "degree"
gordon_deg <- data.frame(degree(huri_g, v = gordon_huri))
names(gordon_deg) <- "degree"
stukalov_deg <- data.frame(degree(huri_g, v = stukalov_huri))
names(stukalov_deg) <- "degree"

hist(husci_deg$degree, xlab = "degree", main = "Degree of HuSCI proteins in HuRI")
addtable2plot(100, 50, summary(husci_deg), vlines = TRUE, bty = "l", cex = 2)

hist(gordon_deg$degree, xlab = "degree", main = "Degree of Gordon proteins in HuRI")
addtable2plot(150, 50, summary(gordon_deg), vlines = TRUE, bty = "l", cex = 2)

hist(stukalov_deg$degree, xlab = "degree", main = "Degree of Stukalov proteins in HuRI")
addtable2plot(70, 100, summary(stukalov_deg), vlines = TRUE, bty = "l", cex = 2)

wb <- createWorkbook()
addWorksheet(wb, "HuSCI")
writeData(wb, "HuSCI", husci_deg, rowNames = TRUE)

addWorksheet(wb, "Gordon et al")
writeData(wb, "Gordon et al", gordon_deg, rowNames = TRUE)

addWorksheet(wb, "Stukalov et al")
writeData(wb, "Stukalov et al", stukalov_deg, rowNames = TRUE)

saveWorkbook(wb, "~/Documents/INET-work/virus_network/statistic_results/GWAS/3dataset_degree_HuRI.xlsx", overwrite = TRUE)
