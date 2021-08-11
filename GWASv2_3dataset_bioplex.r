# COVID-19 GWAS hit analysis with BioPlex 3.0
# meta-analysis: https://www.nature.com/articles/s41586-021-03767-x
# Mapping the human genetic architecture of COVID-19
# plus Gordon and Stukalov dataset
# Lin Chung-wen
# Date: 11.08.2021

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

bioplexRewireMulti <- function(gwas, ctcl, hosp, infct, husci, gordon, stukalov, remove.loops = FALSE) {
    df <- c()
    re <- bioplexRewire(remove.loops) # rewire BioPlex
    merged <- combineNetwork(re, gwas) # get GWAS+1 subnetwork
    df <- c(df, table(V(merged)$name %in% husci)["TRUE"]) # get viral targets from HuSCI in GWAS+1 subnetwork
    df <- c(df, table(V(merged)$name %in% gordon)["TRUE"]) # get viral targets from Gordon in GWAS+1 subnetwork
    df <- c(df, table(V(merged)$name %in% stukalov)["TRUE"]) # get viral targets from Stukalov in GWAS+1 subnetwork
    df <- c(df, gsize(merged)) # get network size of GWAS+1 subnetwork

    merged <- combineNetwork(re, ctcl) # get GWAS+1 subnetwork only from critical illness candidates
    df <- c(df, table(V(merged)$name %in% husci)["TRUE"])
    df <- c(df, table(V(merged)$name %in% gordon)["TRUE"])
    df <- c(df, table(V(merged)$name %in% stukalov)["TRUE"])
    df <- c(df, gsize(merged))

    merged <- combineNetwork(re, hosp) # get GWAS+1 subnetwork only from hospitalized candidates
    df <- c(df, table(V(merged)$name %in% husci)["TRUE"])
    df <- c(df, table(V(merged)$name %in% gordon)["TRUE"])
    df <- c(df, table(V(merged)$name %in% stukalov)["TRUE"])
    df <- c(df, gsize(merged))

    merged <- combineNetwork(re, infct) # get GWAS+1 subnetwork only from reported infection candidates
    df <- c(df, table(V(merged)$name %in% husci)["TRUE"])
    df <- c(df, table(V(merged)$name %in% gordon)["TRUE"])
    df <- c(df, table(V(merged)$name %in% stukalov)["TRUE"])
    df <- c(df, gsize(merged))
    return(df)
}

plotHist <- function(value, title, phenotype, length, xmax, y1, y2) {
    dens_gwas <- hist(value, breaks = c(0:(max(value) + 1)), plot = FALSE, right = F)
    plot(dens_gwas, xlim = c(0, xmax), col = rgb(0.75, 0.75, 0.75, 1/2), border = NA, las = 1, xaxt = "n", freq = FALSE, xlab = "Number of viral targets", ylab = "Frequency", main = "", cex.sub = 0.5)
    mytitle <- paste0("COVID19 GWAS subnetwork\n(", phenotype, ")\nviral targets in ", title)
    mtext(side = 3, line = 1, cex = 1, mytitle)
    mtext(side = 3, line = 0.2, cex = 0.8, "subnetwork extracted from BioPlex3.0")
    axis(side = 1, at = seq(0, xmax, by = 5) + 0.5, labels = seq(0, xmax, by = 5))
    arrows(length + 0.5, y1, length + 0.5, 0, col = "#922687", lwd = 2, length = 0.1)
    text(median(value) + 4, max(dens_gwas$counts / 10000), paste0("median = ", median(value)), col = "grey", cex = 0.5)
    text(length - 2, y2, paste0("observed = ", length, "\np = ", table(value >= length)["TRUE"]/10000), cex = 0.4, pos = 4)
}
######
# load dataset
bioplex <- read.delim("../data/extended_table/BioPlex.3.0_edge.tsv", header = T)
husci <- read.csv("../data/HuSCI_node.csv", header = TRUE)
gwas <- read.xlsx("~/Documents/INET-work/virus_network/statistic_results/GWAS/COVID_GWAS_hit_inHUSCI_v2.xlsx")
gordon <- read.xlsx("../data/extended_table/Extended_Table_2_PPIs.xlsx", sheet = "Gordon")
stukalov <- read.xlsx("../data/extended_table/Extended_Table_2_PPIs.xlsx", sheet = "Stukalov")

######
# 1. BioPlex graph generation
bioplex_symbol <- bioplex[, c(5:6)]
bioplex_g_ori <- graph_from_data_frame(bioplex_symbol, directed = FALSE) # V:13957, E:118162
bioplex_g <- simplify(bioplex_g_ori, remove.loops = FALSE) # V:13957, E:118162

# GWAS list
ctcl <- gwas[, 1][gwas[, 6] == 1]
ctcl <- unique(ctcl[!is.na(ctcl)])
ctcl_bioplex <- ctcl[ctcl %in% V(bioplex_g)$name] # V:10
ctcl_1st <- combineNetwork(bioplex_g, ctcl_bioplex)

hosp <- gwas[, 1][gwas[, 7] == 1]
hosp <- unique(hosp[!is.na(hosp)])
hosp_bioplex <- hosp[hosp %in% V(bioplex_g)$name] # V:17
hosp_1st <- combineNetwork(bioplex_g, hosp_bioplex)

infct <- gwas[, 1][gwas[, 8] == 1]
infct <- unique(infct[!is.na(infct)])
infct_bioplex <- infct[infct %in% V(bioplex_g)$name] # V:10
infct_1st <- combineNetwork(bioplex_g, infct_bioplex)

gwas_bioplex <- gwas$All.LD[gwas$All.LD %in% V(bioplex_g)$name] # 24 of 42 candidates found in BioPlex3.0

# HuSCI, Gordon and Stukalov in BioPlex
husci_sym <- husci$node # V:171
husci_bioplex <- V(bioplex_g)$name[V(bioplex_g)$name %in% husci_sym] # HuSCI in BioPlex whole, V:132

gordon_sym <- unique(gordon$PreyGene) # V:384
gordon_bioplex <- V(bioplex_g)$name[V(bioplex_g)$name %in% gordon_sym] # V:346

stukalov_sym <- unique(stukalov$human) # V:876
stukalov_bioplex <- V(bioplex_g)$name[V(bioplex_g)$name %in% stukalov_sym] # V:723

######
# 2. interactor of GWAS hit
gwas_hit_1st <- make_ego_graph(bioplex_g, nodes = gwas_bioplex, order = 1, mode = "all") # 24 of 42 GWAS in BioPlex3.0

######
# 3. **rewiring analysis of HuRI**, to see if the HuSCI viral target is significant.
# subnetwork of GWAS hit from BioPlex
# inherit from above code
gwas_all_list_df <- lapply(gwas_hit_1st, as_data_frame)
gwas_all_df <- do.call(rbind, gwas_all_list_df)
gwas_all_g_merge <- graph_from_data_frame(gwas_all_df, directed = FALSE)
# to have interaction between 1st interactors
gwas_all_final <- simplify(induced_subgraph(bioplex_g, names(V(gwas_all_g_merge))), remove.loops = F)

# GWAS hit in HuSCI
gwas_all_husci <- V(gwas_all_final)$name[V(gwas_all_final)$name %in% husci_sym] # V:5
gwas_all_husci_length <- length(gwas_all_husci)
gwas_ctcl_husci <- V(ctcl_1st)$name[V(ctcl_1st)$name %in% husci_sym] # V:2
gwas_ctcl_husci_length <- length(gwas_ctcl_husci)
gwas_hosp_husci <- V(hosp_1st)$name[V(hosp_1st)$name %in% husci_sym] # V:2
gwas_hosp_husci_length <- length(gwas_hosp_husci)
gwas_infct_husci <- V(infct_1st)$name[V(infct_1st)$name %in% husci_sym] # V:4
gwas_infct_husci_length <- length(gwas_infct_husci)

# GWAS hit in Gordon
gwas_all_gordon <- V(gwas_all_final)$name[V(gwas_all_final)$name %in% gordon_sym] # V:10
gwas_all_gordon_length <- length(gwas_all_gordon)
gwas_ctcl_gordon <- V(ctcl_1st)$name[V(ctcl_1st)$name %in% gordon_sym]
gwas_ctcl_gordon_length <- length(gwas_ctcl_gordon)
gwas_hosp_gordon <- V(hosp_1st)$name[V(hosp_1st)$name %in% gordon_sym]
gwas_hosp_gordon_length <- length(gwas_hosp_gordon)
gwas_infct_gordon <- V(infct_1st)$name[V(infct_1st)$name %in% gordon_sym]
gwas_infct_gordon_length <- length(gwas_infct_gordon)

# GWAS hit in Stukalov
gwas_all_stukalov <- V(gwas_all_final)$name[V(gwas_all_final)$name %in% stukalov_sym] # V:19
gwas_all_stukalov_length <- length(gwas_all_stukalov)
gwas_ctcl_stukalov <- V(ctcl_1st)$name[V(ctcl_1st)$name %in% stukalov_sym]
gwas_ctcl_stukalov_length <- length(gwas_ctcl_stukalov)
gwas_hosp_stukalov <- V(hosp_1st)$name[V(hosp_1st)$name %in% stukalov_sym]
gwas_hosp_stukalov_length <- length(gwas_hosp_stukalov)
gwas_infct_stukalov <- V(infct_1st)$name[V(infct_1st)$name %in% stukalov_sym]
gwas_infct_stukalov_length <- length(gwas_infct_stukalov)

######
# permutation analysis
gwas_rand_r2 <- c()
gwas_rand_r2 <- c(gwas_rand_r2, mcreplicate(10000, bioplexRewireMulti(gwas_bioplex, ctcl_bioplex, hosp_bioplex, infct_bioplex, husci_sym, gordon_sym, stukalov_sym), mc.cores = detectCores()))
gwas_rand_r2[is.na(gwas_rand_r2)] <- 0
gwas_rand_df_r2 <- data.frame(matrix(gwas_rand_r2, ncol = 16, byrow = T))
names(gwas_rand_df_r2) <- c(
    "allGWAS_viral_target_inHuSCI",
    "allGWAS_viral_target_inGordon",
    "allGWAS_viral_target_inStukalov",
    "allGWAS_subnetworkSize",
    "ctclGWAS_viral_target_inHuSCI",
    "ctclGWAS_viral_target_inGordon",
    "ctclGWAS_viral_target_inStukalov",
    "ctclGWAS_subnetworkSize",
    "hospGWAS_viral_target_inHuSCI",
    "hospGWAS_viral_target_inGordon",
    "hospGWAS_viral_target_inStukalov",
    "hospGWAS_subnetworkSize",
    "infctGWAS_viral_target_inHuSCI",
    "infctGWAS_viral_target_inGordon",
    "infctGWAS_viral_target_inStukalov",
    "infctGWAS_subnetworkSize"
)

######
# plot
pdf(file = "~/Documents/INET-work/virus_network/figure_results/GWAS/GWASv2_3dataset_BioPlex3.pdf", width = 3, height = 3)
par(mgp = c(2, 0.7, 0), ps = 8)
# HuSCI in GWAS+1 all
plotHist(gwas_rand_df_r2$allGWAS_viral_target_inHuSCI, "HuSCI", "all 24 genes and 1st interactors", gwas_all_husci_length, 20, 0.03, 0.05)
# Gordon in GWAS+1 all
plotHist(gwas_rand_df_r2$allGWAS_viral_target_inGordon, "Gordon et al", "all 24 genes and 1st interactors", gwas_all_gordon_length, 25, 0.03, 0.05)
# Stukalov in GWAS+1 all
plotHist(gwas_rand_df_r2$allGWAS_viral_target_inStukalov, "Stukalov et al", "all 24 genes and 1st interactors", gwas_all_stukalov_length, 45, 0.03, 0.05)

# HuSCI in GWAS+1 critical illness
plotHist(gwas_rand_df_r2$ctclGWAS_viral_target_inHuSCI, "HuSCI", "critical illness, 10 genes and 1st interactors", gwas_ctcl_husci_length, 15, 0.03, 0.05)
# Gordon in GWAS+1 critical illness
plotHist(gwas_rand_df_r2$ctclGWAS_viral_target_inGordon, "Gordon et al", "critical illness, 10 genes and 1st interactors", gwas_ctcl_gordon_length, 20, 0.03, 0.05)
# Stukalov in GWAS+1 critical illness
plotHist(gwas_rand_df_r2$ctclGWAS_viral_target_inStukalov, "Stukalov et al", "critical illness, 10 genes and 1st interactors", gwas_ctcl_stukalov_length, 30, 0.03, 0.05)

# HuSCI in GWAS+1 hospitalization
plotHist(gwas_rand_df_r2$hospGWAS_viral_target_inHuSCI, "HuSCI", "hospitalization, 17 genes and 1st interactors", gwas_hosp_husci_length, 15, 0.03, 0.05)
# Gordon in GWAS+1 hospitalization
plotHist(gwas_rand_df_r2$hospGWAS_viral_target_inGordon, "Gordon et al", "hospitalization, 17 genes and 1st interactors", gwas_hosp_gordon_length, 20, 0.03, 0.05)
# Stukalov in GWAS+1 hospitalization
plotHist(gwas_rand_df_r2$hospGWAS_viral_target_inStukalov, "Stukalov et al", "hospitalization, 17 genes and 1st interactors", gwas_hosp_stukalov_length, 30, 0.03, 0.05)

# HuSCI in GWAS+1 reported infection
plotHist(gwas_rand_df_r2$infctGWAS_viral_target_inHuSCI, "HuSCI", "reported infection, 10 genes and 1st interactors", gwas_infct_husci_length, 15, 0.03, 0.05)
# Gordon in GWAS+1 reported infection
plotHist(gwas_rand_df_r2$infctGWAS_viral_target_inGordon, "Gordon et al", "reported infection, 10 genes and 1st interactors", gwas_infct_gordon_length, 20, 0.03, 0.05)
# Stukalov in GWAS+1 reported infection
plotHist(gwas_rand_df_r2$infctGWAS_viral_target_inStukalov, "Stukalov et al", "reported infection, 10 genes and 1st interactor", gwas_infct_stukalov_length, 30, 0.03, 0.05)
dev.off()

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

######
# save workarea data
save.image("~/Documents/INET-work/virus_network/statistic_results/3dataset_degree_BioPlex.RData")
