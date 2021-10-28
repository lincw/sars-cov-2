# HuRI communities jaccard similarity analysis. The communities were used in GWAS trait analysis
# Lin Chung-wen
# 10.06.2021
# 23.09.2021

#######
# load package
library(openxlsx)
library(linkcomm)
library(pheatmap)
# library(ComplexHeatmap)
library(circlize) # color used for ComplexHeatmap
library(seriation)
library(dendextend)
# library(ctc) # output hclust data as .cdt for java treeview
######
# plot heatmap
pheatmapPlotBasic <- function(data, main) {
    pheatmap(t(data),
    color = colorRampPalette(c("blue", "white", "red"))(50), breaks = myBreaks,
        fontsize_number = 10, fontsize_col = 8, fontsize_row = 8, number_color = "black",
    legend_breaks = c(-3, -2, -1, 0, 1, 2, 3, 4, 4.6),
    legend_labels = c("-3", "-2", "-1", "0", "1", "2", "3", "4", "BETA/SE\n"),
    legend = TRUE,
    main = main)
}
######
# load OCG data
load("~/Documents/INET-work/virus_network/statistic_results/community/HuRI_ocg.RData")

######
# function, define Jaccard Similarity function
jaccard <- function(a, b) {
    intersection <- length(intersect(a, b))
    union <- length(a) + length(b) - intersection
    return (intersection/union)
}

######
# The communities used in GWAS trait analysis
community <- c(3682, 2769, 731, 316, 4227, 1353, 4160, 2831, 1833, 2545, 692, 926, 2398, 525, 1882, 1652, 3941, 415, 377, 571, 7, 2563, 3442, 2239, 4215, 479, 1900, 1046, 592, 64, 895)

comm_member <- list()
for (i in 1:length(community)) {
    comm_member[[i]] <- getNodesIn(huri_ocg, clusterids = community[i])
}
write.xlsx(unlist(lapply(comm_member, function(x) paste0(x, collapse = ","))), file = "~/Documents/INET-work/virus_network/statistic_results/community/gwas_trait_members.xlsx")

######
# Jaccard similarity calculation
jac_result <- c()
for (i in 1:length(community)) {
    for (j in 1:length(community)) {
        jac_result <- c(jac_result, jaccard(comm_member[[i]], comm_member[[j]]))
    }
}

jac_result_df <- matrix(jac_result, ncol = 31)
rownames(jac_result_df) <- community
colnames(jac_result_df) <- community

write.xlsx(as.data.frame(jac_result_df), file = "~/Documents/INET-work/virus_network/statistic_results/GWAS/GWAS_trait_jaccard.xlsx", row.names = T, overwrite = T)

paletteLength <- 50
myColor_jac <- colorRampPalette(c("white", "red"))(paletteLength)
myBreaks_jac <- c(seq(min(jac_result_df), 0, length.out = ceiling(paletteLength/2) + 1), seq(max(jac_result_df)/paletteLength, max(jac_result_df), length.out = floor(paletteLength/2)))

pheatmap::pheatmap(jac_result_df, cutree_rows = 7, cutree_cols = 7, main = "Jaccard similarity of community membership",
    color = myColor_jac,
    # breaks = myBreaks_jac,
    legend_breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
    legend_labels = c("0", "0.2", "0.4", "0.6", "0.8", "score\n"),
    legend = TRUE,
    filename = "~/Documents/INET-work/virus_network/figure_results/gwas_trait_jaccard_20211028.pdf",
    width = 5, height = 5)

######
# extract heatmap from Matthias
# custom script for GWAS hits from Matthias
file <- "~/Documents/INET-work/virus_network/statistic_results/GWAS_Mattias/GWAS_Mattias/trait_matrix.xlsx"
data <- read.xlsx(file, sheet = "BETA_SE")
data_fdr <- read.xlsx(file, sheet = "FDR")
data_p <- read.xlsx(file, sheet = "pvalue")

plot_data <- data[, c(2:22)]
plot_fdr <- data_fdr[, c(2:22)]
plot_p <- data_p[, c(2:22)]

rownames(plot_fdr) <- data_fdr[, 1]
rownames(plot_data) <- data[, 1]
rownames(plot_p) <- data_p[, 1]

paletteLength <- 50
myColor <- colorRampPalette(c("red", "white", "blue"))(paletteLength)
myBreaks <- c(seq(min(plot_data), 0, length.out = ceiling(paletteLength/2) + 1), seq(max(plot_data)/paletteLength, max(plot_data), length.out = floor(paletteLength/2)))

pdf("~/Documents/INET-work/virus_network/figure_results/gwas_trait_jaccard_BETA_0.05.pdf", width = 6, height = 5)
# group by @Pascal, fdr < 0.05, significant community only
pheatmap(t(plot_data)[c(9, 2, 5, 17, 21, 7, 10, 11, 16, 20, 12, 13, 18, 1, 8, 14, 3, 4, 6, 15, 19), c(5, 6, 8, 16, 27, 3, 17, 21, 26, 13, 28, 7, 9, 12, 18:20, 29, 10, 22)],
    cluster_cols = F, cluster_rows = F,
    color = colorRampPalette(c("blue", "white", "red"))(50), breaks = myBreaks, display_numbers = t(plot_fdr_005)[c(9, 2, 5, 17, 21, 7, 10, 11, 16, 20, 12, 13, 18, 1, 8, 14, 3, 4, 6, 15, 19), c(5, 6, 8, 16, 27, 3, 17, 21, 26, 13, 28, 7, 9, 12, 18:20, 29, 10, 22)],
        fontsize_number = 10, fontsize_col = 8, fontsize_row = 8, number_color = "black",
    legend_breaks = c(-3, -2, -1, 0, 1, 2, 3, 4, 4.6),
    legend_labels = c("-3", "-2", "-1", "0", "1", "2", "3", "4", "BETA/SE\n"),
    legend = TRUE,
    gaps_col = c(5, 9, 11, 17, 18),
    gaps_row = c(1, 3, 5, 10, 11, 13, 15, 16),
    main = "Community GWAS hits (FDR < 0.05)\n(all sig & trait grouping)")
# group by Pascal, fdr < 0.05
pheatmap(t(plot_data)[c(9, 2, 5, 17, 21, 7, 10, 11, 16, 20, 12, 13, 18, 1, 8, 14, 3, 4, 6, 15, 19), c(1, 24, 25, 5, 6, 8, 16, 27, 3, 17, 21, 26, 13, 28, 2, 14, 15, 23, 7, 9, 12, 18:20, 11, 29, 4, 10, 22)],
    cluster_cols = F, cluster_rows = F,
    color = colorRampPalette(c("blue", "white", "red"))(50), breaks = myBreaks, display_numbers = t(plot_fdr_005)[c(9, 2, 5, 17, 21, 7, 10, 11, 16, 20, 12, 13, 18, 1, 8, 14, 3, 4, 6, 15, 19), c(1, 24, 25, 5, 6, 8, 16, 27, 3, 17, 21, 26, 13, 28, 2, 14, 15, 23, 7, 9, 12, 18:20, 11, 29, 4, 10, 22)],
        fontsize_number = 10, fontsize_col = 8, fontsize_row = 8, number_color = "black",
    legend_breaks = c(-3, -2, -1, 0, 1, 2, 3, 4, 4.6),
    legend_labels = c("-3", "-2", "-1", "0", "1", "2", "3", "4", "BETA/SE\n"),
    legend = TRUE,
    gaps_col = c(3, 8, 12, 14, 18, 24, 26),
    gaps_row = c(1, 3, 5, 10, 11, 13, 15, 16),
    main = "Community GWAS hits (FDR < 0.05)\n(trait grouping)")
# original Matthias plot
pheatmap(t(plot_data)[, c(1, 24, 25, 5, 6, 8, 16, 27, 3, 17, 21, 26, 13, 28, 2, 14, 15, 23, 7, 9, 12, 18:20, 11, 29, 4, 10, 22)], cutree_cols = 8, cutree_rows = 3,
    cluster_cols = F,
    color = colorRampPalette(c("blue", "white", "red"))(50), breaks = myBreaks, display_numbers = t(plot_fdr)[, c(1, 24, 25, 5, 6, 8, 16, 27, 3, 17, 21, 26, 13, 28, 2, 14, 15, 23, 7, 9, 12, 18:20, 11, 29, 4, 10, 22)],
        fontsize_number = 10, fontsize_col = 8, fontsize_row = 8, number_color = "black",
    legend_breaks = c(-3, -2, -1, 0, 1, 2, 3, 4, 4.6),
    legend_labels = c("-3", "-2", "-1", "0", "1", "2", "3", "4", "BETA/SE\n"),
    legend = TRUE,
    gaps_col = c(3, 8, 12, 14, 18, 24, 26),
    main = "Community GWAS hits (FDR < 0.1)")
# original Matthias plot, fdr < 0.05
pheatmap(t(plot_data)[, c(1, 24, 25, 5, 6, 8, 16, 27, 3, 17, 21, 26, 13, 28, 2, 14, 15, 23, 7, 9, 12, 18:20, 11, 29, 4, 10, 22)], cutree_cols = 8, cutree_rows = 3,
    cluster_cols = F,
    color = colorRampPalette(c("blue", "white", "red"))(50), breaks = myBreaks, display_numbers = t(plot_fdr_005)[, c(1, 24, 25, 5, 6, 8, 16, 27, 3, 17, 21, 26, 13, 28, 2, 14, 15, 23, 7, 9, 12, 18:20, 11, 29, 4, 10, 22)],
        fontsize_number = 10, fontsize_col = 8, fontsize_row = 8, number_color = "black",
    legend_breaks = c(-3, -2, -1, 0, 1, 2, 3, 4, 4.6),
    legend_labels = c("-3", "-2", "-1", "0", "1", "2", "3", "4", "BETA/SE\n"),
    legend = TRUE,
    gaps_col = c(3, 8, 12, 14, 18, 24, 26),
    main = "Community GWAS hits (FDR < 0.05)")
dev.off()

######
# play around with clustering
col_fun <- colorRamp2(c(-2, 0, 4), c("blue", "white", "red"))
col_fun_p <- colorRamp2(c(0, 1), c("white", "red"))
col_fun_p10 <- colorRamp2(c(0, 6), c("white", "red"))
# 1. cluster with BETA/SE score
# 2. collapse everything between 1 and -1 as 0
plot_data_2 <- plot_data
plot_data_2[abs(plot_data_2) < 1]  <- 0
# 1.1 & 2.1 flip 90 degree of the table
# 1.3 & 2.3 from 1, with cell gaps
# 3. cluster based on p-value
plot_p10 <- -log(plot_p, 10)
plot_fdr10 <- -log(plot_fdr, 10)
# 4. filter all matrix by FDR < 0.05
row_dim <- c(3, 5:13, 16:22, 26:29)
col_dim <- c(1:14, 16)
# 5. binary
plot_p_binary0.1 <- ifelse(plot_fdr < 0.1, 1, 0)
source("~/Documents/INET-work/virus_network/src/gwasHit_cluster_plot.r")

# reorder dendrogram
# custom trait order (FDR might not < 0.05?)
row_dim2 <- c(16, 9, 11, 4, 10,  2, 5, 7, 12, 8, 13, 14, 6, 3, 1)
non_community_list <- c(891, 3205, 8, 3564, 364, 145, 1515, 3035, 1627)
source("~/Documents/INET-work/virus_network/src/gwasHit_cluster_plot_reorder.r")

gwas_comm <- read.xlsx("~/workplace/GWAS_list.xlsx", sheet = "community")
gwas_trait <- read.xlsx("~/workplace/GWAS_list.xlsx", sheet = "trait")
comm_position <- gwas_comm$index[!gwas_comm$community %in% non_community_list]

comm_cluster <- data.frame(
    community_cluster = c("A", "A",
        "A", "A",
        "A", "D",
        "D", "D",
        "D", "D",
        "D",
        "B",
        "B", "B",
        "B", "E", "E",
        "E", "C", "C"
    )
)
rownames(comm_cluster) <- c("415", "3941",
        "1652", "377",
        "571", "2545",
        "2398", "525",
        "926", "1882",
        "692",
        "1833",
        "3682", "2769",
        "316", "4227", "731",
        "2831", "4160", "1353")
trait_cluster <- data.frame(
    trait_cluster = c("I",
        "I",
        "I", "II",
        "III",
        "IV", "IV",
        "IV", "II",
        "IV", "IV",
        "IV",
        "III",
        "II", "II"
    )
)
rownames(trait_cluster) <- c("IBD_UKBS",
        "T2D_UKBS",
        "OST_UKBS", "BMIA",
        "NEUROT_UKB",
        "HRET", "RET",
        "ADPN", "PHF",
        "HC_UKBS", "FAT_UKB",
        "HYPOTHY_UKBS",
        "SCZ_UKBS",
        "HIP",
        "HEIGHT"
)
cluster_color <- list(trait_cluster = c(
    "I" = "#E6194B", "II" = "#F58231", "III" = "#FFE119", "IV" = "#42D4F4"),
    community_cluster = c(
    A = "#800000", B = "#9A6324", C = "#808000", D = "#469990", E = "#000075"))
