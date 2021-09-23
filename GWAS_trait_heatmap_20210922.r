# Heatmap for GWAS hit traits and community correlation
# Lin Chung-wen

######
# environment
google <- "/Volumes/GoogleDrive/My Drive/VirHostome_CW/GWAS"
stats <- "~/Documents/INET-work/virus_network/statistic_results/community"
######
# load packages
library(openxlsx)
library(pheatmap)
library(circlize)
library(seriation)
library(dendextend)

######
# function
filterFDR <- function(query, value) {
    row <- apply(query, 1, function(x) any(x < value))
    rowName <- rownames(query)[row]
    rowName <- rowName[!is.na(rowName)]
    col <- apply(query, 2, function(x) any(x < value))
    colName <- names(query)[col]
    colName <- colName[!is.na(colName)]
    community <- query$VARIABLE[row]
    community <- community[!is.na(community)]
    df <- query[rowName, colName]
    rownames(df) <- community
    return(df)
}
######
# load data
data <- read.xlsx(file.path(google, "/data/magma_gtex_gwas_by_community_with_gwas_info_20210921.xlsx"), sheet = "BETA_SE")
data_fdr <- read.xlsx(file.path(google, "/data/magma_gtex_gwas_by_community_with_gwas_info_20210921.xlsx"), sheet = "FDR")
data_p <- read.xlsx(file.path(google, "/data/magma_gtex_gwas_by_community_with_gwas_info_20210921.xlsx"), sheet = "P")
rownames(data) <- data$VARIABLE

annotation <- read.xlsx(file.path(stats, "HuRI_GO_annotation.xlsx"))
df <- data.frame(annotation[, c(2, 3)])
rownames(df) <- annotation$community
######
# data process, FDR < 0.05
fdr005 <- filterFDR(data_fdr, 0.05)
beta005 <- data[rownames(fdr005), names(fdr005)]

######
# plotting
paletteLength <- 50
myColor <- colorRampPalette(c("red", "white", "blue"))(paletteLength)
myBreaks <- c(seq(min(beta005), 0, length.out = ceiling(paletteLength/2) + 1), seq(max(beta005)/paletteLength, max(beta005), length.out = floor(paletteLength/2)))
col_fun <- colorRamp2(c(-2, 0, 4), c("blue", "white", "red"))
col_fun_p <- colorRamp2(c(0, 1), c("white", "red"))
col_fun_p10 <- colorRamp2(c(0, 6), c("white", "red"))

######
# plotting V1
phtmap <- pheatmap::pheatmap(t(beta005))
col_dend <- phtmap[[2]]
row_dend <- phtmap[[1]]
# reorder columns
col_dend <- rotate(col_dend, order = c("415", "3941", "1652", "377", "571", "2545", "692", "926", "2398", "525", "1882", "1833", "731", "316", "4227", "3682", "2769", "7", "2831", "4160", "1353"))
# reorder rows
row_dend <- rotate(row_dend, order = c("IBD_UKBS", "OST_UKBS", "BMIA", "NEUROT_UKB", "T2D_UKBS", "HRET", "RET", "ADPN", "PHF", "HEIGHT", "HIP", "FAT_UKB", "HC_UKBS", "HYPOTHY_UKBS", "SCZ_UKBS"))
pdf(file.path(google, "figures/gwas_beta_fdr005.pdf"), height = 4, width = 6)
pheatmap::pheatmap(t(beta005),
    cluster_cols = as.hclust(col_dend), cluster_row = as.hclust(row_dend),
    color = colorRampPalette(c("blue", "white", "red"))(50), breaks = myBreaks, border_color = "white",
    display_numbers = ifelse(t(fdr005) < 0.05, "*", ""),
    fontsize_number = 10, fontsize_col = 10, fontsize_row = 10, number_color = "black",
    angle_col = 90,
    treeheight_col = 20,
    treeheight_row = 20,
    legend_breaks = c(-3, -2, -1, 0, 1, 2, 3, 4, 4.6),
    legend_labels = c("-3", "-2", "-1", "0", "1", "2", "3", "4", "BETA/SE\n"),
    legend = TRUE,
    main = "Community GWAS hits (FDR < 0.05)",
    annotation = df, annotation_legend = FALSE)
dev.off()
