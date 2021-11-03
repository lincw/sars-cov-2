# Heatmap for GWAS hit traits and community correlation
# Lin Chung-wen

######
# environment
google <- "/Volumes/GoogleDrive/My Drive/VirHostome_CW/GWAS"
stats <- "~/Documents/INET-work/virus_network/statistic_results"
######
# load packages
library(openxlsx)
library(pheatmap)
library(circlize)
library(seriation)
library(dendextend)
library(seriation)

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

## COVID19
data_covid <- read.xlsx(file.path(google, "/data/magma_gtex_gwas_by_community_with_gwas_info_COVID.xlsx"), sheet = "BETA_SE")
data_covid_fdr <- read.xlsx(file.path(google, "/data/magma_gtex_gwas_by_community_with_gwas_info_COVID.xlsx"), sheet = "FDR")
data_covid_p <- read.xlsx(file.path(google, "/data/magma_gtex_gwas_by_community_with_gwas_info_COVID.xlsx"), sheet = "P")
rownames(data_covid) <- data_covid$VARIABLE
rownames(data_covid_fdr) <- data_covid_fdr$VARIABLE
rownames(data_covid_p) <- data_covid_p$VARIABLE

annotation <- read.xlsx(file.path(stats, "community/HuRI_GO_annotation.xlsx"))
df <- data.frame(annotation[, c(2, 3)])
rownames(df) <- annotation$community
colnames(df) <- c("Shared functional term", "Membership similarity")
######
# data process, FDR < 0.05
fdr005 <- filterFDR(data_fdr, 0.05)
beta005 <- data[rownames(fdr005), names(fdr005)]

## data process with COVID, FDR < 0.05
fdr005_covid <- filterFDR(data_covid_fdr, 0.05)
beta005_covid <- data_covid[rownames(fdr005_covid), colnames(fdr005_covid)]
p005_covid <- data_covid_p[rownames(fdr005_covid), colnames(fdr005_covid)]

######
# plotting
paletteLength <- 50
myColor <- colorRampPalette(c("red", "white", "blue"))(paletteLength)
myBreaks <- c(seq(min(beta005), 0, length.out = ceiling(paletteLength / 2) + 1), seq(max(beta005) / paletteLength, max(beta005), length.out = floor(paletteLength / 2)))

myBreaks_covid <- c(seq(min(beta005_covid), 0, length.out = ceiling(paletteLength / 2) + 1), seq(max(beta005_covid) / paletteLength, max(beta005_covid), length.out = floor(paletteLength / 2)))

col_fun <- colorRamp2(c(-2, 0, 4), c("blue", "white", "red"))
col_fun_p <- colorRamp2(c(0, 1), c("white", "red"))
col_fun_p10 <- colorRamp2(c(0, 6), c("white", "red"))

######
# plotting V1
phtmap <- pheatmap::pheatmap(t(beta005), silent = T)
col_dend <- phtmap[[2]]
row_dend <- phtmap[[1]]

phtmap_covid <- pheatmap::pheatmap(t(beta005_covid), silent = T)
col_dend_covid <- phtmap_covid[[2]]
row_dend_covid <- phtmap_covid[[1]]

# reorder columns
col_dend <- rotate(col_dend, order = c("415", "3941", "1652", "377", "571", "2545", "692", "926", "2398", "525", "1882", "4227", "316", "731", "1833", "3682", "2769", "7", "2831", "4160", "1353"))
# reorder rows
row_dend <- rotate(row_dend, order = c("IBD_UKBS", "OST_UKBS", "BMIA", "NEUROT_UKB", "T2D_UKBS", "HRET", "RET", "ADPN", "PHF", "HEIGHT", "HIP", "FAT_UKB", "HC_UKBS", "HYPOTHY_UKBS", "SCZ_UKBS"))

# reorder columns, COVID19
col_dend_covid_v2 <- rotate(col_dend_covid, order = c("2563", "3442", "2239", "731", "1833", "4215", "479", "1900", "1046", "592", "64", "895", "1353", "4160", "2831", "7", "2769", "3682", "316", "4227", "571", "415", "3941", "1652", "377", "2398", "2545", "692", "926", "525", "1882"))
col_dend_covid_v3 <- rotate(col_dend_covid, order = c("2769", "3682", "7", "2831", "4160", "731", "2239", "2563", "3442", "895", "64", "1900", "1046", "592", "4215", "479", "1833", "1353", "2398", "2545", "692", "926", "525", "1882", "316", "4227", "571", "415", "3941", "1652", "377"))
col_dend_covid_v4 <- rotate(col_dend_covid, order = c("731", "2239", "2563", "3442", "895", "64", "1900", "1046", "592", "4215", "479", "1833", "2769", "3682", "7", "2831", "4160", "1353", "2398", "2545", "692", "926", "525", "1882", "316", "4227", "571", "415", "3941", "1652", "377"))
col_dend_covid_v6 <- rotate(col_dend_covid, order = c("2563", "3442", "2239", "731", "1833", "4215", "479", "1900", "1046", "592", "64", "895", "1353", "4160", "2831", "7", "2769", "3682", "316", "4227", "571", "415", "3941", "1652", "377", "2398", "2545", "692", "926", "525", "1882"))
# reorder rows, COVID19
row_dend_covid <- rotate(row_dend_covid, order = c("COVID19", "IBD_UKBS", "OST_UKBS", "BMIA", "NEUROT_UKB", "T2D_UKBS", "HRET", "RET", "ADPN", "PHF", "HEIGHT", "HIP", "FAT_UKB", "HC_UKBS", "HYPOTHY_UKBS", "SCZ_UKBS"))

col_re_d <- dist(beta005_covid)
col_re_hc <- hclust(col_re_d)
col_re_GW <- reorder(col_re_hc, col_re_d, method = "GW")
col_re_olo <- reorder(col_re_hc, col_re_d, method = "olo")
col_dend_covid_GW <- rotate(col_dend_covid, order = col_re_GW$labels[col_re_GW$order])
col_dend_covid_olo <- rotate(col_dend_covid, order = col_re_olo$labels[col_re_olo$order])
col_dend_covid_olo2 <- rotate(col_dend_covid, order = col_re_olo$labels[col_re_olo$order][c(1:18, 25:31, 19:24)])

pheatmap::pheatmap(t(beta005_covid),
    cluster_cols = as.hclust(col_dend_covid_olo2), cluster_row = as.hclust(row_dend_covid),
    color = colorRampPalette(c("blue", "white", "red"))(50), breaks = myBreaks_covid, border_color = "white",
    display_numbers = ifelse(t(fdr005_covid) < 0.05, "*", ""),
    fontsize_number = 10, fontsize_col = 10, fontsize_row = 10, number_color = "black",
    angle_col = 90,
    treeheight_col = 20,
    treeheight_row = 20,
    legend_breaks = c(-3, -2, -1, 0, 1, 2, 3, 4, 4.6),
    legend_labels = c("-3", "-2", "-1", "0", "1", "2", "3", "4", "BETA/SE\n"),
    legend = TRUE,
    main = "Community GWAS hits (FDR < 0.05)",
    annotation = df, annotation_legend = TRUE,
    filename = file.path(stats, "GWAS/GWAS_trait_heatmap_olo_v2.pdf"),
    width = 10, height = 4
)

######
# save raw value
write.xlsx(list(BETA_SE = beta005_covid, P = p005_covid, FDR = fdr005_covid), file.path(stats, "GWAS/GWAS_trait_stats.xlsx"), row.names = T, overwrite = T)

save.image(file.path(stats, "GWAS/GWAS_trait.RData"))
