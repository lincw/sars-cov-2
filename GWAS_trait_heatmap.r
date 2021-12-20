# Heatmap for GWAS hit traits and community correlation
# Lin Chung-wen

######
# environment
google <- "/Volumes/GoogleDrive/My Drive/VirHostome_CW"
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
# COVID19
data_covid <- read.xlsx(file.path(google, "GWAS/data/magma_gtex_gwas_by_community_with_gwas_info_COVID.xlsx"), sheet = "BETA_SE")
data_covid_fdr <- read.xlsx(file.path(google, "GWAS/data/magma_gtex_gwas_by_community_with_gwas_info_COVID.xlsx"), sheet = "FDR")
data_covid_p <- read.xlsx(file.path(google, "GWAS/data/magma_gtex_gwas_by_community_with_gwas_info_COVID.xlsx"), sheet = "P")
rownames(data_covid) <- data_covid$VARIABLE
rownames(data_covid_fdr) <- data_covid_fdr$VARIABLE
rownames(data_covid_p) <- data_covid_p$VARIABLE

annotation <- read.xlsx(file.path(google, "communities/data/HuRI_GO_annotation.xlsx"), sheet = "go")
df <- data.frame(annotation[, c(2, 3)])
rownames(df) <- annotation$community
colnames(df) <- c("Shared functional term", "Membership overlap")

id_table <- read.xlsx(file.path(google, "communities/data/HuRI_GO_annotation.xlsx"), sheet = "id")

######
# data process with COVID, FDR < 0.05
fdr005_covid <- filterFDR(data_covid_fdr, 0.05)
beta005_covid <- data_covid[rownames(fdr005_covid), colnames(fdr005_covid)]
p005_covid <- data_covid_p[rownames(fdr005_covid), colnames(fdr005_covid)]

######
# plotting
paletteLength <- 50
myColor <- colorRampPalette(c("red", "white", "blue"))(paletteLength)

myBreaks_covid <- c(seq(min(beta005_covid), 0, length.out = ceiling(paletteLength / 2) + 1), seq(max(beta005_covid) / paletteLength, max(beta005_covid), length.out = floor(paletteLength / 2)))

col_fun <- colorRamp2(c(-2, 0, 4), c("blue", "white", "red"))
col_fun_p <- colorRamp2(c(0, 1), c("white", "red"))
col_fun_p10 <- colorRamp2(c(0, 6), c("white", "red"))

######
# plotting V1
phtmap_covid <- pheatmap::pheatmap(t(beta005_covid), silent = T)
col_dend_covid <- phtmap_covid[[2]]
row_dend_covid <- phtmap_covid[[1]]

# reorder rows, COVID19
row_dend_covid <- rotate(row_dend_covid, order = c("COVID19", "IBD_UKBS", "OST_UKBS", "BMIA", "NEUROT_UKB", "T2D_UKBS", "HRET", "RET", "ADPN", "PHF", "HEIGHT", "HIP", "FAT_UKB", "HC_UKBS", "HYPOTHY_UKBS", "SCZ_UKBS"))

# reorder columns, COVID19
col_re_d <- dist(beta005_covid)
col_re_hc <- hclust(col_re_d)
col_re_olo <- reorder(col_re_hc, col_re_d, method = "olo")
col_dend_covid_olo <- rotate(col_dend_covid, order = col_re_olo$labels[col_re_olo$order][c(1:8, 10, 9, 14, 13, 15:18, 11:12, 25:27, 31, 30, 28:29, 19, 21, 20, 22, 24, 23)])

pheatmap::pheatmap(t(beta005_covid),
    labels_col = id_table[match(rownames(beta005_covid), id_table[, 1]), 2],
    cluster_cols = as.hclust(col_dend_covid_olo), cluster_row = as.hclust(row_dend_covid),
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
    filename = "GWAS_trait_heatmap_olo_v3.pdf",
    width = 10, height = 4
)
