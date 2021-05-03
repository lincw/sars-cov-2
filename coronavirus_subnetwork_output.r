library(pheatmap)
library(openxlsx)
library(igraph)
library(reshape2)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
###########################################################
# set work directory
setwd("/tmp")

###########################################################
# custom function
grab1 <- function(network, id) {
    make_ego_graph(network, order = 1, nodes = id, mode = "all", mindist = 0)
}

###########################################################
# Human genome, GRCh38.p13
human <- read.csv("~/workplace/Human_gene_GRCh38.p13.csv", header = T)
names(human) <- c("ensemblID", "description", "symbol", "entrezID", "family")

###########################################################
# Human reference interactome
node <- read.csv("~/Documents/INET-work/references/HuRI_binaryPPI/HuRI_node.csv", header = T, as.is = T)
edge <- read.csv("~/Documents/INET-work/references/HuRI_binaryPPI/HuRI_edge.csv", header = T, as.is = T)
net <- graph_from_data_frame(d = edge, vertices = node, directed = F)

###########################################################
# load all coronavirus candidators from INET screening
corona_interactor <- read.xlsx("~/Documents/INET-work/virus_network/Y2H_screening/INET_7coronavirus_node.xlsx")

###########################################################
# reveal subnetwork of SARS-CoV-2 interactor
sars2_interactor <- read.xlsx("~/Documents/INET-work/virus_network/Y2H_screening/INET_R1.xlsx")
inet_s2r1 <- unique(sars2_interactor[sars2_interactor$lab == "INET_R1", c("target_ensemblID", "target")])
for (i in 1:length(inet_s2r1$target_ensemblID)) {
    if(inet_s2r1$target_ensemblID[i] %in% V(net)$name) {
        id1 <- grab1(net, inet_s2r1$target_ensemblID[i])
        V(id1[[1]])$color <- "orange"
        V(id1[[1]])[inet_s2r1$target_ensemblID[i]]$color <- "pink"
        pdf(paste0(inet_s2r1$target_ensemblID[i], ".pdf"), width = 10, height = 10)
        plot(id1[[1]], vertex.label = V(id1[[1]])$symbol)
        text(x = 0.9, y = 0.9, paste0("degree = ", rev(sort(degree(id1[[1]])))[1]))
        dev.off()
    } else {
        pdf(paste0(inet_s2r1$target_ensemblID[i], ".pdf"), width = 10, height = 10)
        par(mar = c(0, 0, 0, 0))
        plot(c(0, 1), c(0, 1), ann = F, bty = "n", type = "n", xaxt = "n", yaxt = "n")
        text(x = 0.5, y = 0.5, paste0(inet_s2r1$target_ensemblID[i], " was not found in the Human Reference Interactome"))
        dev.off()
    }
}

## 24.09.2020 re-defined INET SARS-CoV-2 interactome
sars2_re <- read.xlsx("~/Documents/INET-work/virus_network/Y2H_screening/20200923_verified/INET_All_Interactions_200924.xlsx") # update 0924
names(sars2_re) <- c("symbol", "viralProtein", "sumABS", "description", "growth")
inet_re <- unique(human[human$symbol %in% sars2_re$"symbol", c("ensemblID", "symbol")])
mainDir <- getwd()
subDir <- "sars2_re"
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
setwd(file.path(mainDir, subDir))
for (i in 1:length(inet_re$ensemblID)) {
    if(inet_re$ensemblID[i] %in% V(net)$name) {
        id1 <- grab1(net, inet_re$ensemblID[i])
        V(id1[[1]])$color <- "orange"
        V(id1[[1]])[inet_re$ensemblID[i]]$color <- "pink"
        pdf(paste0(inet_re$ensemblID[i], ".pdf"), width = 10, height = 10)
        plot(id1[[1]], vertex.label = V(id1[[1]])$symbol, vertex.label.degree = -pi/2, vertex.label.cex = .6)
        # text(x = 0.9, y = 0.9, paste0("degree = ", rev(sort(degree(id1[[1]])))[1]))
        dev.off()
    } else {
        pdf(paste0(inet_re$ensemblID[i], ".pdf"), width = 10, height = 10)
        par(mar = c(0, 0, 0, 0))
        plot(c(0, 1), c(0, 1), ann = F, bty = "n", type = "n", xaxt = "n", yaxt = "n")
        text(x = 0.5, y = 0.5, paste0(inet_re$ensemblID[i], " was not found in the Human Reference Interactome"))
        dev.off()
    }
}

###########################################################
# integrate with other resource
## 1. RNA-seq and proteomoe of A549 cell under SARS-CoV-2 treatment
## ref: https://doi.org/10.1101/2020.06.17.156455
## RNA
sars2_rna <- read.xlsx("~/Documents/INET-work/virus_network/references/SARS-CoV2-4omics/STable4.xlsx", sheet = 2, na.strings = "", startRow = 3, colNames = FALSE)
sars2_rna_tmp <- sars2_rna[, c(2, 7, 9, 10, 15, 17, 18, 23, 25, 26, 31, 33, 34, 39, 41, 42, 47, 49, 50)]
names(sars2_rna_tmp) <- c("ensemblID", "FC_3h", "pvalue_3h", "adjP_3h", "FC_6h", "pvalue_6h", "adjP_6h", "FC_12h", "pvalue_12h", "adjP_12h", "FC_18h", "pvalue_18h", "adjP_18h", "FC_24h", "pvalue_24h", "adjP_24h", "FC_30h", "pvalue_30h", "adjP_30h")
inet_s2r1_rna <- sars2_rna_tmp[sars2_rna_tmp$ensemblID %in% inet_s2r1$target_ensemblID, ]
row.names(inet_s2r1_rna) <- inet_s2r1[match(inet_s2r1_rna$ensemblID, inet_s2r1$target_ensemblID), "target"]
inet_cor_rna <- sars2_rna_tmp[sars2_rna_tmp$ensemblID %in% corona_interactor$ensemblID, ]
row.names(inet_cor_rna) <- corona_interactor[match(inet_cor_rna$ensemblID, corona_interactor$ensemblID), "Symbol"]

## 24.09.2020 re-defined INET SARS-CoV-2 interactome
### RNA
inet_re_rna <- sars2_rna_tmp[sars2_rna_tmp$ensemblID %in% inet_re$ensemblID, ]
row.names(inet_re_rna) <- inet_re[match(inet_re_rna$ensemblID, inet_re$ensemblID), "symbol"]

## protein
sars2_prot <- read.xlsx("~/Documents/INET-work/virus_network/references/SARS-CoV2-4omics/STable5.xlsx", sheet = 2, na.strings = "", startRow = 3, colNames = FALSE)
sars2_prot_tmp <- sars2_prot[, c(2, 4, 9, 7, 8, 13, 11, 12)]
names(sars2_prot_tmp) <- c("geneName", "description", "FC_6h", "pvalue_6h", "adjP_6h", "FC_24h", "pvalue_24h", "adjP_24h")

inet_s2r1_prot <- sars2_prot_tmp[sars2_prot_tmp$geneName %in% row.names(inet_s2r1_rna), ]
row.names(inet_s2r1_prot) <- inet_s2r1_prot$geneName
inet_cor_prot <- sars2_prot_tmp[sars2_prot_tmp$ensemblID %in% corona_interactor$ensemblID, ]
row.names(inet_cor_prot) <- corona_interactor[match(inet_cor_prot$ensemblID, corona_interactor$ensemblID), "Symbol"]

## 24.09.2020 re-defined INET SARS-CoV-2 interactome
### protein
inet_re_prot <- sars2_prot_tmp[sars2_prot_tmp$geneName %in% row.names(inet_re_rna), ]
row.names(inet_re_prot) <- inet_re_prot$geneName

## plot
pdf("s2r1.pdf")
pheatmap(inet_s2r1_rna[, c(2, 5, 8, 11, 14, 17)])
pheatmap(inet_s2r1_prot[, c(3, 6)])
dev.off()

## 24.09.2020 re-defined INET SARS-CoV-2 interactome
pdf("inet_re.pdf")
pheatmap(inet_s2r1_rna[, c(2, 5, 8, 11, 14, 17)])
pheatmap(inet_s2r1_prot[, c(3, 6)])
dev.off()

###########################################################
# integrate with other resource
## 2. The lung cells from human COVID19 confirmed case
## ref: https://doi.org/10.1016/j.cell.2020.04.035
## !!!GRCh37.p13 was applied for assembly, not GRCh38.p13!!!
cell <- read.xlsx("~/Documents/INET-work/virus_network/references/SARS-CoV-2_ACE2_TMPRSS2_expression_Ziegler/1-s2.0-S0092867420305006-mmc3.xlsx", sheet = 1, start = 6)
grch37 <- read.csv("~/workplace/grch37.csv", header = F)
cell$ensemblID <- grch37[match(cell$Gene, grch37$V1), "V2"] # match ensemblID from GRCh37 version
cell$grch38 <- human[match(cell$ensemblID, human$ensemblID), "symbol"] # match symbol from GRCh38 version
inet_s2r1_cell <- cell[cell$grch38 %in% inet_s2r1$target, ]
names(inet_s2r1_cell) <- c("gene", "cell", "pvalue", "adjP", "FC", "expresPercent", "expressPercentRest")
inet_s2r1_cell_sub <- inet_s2r1_cell[, c(1, 2, 6)]
inet_s2r1_cell_melt <- melt(inet_s2r1_cell_sub)
inet_s2r1_cell_melt$cell <- factor(inet_s2r1_cell$cell, levels = c("Neutrophil", "Monocyte1", "Monocyte2", "Macrophage1", "Macrophage3", "Proliferating", "Fibroblast2", "Endothelial/Lymphatics", "Ciliated Cells", "Type I Pneumocytes", "Type II Pneumocytes"))
pdf("INET_SARS2_R1screening.pdf", width = 15, height = 10)
ggplot(inet_s2r1_cell_melt, aes(x = gene, y = cell)) +
    geom_point(aes(size = value * 100)) +
    labs(x = "interactor", y = "cell") +
    guides(size = guide_legend(title = "% Expression")) +
    theme_bw(base_size = 20) +
    theme(axis.text.x = element_text(angle = 45, color = "black", hjust = 1, vjust = 1),
        axis.text.y = element_text(color = "black")) +
    theme(plot.margin = unit(c(1, 4, 1, 4), "cm"))

inet_cor_cell <- cell[cell$Gene %in% corona_interactor$Symbol, ]
names(inet_cor_cell) <- c("gene", "cell", "pvalue", "adjP", "FC", "expresPercent", "expressPercentRest")
inet_cor_cell_sub <- inet_cor_cell[, c(1, 2, 6)]
inet_cor_cell_melt <- melt(inet_cor_cell_sub)
inet_cor_cell_melt$cell <- factor(inet_cor_cell$cell, levels = c("Mast", "Neutrophil", "Monocyte1", "Monocyte2", "Macrophage1", "Macrophage2", "Macrophage3", "Proliferating", "Fibroblast1", "Fibroblast2", "Endothelial/Lymphatics", "T", "Ciliated Cells", "Type I Pneumocytes", "Type II Pneumocytes"))
ggplot(inet_cor_cell_melt, aes(x = gene, y = cell)) +
    geom_point(aes(size = value * 100)) +
    labs(x = "interactor", y = "cell") +
    guides(size = guide_legend(title = "% Expression")) +
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(angle = 90, color = "black", hjust = 1, vjust = .5),
        axis.text.y = element_text(color = "black")) +
    theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
dev.off()

## 24.09.2020 re-defined INET SARS-CoV-2 interactome
inet_re_cell <- cell[cell$grch38 %in% inet_re$symbol, ]
names(inet_re_cell) <- c("gene", "cell", "pvalue", "adjP", "FC", "expresPercent", "expressPercentRest")
inet_re_cell_sub <- inet_re_cell[, c(1, 2, 6)]
inet_re_cell_melt <- melt(inet_re_cell_sub)
inet_re_cell_melt$cell <- factor(inet_re_cell$cell, levels = names(table(inet_re_cell_melt$cell)[c(2:13, 1, 14, 15)]))
pdf("INET_RE_screening.pdf", width = 20, height = 10)
ggplot(inet_re_cell_melt, aes(x = gene, y = cell)) +
    geom_point(aes(size = value * 100)) +
    labs(x = "interactor", y = "cell") +
    guides(size = guide_legend(title = "% Expression")) +
    theme_bw(base_size = 28) +
    theme(axis.text.x = element_text(angle = 45, color = "black", hjust = 1, vjust = 1),
        axis.text.y = element_text(color = "black")) +
    theme(plot.margin = unit(c(1, 4, 1, 4), "cm"))
dev.off()

## check the interactors of inet_s2r1 also expressed in SARS2 infected cells or not
### KRT8, ENSG00000170421
### USO1, ENSG00000138768
krt8 <- grab1(net, "ENSG00000170421")
krt8_interactors <- cell[cell$Gene %in% node[node$ensemblID %in% as_ids(V(krt8[[1]])), ]$symbol, ]
names(krt8_interactors) <- c("gene", "cell", "pvalue", "adjP", "FC", "expresPercent", "expressPercentRest")
krt8_interactors_sub <- krt8_interactors[, c(1, 2, 6)]
krt8_interactors_melt <- melt(krt8_interactors_sub)
krt8_interactors_melt$cell <- factor(krt8_interactors$cell, levels = c("Macrophage1", "Macrophage3", "Proliferating", "Ciliated Cells", "Type I Pneumocytes", "Type II Pneumocytes"))
uso1 <- grab1(net, "ENSG00000138768")
uso1_interactors <- cell[cell$Gene %in% node[node$ensemblID %in% as_ids(V(uso1[[1]])), ]$symbol, ]
names(uso1_interactors) <- c("gene", "cell", "pvalue", "adjP", "FC", "expresPercent", "expressPercentRest")
uso1_interactors_sub <- uso1_interactors[, c(1, 2, 6)]
uso1_interactors_melt <- melt(uso1_interactors_sub)
uso1_interactors_melt$cell <- factor(uso1_interactors$cell, levels = c("Mast", "Macrophage1", "Macrophage3", "Proliferating", "Fibroblast2", "Ciliated Cells", "Type II Pneumocytes"))

pdf("inet_s2r1_interactors.pdf", width = 10, height = 8)
ggplot(krt8_interactors_melt, aes(x = gene, y = cell)) +
    geom_point(aes(size = value * 100)) +
    labs(x = "interactor", y = "cell") +
    guides(size = guide_legend(title = "% Expression")) +
    theme_bw(base_size = 20) +
    theme(axis.text.x = element_text(angle = 45, color = "black", hjust = 1, vjust = 1),
        axis.text.y = element_text(color = "black")) +
    theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
ggplot(uso1_interactors_melt, aes(x = gene, y = cell)) +
    geom_point(aes(size = value * 100)) +
    labs(x = "interactor", y = "cell") +
    guides(size = guide_legend(title = "% Expression")) +
    theme_bw(base_size = 20) +
    theme(axis.text.x = element_text(angle = 45, color = "black", hjust = 1, vjust = 1),
        axis.text.y = element_text(color = "black")) +
    theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
dev.off()

###########################################################
# check Roth fBFG screen result with previous dataset
## 1. The lung cells from human COVID19 confirmed case
## ref: https://doi.org/10.1016/j.cell.2020.04.035

roth <- read.xlsx("~/Documents/INET-work/virus_network/references/Roth_fBFG/20200811_DK_VirHostome.xlsx", start = 6)
names(roth) <- c("source", "target")
roth_cell <- cell[cell$Gene %in% roth$target, ]
names(roth_cell) <- c("gene", "cell", "pvalue", "adjP", "FC", "expresPercent", "expressPercentRest")
roth_cell_sub <- roth_cell[, c(1, 2, 6)]
roth_cell_melt <- melt(roth_cell_sub)
roth_cell_melt$cell <- factor(roth_cell$cell, levels = c("Mast", "Monocyte1", "Monocyte2", "Neutrophil", "Proliferating", "Endothelial/Lymphatics", "Fibroblast1", "Fibroblast2", "Macrophage1", "Macrophage3", "T", "Ciliated Cells", "Type I Pneumocytes", "Type II Pneumocytes"))
pdf("Roth_cell_profile.pdf", width = 15, height = 10)
ggplot(roth_cell_melt, aes(x = gene, y = cell)) +
    geom_point(aes(size = value * 100)) +
    labs(x = "interactor", y = "cell", title = "The Roth lab's fBFG screen interactors") +
    guides(size = guide_legend(title = "% Expression")) +
    theme_bw(base_size = 16) +
    theme(axis.text.x = element_text(angle = 45, color = "black", hjust = 1, vjust = 1),
        axis.text.y = element_text(color = "black")) +
    theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
dev.off()

## 2. tissue preferentially expression
roth$target_ensemblID <- node[match(roth$target, node$symbol), "ensemblID"]
roth_tissue <- huri_tissue[huri_tissue$Ensembl_gene_id %in% roth$target_ensemblID, ]
row.names(roth_tissue) <- roth[match(roth_tissue$Ensembl_gene_id, roth$target_ensemblID), "target"]
roth_tissue <- roth_tissue[, -1]
pdf("roth_tissue.pdf", width = 10, height = 7)
Heatmap(as.matrix(roth_tissue[c("CD9", "AKAP9", "CCDC113"), ]),
    name = "z-score",
    colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
    column_title = "tissue preferentially expression",
    width = unit(8, "inch"), height = unit(2, "inch"))
dev.off()

## 3. check the interactors of roth also expressed in SARS2 infected cells or not
### ARRDC3, ENSG00000113369
### CD9, ENSG00000010278
arrdc3 <- grab1(net, "ENSG00000113369")
arrdc3_interactors <- cell[cell$Gene %in% node[node$ensemblID %in% as_ids(V(arrdc3[[1]])), ]$symbol, ]
names(arrdc3_interactors) <- c("gene", "cell", "pvalue", "adjP", "FC", "expresPercent", "expressPercentRest")
arrdc3_interactors_sub <- arrdc3_interactors[, c(1, 2, 6)]
arrdc3_interactors_melt <- melt(arrdc3_interactors_sub)
arrdc3_interactors_melt$cell <- factor(arrdc3_interactors_melt$cell, levels = names(table(arrdc3_interactors_melt$cell))[c(2:11, 1, 12, 13)])
cd9 <- grab1(net, "ENSG00000010278")
cd9_interactors <- cell[cell$Gene %in% node[node$ensemblID %in% as_ids(V(cd9[[1]])), ]$symbol, ]
names(cd9_interactors) <- c("gene", "cell", "pvalue", "adjP", "FC", "expresPercent", "expressPercentRest")
cd9_interactors_sub <- cd9_interactors[, c(1, 2, 6)]
cd9_interactors_melt <- melt(cd9_interactors_sub)
cd9_interactors_melt$cell <- factor(cd9_interactors$cell, levels = names(table(cd9_interactors_melt$cell))[c(2:5, 1, 6, 7)])

pdf("roth_interactors.pdf", width = 10, height = 8)
ggplot(arrdc3_interactors_melt, aes(x = gene, y = cell)) +
    geom_point(aes(size = value * 100)) +
    labs(x = "interactor", y = "cell") +
    guides(size = guide_legend(title = "% Expression")) +
    theme_bw(base_size = 20) +
    theme(axis.text.x = element_text(angle = 45, color = "black", hjust = 1, vjust = 1),
        axis.text.y = element_text(color = "black")) +
    theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
# ggplot(cd9_interactors_melt, aes(x = gene, y = cell)) +
#     geom_point(aes(size = value * 100)) +
#     labs(x = "interactor", y = "cell") +
#     guides(size = guide_legend(title = "% Expression")) +
#     theme_bw(base_size = 20) +
#     theme(axis.text.x = element_text(angle = 45, color = "black", hjust = 1, vjust = 1),
#         axis.text.y = element_text(color = "black")) +
#     theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
dev.off()

###########################################################
# check Gordon APMS screen result with previous dataset
## 1. The lung cells from human COVID19 confirmed case
## ref: https://doi.org/10.1016/j.cell.2020.04.035

gordon <- read.xlsx("~/Documents/INET-work/virus_network/references/Gordon_APMS/41586_2020_2286_MOESM6_ESM.xlsx", start = 2)
gordon <- gordon[, c(1:3)]
names(gordon) <- c("source", "target_ID", "target")
gordon_cell <- cell[cell$Gene %in% gordon$target, ]
## 110 of 332 shown in SARS-CoV-2 treated tissue / cell
names(gordon_cell) <- c("gene", "cell", "pvalue", "adjP", "FC", "expresPercent", "expressPercentRest")
gordon_cell_sub <- gordon_cell[, c(1, 2, 6)]
gordon_cell_melt <- melt(gordon_cell_sub)
gordon_cell_melt$cell <- factor(gordon_cell$cell, levels = names(table(gordon_cell_melt$cell))[c(2:13, 1, 14:15)])
pdf("gordon_cell_profile.pdf", width = 16, height = 10)
ggplot(gordon_cell_melt, aes(x = gene, y = cell)) +
    geom_point(aes(size = value * 100)) +
    labs(x = "interactor", y = "cell", title = "The interactors of viral-host interactome from Gordon's AP-MS screening") +
    guides(size = guide_legend(title = "% Expression")) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 90, color = "black", hjust = 1, vjust = .5),
        axis.text.y = element_text(color = "black")) +
    theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
dev.off()

## 2. tissue preferentially expression
## ! none of the ATP1B1 or HECDT1 shown in tissue preferentially expression data.

###########################################################
# check Li APMS screen result with previous dataset
## 1. The lung cells from human COVID19 confirmed case
## ref: https://doi.org/10.1016/j.cell.2020.04.035

li <- read.xlsx("~/Documents/INET-work/virus_network/references/Li_virus-host/1-s2.0-S2666634020300155-mmc2.xlsx", start = 4)
li <- li[, c(1, 4)]
names(li) <- c("source", "target")
li_cell <- cell[cell$Gene %in% li$target, ]
## 109 of 285 shown in SARS-CoV-2 treated tissue / cell
names(li_cell) <- c("gene", "cell", "pvalue", "adjP", "FC", "expresPercent", "expressPercentRest")
li_cell_sub <- li_cell[, c(1, 2, 6)]
li_cell_melt <- melt(li_cell_sub)
li_cell_melt$cell <- factor(li_cell$cell, levels = names(table(li_cell_melt$cell))[c(2:13, 1, 14:15)])
pdf("li_cell_profile.pdf", width = 16, height = 10)
ggplot(li_cell_melt, aes(x = gene, y = cell)) +
    geom_point(aes(size = value * 100)) +
    labs(x = "interactor", y = "cell", title = "The interactors of viral-host interactome from Li's AP-MS screening") +
    guides(size = guide_legend(title = "% Expression")) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 90, color = "black", hjust = 1, vjust = .5),
        axis.text.y = element_text(color = "black")) +
    theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
dev.off()

## 2. tissue preferentially expression
## ! none of the ATP1B1 or HECDT1 shown in tissue preferentially expression data.
li$target_ensemblID <- node[match(li$target, node$symbol), "ensemblID"]
li_tissue <- huri_tissue[huri_tissue$Ensembl_gene_id %in% li$target_ensemblID, ]
row.names(li_tissue) <- li[match(li_tissue$Ensembl_gene_id, li$target_ensemblID), "target"]
li_tissue <- li_tissue[, -1]
pdf("li_tissue.pdf", width = 7, height = 7)
Heatmap(as.matrix(li_tissue["HECTD1", ]),
    name = "z-score",
    colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
    column_title = "tissue preferentially expression",
    width = unit(5, "inch"), height = unit(.5, "inch"))
dev.off()
