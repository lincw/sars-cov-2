# HuSCI custom GO analysis
# Lin Chung-wen

######
# load package
library(gprofiler2)
library(ggplot2)
library(openxlsx)
library(dplyr)

######
# set up function
## GO enrichment analysis
funcEnrich <- function(query_list, organism, correction, bg = NULL, domain) {
    goquery <- gost(query = query_list, organism = organism, correction_method = correction, evcodes = TRUE, custom_bg = bg, domain_scope = domain)
    goquery$result$inCommunity <- paste0( goquery$meta$query_metadata$queries$query_1, collapse = ",")
    goquery$result$annotatedInCommunity <- paste0(goquery$meta$genes_metadata$query$query_1$ensgs, collapse = ",")
    target <- query_list[query_list %in% husci_node]
    goquery$result$viralTarget <- paste0(target, collapse = ",")
    return(goquery$result)
}

## merge into data frame
toDataFrame <- function(x) {
    husci_all_plot <- data.frame(
        "term" = x[, "term_name"],
        "term_size" = x[, "term_size"],
        "observed" = x[, "precision"],
        "background" = x[, "term_size"] / x[, "effective_domain_size"],
        "observeRatio" = x[, "precision"] / (x[, "term_size"] / x[, "effective_domain_size"]))
    return(husci_all_plot)
}
# plot
metaPlot <- function(x, main, size = 8) {
  x$term2 <- factor(x$term2, levels = x$term2[order(x$observeRatio)])
  ggplot(x, aes(x = term2, y = observeRatio)) +
        geom_bar(stat = "identity", fill = rgb(128, 41, 227, maxColorValue = 255)) +
        coord_flip() +
        labs(y = "Effect size", x = "Functional term (p < 0.05)", title = main) +
        theme_bw() +
        theme(axis.text = element_text(color = "black", size = size, hjust = 0),
        axis.title = element_text(size = size, color = "black", face = "bold"))
}

######
# load different GO background
# for HuSCI
horf_bp_all <- upload_GMT_file(gmtfile = "~/workplace/database/Homo_sapiens/GO/hORFeome_gobp_2706.gmt")
horf_mf_all <- upload_GMT_file(gmtfile = "~/workplace/database/Homo_sapiens/GO/hORFeome_gomf_2706.gmt")
horf_cc_all <- upload_GMT_file(gmtfile = "~/workplace/database/Homo_sapiens/GO/hORFeome_gocc_2706.gmt")

# for all others
bp_all <- upload_GMT_file(gmtfile = "~/workplace/database/Homo_sapiens/GO/gobp_2706.gmt")
mf_all <- upload_GMT_file(gmtfile = "~/workplace/database/Homo_sapiens/GO/gomf_2706.gmt")
cc_all <- upload_GMT_file(gmtfile = "~/workplace/database/Homo_sapiens/GO/gocc_2706.gmt")

# for all others, GMT from gProfiler
bp_gprofiler <- upload_GMT_file(gmtfile = "~/workplace/database/Homo_sapiens/GO/gprofiler_hsapiens.name_29062021/hsapiens.GO:BP.name.gmt")
mf_gprofiler <- upload_GMT_file(gmtfile = "~/workplace/database/Homo_sapiens/GO/gprofiler_hsapiens.name_29062021/hsapiens.GO:MF.name.gmt")
cc_gprofiler <- upload_GMT_file(gmtfile = "~/workplace/database/Homo_sapiens/GO/gprofiler_hsapiens.name_29062021/hsapiens.GO:CC.name.gmt")

######
# load screening space of HuSCI
space <- read.xlsx("/Volumes/GoogleDrive/My\ Drive/VirHostome_CW/GitHub/data/extended_table/Extended_Table_1_search_space.xlsx", sheet = 2)

######
# load PPI data
husci <- read.csv("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/binary_node_1126_simple.csv", header = T)
husci_node <- unique(husci[husci$category == "human", "node"])

ext_table2 <- "/Volumes/GoogleDrive/My\ Drive/VirHostome_CW/GitHub/data/extended_table/Extended_Table_2_PPIs.xlsx"
gordon <- read.xlsx(ext_table2, sheet = 3)
gordon_node <- unique(gordon$PreyGene)

stukalov <- read.xlsx(ext_table2, sheet = 4)
stukalov_node <- unique(stukalov$human)

li <- read.xlsx(ext_table2, sheet = 5)
li_node <- unique(li$human)

nabeel <- read.xlsx(ext_table2, sheet = 6)
nabeel_node <- unique(nabeel$PreyGene)

laurent <- read.xlsx(ext_table2, sheet = 7)
laurent_node <- unique(laurent$human)

stgermain <- read.xlsx(ext_table2, sheet = 8)
stgermain_node <- unique(stgermain$PreyGene)

samavarchi <- read.xlsx(ext_table2, sheet = 9)
samavarchi_node <- unique(samavarchi$PreyGene)

######
# GO analysis
# 1. HuSCI
# based on all evidence codes
husci_bp <- funcEnrich(husci_node, bp_gprofiler, "fdr", space$ensembl_gene_name)
husci_mf <- funcEnrich(husci_node, mf_gprofiler, "fdr", space$ensembl_gene_name)
husci_cc <- funcEnrich(husci_node, cc_gprofiler, "fdr", space$ensembl_gene_name)
husci_all <- list(BP = husci_bp, MF = husci_mf, CC = husci_cc)
husci_all_df <- do.call(rbind, husci_all)[, c(1:13, 15:19)]
write.csv(husci_all_df, file = "~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/GO/husci_go_fdr_2906.csv", row.names = F)
## plot
husci_all_plot <- toDataFrame(husci_all_df)
metaPlot(husci_all_plot, "HuSCI")
ggsave("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/GO/husci_go_fdr_2706.pdf")

# 2. Gordon
# based on Experimental evidence codes
gordon_bp <- funcEnrich(gordon_node, bp_gprofiler, "fdr", bg = NULL, "annotated")
gordon_mf <- funcEnrich(gordon_node, mf_gprofiler, "fdr", bg = NULL, "annotated")
gordon_cc <- funcEnrich(gordon_node, cc_gprofiler, "fdr", bg = NULL, "annotated")
gordon_go <- list(BP = gordon_bp, MF = gordon_mf, CC = gordon_cc)
gordon_go_df <- do.call(rbind, gordon_go)[, c(1:13, 15:19)]
write.csv(do.call(rbind, gordon_go)[, c(1:13, 15:19)], file = "~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/GO/gordon_go_fdr_2706.csv", row.names = F)

## plot
gordon_go_plot <- toDataFrame(gordon_go_df)
metaPlot(gordon_go_plot, "Gordon et al")
ggsave("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/GO/gordon_go_fdr_2706.pdf")

# 3. Stukalov
# based on Experimental evidence codes
stukalov_bp <- funcEnrich(stukalov_node, bp_all, "fdr")
stukalov_mf <- funcEnrich(stukalov_node, mf_all, "fdr")
stukalov_cc <- funcEnrich(stukalov_node, cc_all, "fdr")
stukalov_go <- list(BP = stukalov_bp, MF = stukalov_mf, CC = stukalov_cc)
stukalov_go_df <- do.call(rbind, stukalov_go)[, c(1:13, 15:19)]
write.csv(do.call(rbind, stukalov_go)[, c(1:13, 15:19)], file = "~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/GO/stukalov_go_fdr_2706.csv", row.names = F)

## plot
stukalov_go_plot <- toDataFrame(stukalov_go_df)
metaPlot(stukalov_go_plot, "Stukalov et al")
ggsave("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/GO/stukalov_go_fdr_2706.pdf")

# 4. Li
# based on Experimental evidence codes
li_bp <- funcEnrich(li_node, bp_all, "fdr")
li_mf <- funcEnrich(li_node, mf_all, "fdr")
li_cc <- funcEnrich(li_node, cc_all, "fdr")
li_go <- list(BP = li_bp, MF = li_mf, CC = li_cc)
li_go_df <- do.call(rbind, li_go)[, c(1:13, 15:19)]
write.csv(do.call(rbind, li_go)[, c(1:13, 15:19)], file = "~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/GO/li_go_fdr_2706.csv", row.names = F)

## plot
li_go_plot <- toDataFrame(li_go_df)
metaPlot(li_go_plot, "Li et al")
ggsave("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/GO/li_go_fdr_2706.pdf")

# 5. Nabeel
# based on Experimental evidence codes
nabeel_bp <- funcEnrich(nabeel_node, bp_all, "fdr")
nabeel_mf <- funcEnrich(nabeel_node, mf_all, "fdr")
nabeel_cc <- funcEnrich(nabeel_node, cc_all, "fdr")
nabeel_go <- list(BP = nabeel_bp, MF = nabeel_mf, CC = nabeel_cc)
nabeel_go_df <- do.call(rbind, nabeel_go)[, c(1:13, 15:19)]
write.csv(do.call(rbind, nabeel_go)[, c(1:13, 15:19)], file = "~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/GO/nabeel_go_fdr_2706.csv", row.names = F)

## plot
nabeel_go_plot <- toDataFrame(nabeel_go_df)
metaPlot(nabeel_go_plot, "Nabeel et al")
ggsave("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/GO/nabeel_go_fdr_2706.pdf")

# 6. Laurent
# based on Experimental evidence codes
laurent_bp <- funcEnrich(laurent_node, bp_all, "fdr")
laurent_mf <- funcEnrich(laurent_node, mf_all, "fdr")
laurent_cc <- funcEnrich(laurent_node, cc_all, "fdr")
laurent_go <- list(BP = laurent_bp, MF = laurent_mf, CC = laurent_cc)
laurent_go_df <- do.call(rbind, laurent_go)[, c(1:13, 15:19)]
write.csv(do.call(rbind, laurent_go)[, c(1:13, 15:19)], file = "~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/GO/laurent_go_fdr_2706.csv", row.names = F)

## plot
laurent_go_plot <- toDataFrame(laurent_go_df)
metaPlot(laurent_go_plot, "Laurent et al")
ggsave("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/GO/laurent_go_fdr_2706.pdf")

# 7. St-Germain
# based on Experimental evidence codes
stgermain_bp <- funcEnrich(stgermain_node, bp_all, "fdr")
stgermain_mf <- funcEnrich(stgermain_node, mf_all, "fdr")
stgermain_cc <- funcEnrich(stgermain_node, cc_all, "fdr")
stgermain_go <- list(BP = stgermain_bp, MF = stgermain_mf, CC = stgermain_cc)
stgermain_go_df <- do.call(rbind, stgermain_go)[, c(1:13, 15:19)]
write.csv(do.call(rbind, stgermain_go)[, c(1:13, 15:19)], file = "~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/GO/stgermain_go_fdr_2706.csv", row.names = F)

## plot
stgermain_go_plot <- toDataFrame(stgermain_go_df)
metaPlot(stgermain_go_plot, "St-Germain et al")
ggsave("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/GO/stgermain_go_fdr_2706.pdf")

# 8. Samavarchi
# based on Experimental evidence codes
samavarchi_bp <- funcEnrich(samavarchi_node, bp_all, "fdr")
samavarchi_mf <- funcEnrich(samavarchi_node, mf_all, "fdr")
samavarchi_cc <- funcEnrich(samavarchi_node, cc_all, "fdr")
samavarchi_go <- list(BP = samavarchi_bp, MF = samavarchi_mf, CC = samavarchi_cc)
samavarchi_go_df <- do.call(rbind, samavarchi_go)[, c(1:13, 15:19)]
write.csv(do.call(rbind, samavarchi_go)[, c(1:13, 15:19)], file = "~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/GO/samavarchi_go_fdr_2706.csv", row.names = F)

## plot
samavarchi_go_plot <- toDataFrame(samavarchi_go_df)
metaPlot(samavarchi_go_plot, "Samavarchi et al")
ggsave("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/GO/samavarchi_go_fdr_2706.pdf")

# save RData
save.image("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/GO/custom_GO_2706.RData")

######
# exclude above, using gProfiler web server, via Benny
husci_gprofiler <- read.csv("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/GO/gProfiler_hsapiens_9-20-2021_11-03-20 AM__intersections.csv", skip = 17, header = T)
husci_gprofiler <- husci_gprofiler[husci_gprofiler$adjusted_p_value < 0.05, c(1:10)]

husci_toPlot <- data.frame(
    term = husci_gprofiler[, 2],
    term2 = husci_gprofiler[, 2],
    intersection = husci_gprofiler[, 8],
    query = husci_gprofiler[, 8] / husci_gprofiler[, 7],
    background = husci_gprofiler[, 6] / husci_gprofiler[, 9],
    observeRatio = (husci_gprofiler[, 8] / husci_gprofiler[, 7]) / (husci_gprofiler[, 6] / husci_gprofiler[, 9])
)


pdf("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/GO/HuSCI_GO_gProfiler_20210920.pdf", height = 4)
# metaPlot(husci_toPlot, "no IEA")
# metaPlot(husci_toPlot[grepl("GO", husci_toPlot$term), ], "no IEA, GO only")
# metaPlot(husci_toPlot[husci_gprofiler[, 6] > 2,], "no IEA, term_size > 2")
# metaPlot(husci_toPlot[grepl("GO", husci_toPlot$term) & husci_gprofiler[, 6] > 4, ], "no IEA, GO only\nterm_size > 4")
metaPlot(husci_toPlot[husci_gprofiler[, 6] > 4,], "HuSCI functional analysis")
# metaPlot(husci_toPlot[grepl("GO", husci_toPlot$term) & husci_gprofiler[, 6] > 4,], "no IEA, GO only\nterm_size > 4")
dev.off()

######
# 2nd GO figure for 20 GWAS hits
metaPlot2 <- function(x, main, size = 8, p) {
    ggplot(x, aes(x = term2, y = observeRatio)) +
        geom_bar(stat = "identity", fill = rgb(128, 41, 227, maxColorValue = 255)) +
        coord_flip() +
        labs(y = "Effect size", x = paste0("Functional term (p < ", p, ")"), title = main) +
        theme_bw() +
        theme(axis.text = element_text(color = "black", size = size, hjust = 0),
        axis.title = element_text(size = size, color = "black", face = "bold"))
}
gwas20 <- read.csv("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/GO/DY_gProfiler_hsapiens_7-20-2021_5-32-27 PM__intersections.csv", skip = 17, header = T)
gwas20 <- gwas20[, c(1:10)]

gwas20 <- data.frame(
    term = paste0(gwas20, "_", gwas20[, 2]),
    term2 = gwas20[, 2],
    intersection = gwas20[, 8],
    query = gwas20[, 8] / gwas20[, 7],
    background = gwas20[, 6] / gwas20[, 9],
    observeRatio = (gwas20[, 8] / gwas20[, 7]) / (gwas20[, 6] / gwas20[, 9])
)
gwas20$term2 <- factor(gwas20$term2, levels = rev(factor(gwas20$term2)[c(12, 9, 6, 7, 8, 11, 10, 13, 5, 4, 3, 2, 1)]))

pdf("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/GO/HuSCI_GO_gProfiler_GWAS20.pdf", height = 4)
metaPlot2(gwas20, "Viral targets from 20 GWAS hits", p = "0.1")
dev.off()

gwas20_husci <- read.csv("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/GO/DY_gProfiler_hsapiens_7-20-2021_5-32-27 PM__intersections_HuSCI_BG.csv", skip = 17, header = T)
gwas20_husci <- gwas20_husci[, c(1:10)]
gwas20_husci <- data.frame(
    term = paste0(gwas20_husci, "_", gwas20_husci[, 2]),
    term2 = gwas20_husci[, 2],
    intersection = gwas20_husci[, 8],
    query = gwas20_husci[, 8] / gwas20_husci[, 7],
    background = gwas20_husci[, 6] / gwas20_husci[, 9],
    observeRatio = (gwas20_husci[, 8] / gwas20_husci[, 7]) / (gwas20_husci[, 6] / gwas20_husci[, 9])
)
gwas20_husci$term2 <- factor(gwas20_husci$term2, levels = rev(factor(gwas20_husci$term2)[c(5,7,11,9,8,10,6,2,1,3,4)]))

pdf("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/GO/HuSCI_GO_gProfiler_GWAS20_husci.pdf", height = 4)
metaPlot2(gwas20_husci, "Viral targets from 20 GWAS hits\n(background: HuSCI n = 171)", p = "0.1")
dev.off()
