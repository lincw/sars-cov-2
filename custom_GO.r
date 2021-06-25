# HuSCI custom GO analysis
# Lin Chung-wen

######
# load package
library(gprofiler2)
library(ggplot2)
library(openxlsx)

######
# set up function
## GO enrichment analysis
funcEnrich <- function(query_list, organism, correction, significant_q = TRUE) {
    goquery <- gost(query = query_list, organism = organism, correction_method = correction, evcodes = TRUE, significant = significant_q)
    goquery$result$inCommunity <- paste0( goquery$meta$query_metadata$queries$query_1, collapse = ",")
    goquery$result$annotatedInCommunity <- paste0(goquery$meta$genes_metadata$query$query_1$ensgs, collapse = ",")
    target <- query_list[query_list %in% husci_node]
    goquery$result$viralTarget <- paste0(target, collapse = ",")
    return(goquery$result)
}

## merge into data frame
toDataFrame <- function(x) {
    rows <- which(x[, "term_size"] > 3)
    husci_all_plot <- data.frame(
        "term" = x[rows, "term_name"],
        "observed" = x[rows, "precision"],
        "background" = x[rows, "term_size"] / x[rows, "effective_domain_size"],
        "observeRatio" = x[rows, "precision"] / (x[rows, "term_size"] / x[rows, "effective_domain_size"]))
    return(husci_all_plot)
}
# plot
metaPlot <- function(x, main, size = 8) {
    x$term <- factor(x$term, levels = x$term[order(x$observeRatio)])
    ggplot(x, aes(x = term, y = observeRatio)) +
        geom_bar(stat = "identity", fill = rgb(128, 41, 227, maxColorValue = 255)) +
        coord_flip() +
        labs(y = "Effect size (%)", x = "Functional term (p < 0.05)", title = main) +
        theme_bw() +
        theme(axis.text = element_text(color = "black", size = size, hjust = 0),
        axis.title = element_text(size = size, color = "black", face = "bold"))
}

######
# load different GO background
# for HuSCI
horf_bp_all <- upload_GMT_file(gmtfile = "~/workplace/database/Homo_sapiens/GO/hORFeome_gobp.gmt")
horf_mf_all <- upload_GMT_file(gmtfile = "~/workplace/database/Homo_sapiens/GO/hORFeome_gomf.gmt")
horf_cc_all <- upload_GMT_file(gmtfile = "~/workplace/database/Homo_sapiens/GO/hORFeome_gocc.gmt")

horf_bp_noIEA <- upload_GMT_file(gmtfile = "~/workplace/database/Homo_sapiens/GO/hORFeome_gobp_noIEA.gmt")
horf_mf_noIEA <- upload_GMT_file(gmtfile = "~/workplace/database/Homo_sapiens/GO/hORFeome_gomf_noIEA.gmt")
horf_cc_noIEA <- upload_GMT_file(gmtfile = "~/workplace/database/Homo_sapiens/GO/hORFeome_gocc_noIEA.gmt")

horf_bp_exp <- upload_GMT_file(gmtfile = "~/workplace/database/Homo_sapiens/GO/hORFeome_gobp_exp.gmt")
horf_mf_exp <- upload_GMT_file(gmtfile = "~/workplace/database/Homo_sapiens/GO/hORFeome_gomf_exp.gmt")
horf_cc_exp <- upload_GMT_file(gmtfile = "~/workplace/database/Homo_sapiens/GO/hORFeome_gocc_exp.gmt")

horf_bp_exp_phy <- upload_GMT_file(gmtfile = "~/workplace/database/Homo_sapiens/GO/hORFeome_gobp_exp_phy.gmt")
horf_mf_exp_phy <- upload_GMT_file(gmtfile = "~/workplace/database/Homo_sapiens/GO/hORFeome_gomf_exp_phy.gmt")
horf_cc_exp_phy <- upload_GMT_file(gmtfile = "~/workplace/database/Homo_sapiens/GO/hORFeome_gocc_exp_phy.gmt")

# for all others
bp_all <- upload_GMT_file(gmtfile = "~/workplace/database/Homo_sapiens/GO/gobp.gmt")
mf_all <- upload_GMT_file(gmtfile = "~/workplace/database/Homo_sapiens/GO/gomf.gmt")
cc_all <- upload_GMT_file(gmtfile = "~/workplace/database/Homo_sapiens/GO/gocc.gmt")

bp_noIEA <- upload_GMT_file(gmtfile = "~/workplace/database/Homo_sapiens/GO/gobp_noIEA.gmt")
mf_noIEA <- upload_GMT_file(gmtfile = "~/workplace/database/Homo_sapiens/GO/gomf_noIEA.gmt")
cc_noIEA <- upload_GMT_file(gmtfile = "~/workplace/database/Homo_sapiens/GO/gocc_noIEA.gmt")

bp_exp <- upload_GMT_file(gmtfile = "~/workplace/database/Homo_sapiens/GO/gobp_exp.gmt")
mf_exp <- upload_GMT_file(gmtfile = "~/workplace/database/Homo_sapiens/GO/gomf_exp.gmt")
cc_exp <- upload_GMT_file(gmtfile = "~/workplace/database/Homo_sapiens/GO/gocc_exp.gmt")

bp_exp_phy <- upload_GMT_file(gmtfile = "~/workplace/database/Homo_sapiens/GO/gobp_exp_phy.gmt")
mf_exp_phy <- upload_GMT_file(gmtfile = "~/workplace/database/Homo_sapiens/GO/gomf_exp_phy.gmt")
cc_exp_phy <- upload_GMT_file(gmtfile = "~/workplace/database/Homo_sapiens/GO/gocc_exp_phy.gmt")

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
husci_bp <- funcEnrich(husci_node, horf_bp_all, "fdr")
husci_mf <- funcEnrich(husci_node, horf_mf_all, "fdr")
husci_cc <- funcEnrich(husci_node, horf_cc_all, "fdr")
husci_all <- list(BP = husci_bp, MF = husci_mf, CC = husci_cc)
husci_all_df <- do.call(rbind, husci_all)[, c(1:13, 15:19)]
write.csv(husci_all_df, file = "~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/GO/husci_all_fdr.csv", row.names = F)
## plot
husci_all_plot <- toDataFrame(husci_all_df)
metaPlot(husci_all_plot, "HuSCI vs all ECs")
ggsave("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/GO/husci_fdr_all.pdf")

# based on all evidence codes, IEA excluded
husci_bp_noIEA <- funcEnrich(husci_node, horf_bp_noIEA, "fdr")
husci_mf_noIEA <- funcEnrich(husci_node, horf_mf_noIEA, "fdr")
husci_cc_noIEA <- funcEnrich(husci_node, horf_cc_noIEA, "fdr")
husci_noIEA <- list(BP = husci_bp_noIEA, MF = husci_mf_noIEA, CC = husci_cc_noIEA)
husci_noIEA_df <- do.call(rbind, husci_noIEA)[, c(1:13, 15:19)]
write.csv(husci_noIEA_df, file = "~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/GO/husci_noIEA_fdr.csv", row.names = F)
## plot
husci_noIEA_plot <- toDataFrame(husci_noIEA_df)
metaPlot(husci_noIEA_plot, "HuSCI vs no IEA")
ggsave("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/GO/husci_fdr_noIEA.pdf")

# based on Experimental evidence codes
husci_bp_exp <- funcEnrich(husci_node, horf_bp_exp, "fdr")
husci_mf_exp <- funcEnrich(husci_node, horf_mf_exp, "fdr")
husci_cc_exp <- funcEnrich(husci_node, horf_cc_exp, "fdr") # no significant terms
husci_exp <- list(BP = husci_bp_exp, MF = husci_mf_exp)
husci_exp_df <- do.call(rbind, husci_exp)[, c(1:13, 15:19)]
write.csv(do.call(rbind, husci_exp)[, c(1:13, 15:19)], file = "~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/GO/husci_exp_fdr.csv", row.names = F)

## plot
husci_exp_plot <- toDataFrame(husci_exp_df)
metaPlot(husci_exp_plot, "HuSCI vs EXP")
ggsave("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/GO/husci_fdr_exp.pdf")

# based on Experimental evidence codes and phylogenetically-inferred  annotations
husci_bp_exp_phy <- funcEnrich(husci_node, horf_bp_exp_phy, "fdr")
husci_mf_exp_phy <- funcEnrich(husci_node, horf_mf_exp_phy, "fdr")
husci_cc_exp_phy <- funcEnrich(husci_node, horf_cc_exp_phy, "fdr")
husci_exp_phy <- list(BP = husci_bp_exp_phy, MF = husci_mf_exp_phy, CC = husci_cc_exp_phy)
husci_exp_phy_df <- do.call(rbind, husci_exp_phy)[, c(1:13, 15:19)]
write.csv(do.call(rbind, husci_exp_phy)[, c(1:13, 15:19)], file = "~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/GO/husci_exp_phy_fdr.csv", row.names = F)

## plot
husci_exp_phy_plot <- toDataFrame(husci_exp_phy_df)
metaPlot(husci_exp_phy_plot, "HuSCI vs EXP PHY")
ggsave("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/GO/husci_exp_phy_fdr.pdf")

# 2. Gordon
# based on Experimental evidence codes
gordon_bp_exp <- funcEnrich(gordon_node, bp_exp, "fdr")
gordon_mf_exp <- funcEnrich(gordon_node, mf_exp, "fdr")
gordon_cc_exp <- funcEnrich(gordon_node, cc_exp, "fdr")
gordon_exp <- list(BP = gordon_bp_exp, MF = gordon_mf_exp, CC = gordon_cc_exp)
gordon_exp_df <- do.call(rbind, gordon_exp)[, c(1:13, 15:19)]
# write.csv(do.call(rbind, gordon_exp)[, c(1:13, 15:19)], file = "~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/GO/gordon_exp_fdr.csv", row.names = F)

## plot
gordon_exp_plot <- toDataFrame(gordon_exp_df)
metaPlot(gordon_exp_plot, "gordon vs EXP")
ggsave("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/GO/gordon_fdr_exp.pdf")

# based on Experimental evidence codes and phylogenetically-inferred  annotations
gordon_bp_exp_phy <- funcEnrich(gordon_node, bp_exp_phy, "fdr")
gordon_mf_exp_phy <- funcEnrich(gordon_node, mf_exp_phy, "fdr")
gordon_cc_exp_phy <- funcEnrich(gordon_node, cc_exp_phy, "fdr")
gordon_exp_phy <- list(BP = gordon_bp_exp_phy, MF = gordon_mf_exp_phy, CC = gordon_cc_exp_phy)
gordon_exp_phy_df <- do.call(rbind, gordon_exp_phy)[, c(1:13, 15:19)]
# write.csv(do.call(rbind, gordon_exp_phy)[, c(1:13, 15:19)], file = "~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/GO/gordon_exp_phy_fdr.csv", row.names = F)

## plot
gordon_exp_phy_plot <- toDataFrame(gordon_exp_phy_df)
metaPlot(gordon_exp_phy_plot, "gordon vs EXP PHY", 6)
ggsave("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/GO/gordon_fdr_exp_phy.pdf")

# 3. Stukalov
# based on Experimental evidence codes
stukalov_bp_exp <- funcEnrich(stukalov_node, bp_exp, "fdr")
stukalov_mf_exp <- funcEnrich(stukalov_node, mf_exp, "fdr")
stukalov_cc_exp <- funcEnrich(stukalov_node, cc_exp, "fdr")
stukalov_exp <- list(BP = stukalov_bp_exp, MF = stukalov_mf_exp, CC = stukalov_cc_exp)
stukalov_exp_df <- do.call(rbind, stukalov_exp)[, c(1:13, 15:19)]
# write.csv(do.call(rbind, stukalov_exp)[, c(1:13, 15:19)], file = "~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/GO/stukalov_exp_fdr.csv", row.names = F)

## plot
stukalov_exp_plot <- toDataFrame(stukalov_exp_df)
metaPlot(stukalov_exp_plot, "stukalov vs EXP", 6)
ggsave("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/GO/stukalov_fdr_exp.pdf")

# based on Experimental evidence codes and phylogenetically-inferred  annotations
stukalov_bp_exp_phy <- funcEnrich(stukalov_node, bp_exp_phy, "fdr")
stukalov_mf_exp_phy <- funcEnrich(stukalov_node, mf_exp_phy, "fdr")
stukalov_cc_exp_phy <- funcEnrich(stukalov_node, cc_exp_phy, "fdr")
stukalov_exp_phy <- list(BP = stukalov_bp_exp_phy, MF = stukalov_mf_exp_phy, CC = stukalov_cc_exp_phy)
stukalov_exp_phy_df <- do.call(rbind, stukalov_exp_phy)[, c(1:13, 15:19)]
# write.csv(do.call(rbind, stukalov_exp_phy)[, c(1:13, 15:19)], file = "~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/GO/stukalov_exp_phy_fdr.csv", row.names = F)

## plot
stukalov_exp_phy_plot <- toDataFrame(stukalov_exp_phy_df)
metaPlot(stukalov_exp_phy_plot, "stukalov vs EXP PHY", 5)
ggsave("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/GO/stukalov_fdr_exp_phy.pdf")

# 4. Li
# based on Experimental evidence codes
li_bp_exp <- funcEnrich(li_node, bp_exp, "fdr")
li_mf_exp <- funcEnrich(li_node, mf_exp, "fdr")
li_cc_exp <- funcEnrich(li_node, cc_exp, "fdr")
li_exp <- list(BP = li_bp_exp, MF = li_mf_exp, CC = li_cc_exp)
li_exp_df <- do.call(rbind, li_exp)[, c(1:13, 15:19)]
# write.csv(do.call(rbind, li_exp)[, c(1:13, 15:19)], file = "~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/GO/li_exp_fdr.csv", row.names = F)

## plot
li_exp_plot <- toDataFrame(li_exp_df)
metaPlot(li_exp_plot, "li vs EXP", 6)
ggsave("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/GO/li_fdr_exp.pdf")

# based on Experimental evidence codes and phylogenetically-inferred  annotations
li_bp_exp_phy <- funcEnrich(li_node, bp_exp_phy, "fdr")
li_mf_exp_phy <- funcEnrich(li_node, mf_exp_phy, "fdr")
li_cc_exp_phy <- funcEnrich(li_node, cc_exp_phy, "fdr")
li_exp_phy <- list(BP = li_bp_exp_phy, MF = li_mf_exp_phy, CC = li_cc_exp_phy)
li_exp_phy_df <- do.call(rbind, li_exp_phy)[, c(1:13, 15:19)]
# write.csv(do.call(rbind, li_exp_phy)[, c(1:13, 15:19)], file = "~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/GO/li_exp_phy_fdr.csv", row.names = F)

## plot
li_exp_phy_plot <- toDataFrame(li_exp_phy_df)
metaPlot(li_exp_phy_plot, "li vs EXP PHY", 6)
ggsave("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/GO/li_fdr_exp_phy.pdf")

# 5. Nabeel
# based on Experimental evidence codes
nabeel_bp_exp <- funcEnrich(nabeel_node, bp_exp, "fdr")
nabeel_mf_exp <- funcEnrich(nabeel_node, mf_exp, "fdr")
nabeel_cc_exp <- funcEnrich(nabeel_node, cc_exp, "fdr")
nabeel_exp <- list(BP = nabeel_bp_exp, MF = nabeel_mf_exp, CC = nabeel_cc_exp)
nabeel_exp_df <- do.call(rbind, nabeel_exp)[, c(1:13, 15:19)]
# write.csv(do.call(rbind, nabeel_exp)[, c(1:13, 15:19)], file = "~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/GO/nabeel_exp_fdr.csv", row.names = F)

## plot
nabeel_exp_plot <- toDataFrame(nabeel_exp_df)
metaPlot(nabeel_exp_plot, "li vs EXP", 6)
ggsave("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/GO/nabeel_fdr_exp.pdf")

# based on Experimental evidence codes and phylogenetically-inferred  annotations
nabeel_bp_exp_phy <- funcEnrich(nabeel_node, bp_exp_phy, "fdr")
nabeel_mf_exp_phy <- funcEnrich(nabeel_node, mf_exp_phy, "fdr")
nabeel_cc_exp_phy <- funcEnrich(nabeel_node, cc_exp_phy, "fdr")
nabeel_exp_phy <- list(BP = nabeel_bp_exp_phy, MF = nabeel_mf_exp_phy, CC = nabeel_cc_exp_phy)
nabeel_exp_phy_df <- do.call(rbind, nabeel_exp_phy)[, c(1:13, 15:19)]
# write.csv(do.call(rbind, nabeel_exp_phy)[, c(1:13, 15:19)], file = "~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/GO/nabeel_exp_phy_fdr.csv", row.names = F)

## plot
nabeel_exp_phy_plot <- toDataFrame(nabeel_exp_phy_df)
metaPlot(nabeel_exp_phy_plot, "nabeel vs EXP PHY", 4)
ggsave("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/GO/nabeel_fdr_exp_phy.pdf")

# 6. Laurent
# based on Experimental evidence codes
laurent_bp_exp <- funcEnrich(laurent_node, bp_exp, "fdr")
laurent_mf_exp <- funcEnrich(laurent_node, mf_exp, "fdr")
laurent_cc_exp <- funcEnrich(laurent_node, cc_exp, "fdr")
laurent_exp <- list(BP = laurent_bp_exp, MF = laurent_mf_exp, CC = laurent_cc_exp)
laurent_exp_df <- do.call(rbind, laurent_exp)[, c(1:13, 15:19)]
# write.csv(do.call(rbind, laurent_exp)[, c(1:13, 15:19)], file = "~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/GO/laurent_exp_fdr.csv", row.names = F)

## plot
laurent_exp_plot <- toDataFrame(laurent_exp_df)
metaPlot(laurent_exp_plot, "li vs EXP", 6)
ggsave("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/GO/laurent_fdr_exp.pdf")

# based on Experimental evidence codes and phylogenetically-inferred  annotations
laurent_bp_exp_phy <- funcEnrich(laurent_node, bp_exp_phy, "fdr")
laurent_mf_exp_phy <- funcEnrich(laurent_node, mf_exp_phy, "fdr")
laurent_cc_exp_phy <- funcEnrich(laurent_node, cc_exp_phy, "fdr")
laurent_exp_phy <- list(BP = laurent_bp_exp_phy, MF = laurent_mf_exp_phy, CC = laurent_cc_exp_phy)
laurent_exp_phy_df <- do.call(rbind, laurent_exp_phy)[, c(1:13, 15:19)]
# write.csv(do.call(rbind, laurent_exp_phy)[, c(1:13, 15:19)], file = "~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/GO/laurent_exp_phy_fdr.csv", row.names = F)

## plot
laurent_exp_phy_plot <- toDataFrame(laurent_exp_phy_df)
metaPlot(laurent_exp_phy_plot, "laurent vs EXP PHY", 4)
ggsave("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/GO/laurent_fdr_exp_phy.pdf")

# 7. St-Germain
# based on Experimental evidence codes
stgermain_bp_exp <- funcEnrich(stgermain_node, bp_exp, "fdr")
stgermain_mf_exp <- funcEnrich(stgermain_node, mf_exp, "fdr")
stgermain_cc_exp <- funcEnrich(stgermain_node, cc_exp, "fdr")
stgermain_exp <- list(BP = stgermain_bp_exp, MF = stgermain_mf_exp, CC = stgermain_cc_exp)
stgermain_exp_df <- do.call(rbind, stgermain_exp)[, c(1:13, 15:19)]
# write.csv(do.call(rbind, stgermain_exp)[, c(1:13, 15:19)], file = "~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/GO/stgermain_exp_fdr.csv", row.names = F)

## plot
stgermain_exp_plot <- toDataFrame(stgermain_exp_df)
metaPlot(stgermain_exp_plot, "stgermain vs EXP", 5)
ggsave("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/GO/stgermain_fdr_exp.pdf")

# based on Experimental evidence codes and phylogenetically-inferred  annotations
stgermain_bp_exp_phy <- funcEnrich(stgermain_node, bp_exp_phy, "fdr")
stgermain_mf_exp_phy <- funcEnrich(stgermain_node, mf_exp_phy, "fdr")
stgermain_cc_exp_phy <- funcEnrich(stgermain_node, cc_exp_phy, "fdr")
stgermain_exp_phy <- list(BP = stgermain_bp_exp_phy, MF = stgermain_mf_exp_phy, CC = stgermain_cc_exp_phy)
stgermain_exp_phy_df <- do.call(rbind, stgermain_exp_phy)[, c(1:13, 15:19)]
# write.csv(do.call(rbind, stgermain_exp_phy)[, c(1:13, 15:19)], file = "~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/GO/stgermain_exp_phy_fdr.csv", row.names = F)

## plot
stgermain_exp_phy_plot <- toDataFrame(stgermain_exp_phy_df)
metaPlot(stgermain_exp_phy_plot, "St-Germain vs EXP PHY", 4)
ggsave("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/GO/stgermain_fdr_exp_phy.pdf")

# 8. Samavarchi
# based on Experimental evidence codes
samavarchi_bp_exp <- funcEnrich(samavarchi_node, bp_exp, "fdr")
samavarchi_mf_exp <- funcEnrich(samavarchi_node, mf_exp, "fdr")
samavarchi_cc_exp <- funcEnrich(samavarchi_node, cc_exp, "fdr")
samavarchi_exp <- list(BP = samavarchi_bp_exp, MF = samavarchi_mf_exp, CC = samavarchi_cc_exp)
samavarchi_exp_df <- do.call(rbind, samavarchi_exp)[, c(1:13, 15:19)]
# write.csv(do.call(rbind, samavarchi_exp)[, c(1:13, 15:19)], file = "~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/GO/samavarchi_exp_fdr.csv", row.names = F)

## plot
samavarchi_exp_plot <- toDataFrame(samavarchi_exp_df)
metaPlot(samavarchi_exp_plot, "samavarchi vs EXP", 5)
ggsave("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/GO/samavarchi_fdr_exp.pdf")

# based on Experimental evidence codes and phylogenetically-inferred  annotations
samavarchi_bp_exp_phy <- funcEnrich(samavarchi_node, bp_exp_phy, "fdr")
samavarchi_mf_exp_phy <- funcEnrich(samavarchi_node, mf_exp_phy, "fdr")
samavarchi_cc_exp_phy <- funcEnrich(samavarchi_node, cc_exp_phy, "fdr")
samavarchi_exp_phy <- list(BP = samavarchi_bp_exp_phy, MF = samavarchi_mf_exp_phy, CC = samavarchi_cc_exp_phy)
samavarchi_exp_phy_df <- do.call(rbind, samavarchi_exp_phy)[, c(1:13, 15:19)]
# write.csv(do.call(rbind, samavarchi_exp_phy)[, c(1:13, 15:19)], file = "~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/GO/samavarchi_exp_phy_fdr.csv", row.names = F)

## plot
samavarchi_exp_phy_plot <- toDataFrame(samavarchi_exp_phy_df)
metaPlot(samavarchi_exp_phy_plot, "Samavarchi vs EXP PHY", 6)
ggsave("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/GO/samavarchi_fdr_exp_phy.pdf")
