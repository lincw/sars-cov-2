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
funcEnrich <- function(query_list, organism, correction) {
    goquery <- gost(query = query_list, organism = organism, correction_method = correction, evcodes = TRUE)
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
metaPlot <- function(x, main) {
    x$term <- factor(x$term, levels = x$term[order(x$observeRatio)])
    ggplot(x, aes(x = term, y = observeRatio)) +
        geom_bar(stat = "identity", fill = rgb(128, 41, 227, maxColorValue = 255)) +
        coord_flip() +
        labs(y = "Effect size (%)", x = "Functional term (p < 0.05)", title = main) +
        theme_bw() +
        theme(axis.text = element_text(color = "black", size = 8, hjust = 0),
        axis.title = element_text(size = 8, color = "black", face = "bold"))
}

######
# load different GO background
# for HuSCI
horf_bp_all <- upload_GMT_file(gmtfile = "~/workplace/database/Homo_sapiens/GO/hORFeome_gobp.gmt")
horf_mf_all <- upload_GMT_file(gmtfile = "~/workplace/database/Homo_sapiens/GO/hORFeome_gomf.gmt")
horf_cc_all <- upload_GMT_file(gmtfile = "~/workplace/database/Homo_sapiens/GO/hORFeome_gocc.gmt")

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
husci_bp <- funcEnrich(husci_node, horf_bp_all, "bonferroni")
husci_mf <- funcEnrich(husci_node, horf_mf_all, "bonferroni")
husci_cc <- funcEnrich(husci_node, horf_cc_all, "bonferroni")
husci_all <- list(BP = husci_bp, MF = husci_mf, CC = husci_cc)
husci_all_df <- do.call(rbind, husci_all)[, c(1:13, 15:19)]
write.csv(husci_all_df, file = "~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/GO/husci_all_bonferroni.csv", row.names = F)
## plot
husci_all_plot <- toDataFrame(husci_all_df)
metaPlot(husci_all_plot, "HuSCI vs all GO evidence codes")
ggsave("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/GO/husci_all.pdf")

# based on Experimental evidence codes
husci_bp_exp <- funcEnrich(husci_node, horf_bp_exp, "bonferroni")
husci_mf_exp <- funcEnrich(husci_node, horf_mf_exp, "bonferroni")
husci_cc_exp <- funcEnrich(husci_node, horf_cc_exp, "bonferroni") # no significant terms
husci_exp <- list(BP = husci_bp_exp, MF = husci_mf_exp)
husci_exp_df <- do.call(rbind, husci_exp)[, c(1:13, 15:19)]
write.csv(do.call(rbind, husci_exp)[, c(1:13, 15:19)], file = "~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/GO/husci_exp_bonferroni.csv", row.names = F)

## plot
husci_exp_plot <- toDataFrame(husci_exp_df)
metaPlot(husci_exp_plot, "HuSCI vs GO EXP codes")
ggsave("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/GO/husci_exp.pdf")

# based on Experimental evidence codes and phylogenetically-inferred  annotations
husci_bp_exp_phy <- funcEnrich(husci_node, horf_bp_exp_phy, "bonferroni")
husci_mf_exp_phy <- funcEnrich(husci_node, horf_mf_exp_phy, "bonferroni")
husci_cc_exp_phy <- funcEnrich(husci_node, horf_cc_exp_phy, "bonferroni")
husci_exp_phy <- list(BP = husci_bp_exp_phy, MF = husci_mf_exp_phy, CC = husci_cc_exp_phy)
husci_exp_phy_df <- do.call(rbind, husci_exp_phy)[, c(1:13, 15:19)]
write.csv(do.call(rbind, husci_exp_phy)[, c(1:13, 15:19)], file = "~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/GO/husci_exp_phy_bonferroni.csv", row.names = F)

## plot
husci_exp_phy_plot <- toDataFrame(husci_exp_phy_df)
metaPlot(husci_exp_phy_plot, "HuSCI vs GO EXP PHY codes")
ggsave("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/GO/husci_exp_phy.pdf")

# 2. Gordon
# based on Experimental evidence codes and phylogenetically-inferred  annotations
gordon_bp_exp_phy <- funcEnrich(gordon_node, bp_exp_phy, "bonferroni")
gordon_mf_exp_phy <- funcEnrich(gordon_node, mf_exp_phy, "bonferroni")
gordon_cc_exp_phy <- funcEnrich(gordon_node, cc_exp_phy, "bonferroni")
gordon_exp_phy <- list(BP = gordon_bp_exp_phy, MF = gordon_mf_exp_phy, CC = gordon_cc_exp_phy)
gordon_exp_phy_df <- do.call(rbind, gordon_exp_phy)[, c(1:13, 15:19)]
# write.csv(do.call(rbind, gordon_exp_phy)[, c(1:13, 15:19)], file = "~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/GO/gordon_exp_phy_bonferroni.csv", row.names = F)

## plot
gordon_exp_phy_plot <- toDataFrame(gordon_exp_phy_df)
metaPlot(gordon_exp_phy_plot, "gordon vs GO EXP PHY codes")
ggsave("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/GO/gordon_exp_phy.pdf")

# 3. Stukalov
# based on Experimental evidence codes and phylogenetically-inferred  annotations
stukalov_bp_exp_phy <- funcEnrich(stukalov_node, bp_exp_phy, "bonferroni")
stukalov_mf_exp_phy <- funcEnrich(stukalov_node, mf_exp_phy, "bonferroni")
stukalov_cc_exp_phy <- funcEnrich(stukalov_node, cc_exp_phy, "bonferroni")
stukalov_exp_phy <- list(BP = stukalov_bp_exp_phy, MF = stukalov_mf_exp_phy, CC = stukalov_cc_exp_phy)
stukalov_exp_phy_df <- do.call(rbind, stukalov_exp_phy)[, c(1:13, 15:19)]
# write.csv(do.call(rbind, stukalov_exp_phy)[, c(1:13, 15:19)], file = "~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/GO/stukalov_exp_phy_bonferroni.csv", row.names = F)

## plot
stukalov_exp_phy_plot <- toDataFrame(stukalov_exp_phy_df)
metaPlot(stukalov_exp_phy_plot, "stukalov vs GO EXP PHY codes")
ggsave("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/GO/stukalov_exp_phy.pdf")

# 4. Li
# based on Experimental evidence codes and phylogenetically-inferred  annotations
li_bp_exp_phy <- funcEnrich(li_node, bp_exp_phy, "bonferroni")
li_mf_exp_phy <- funcEnrich(li_node, mf_exp_phy, "bonferroni")
li_cc_exp_phy <- funcEnrich(li_node, cc_exp_phy, "bonferroni")
li_exp_phy <- list(BP = li_bp_exp_phy, MF = li_mf_exp_phy, CC = li_cc_exp_phy)
li_exp_phy_df <- do.call(rbind, li_exp_phy)[, c(1:13, 15:19)]
# write.csv(do.call(rbind, li_exp_phy)[, c(1:13, 15:19)], file = "~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/GO/li_exp_phy_bonferroni.csv", row.names = F)

## plot
li_exp_phy_plot <- toDataFrame(li_exp_phy_df)
metaPlot(li_exp_phy_plot, "li vs GO EXP PHY codes")
ggsave("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/GO/li_exp_phy.pdf")

# 5. Nabeel
# based on Experimental evidence codes and phylogenetically-inferred  annotations
nabeel_bp_exp_phy <- funcEnrich(nabeel_node, bp_exp_phy, "bonferroni")
nabeel_mf_exp_phy <- funcEnrich(nabeel_node, mf_exp_phy, "bonferroni")
nabeel_cc_exp_phy <- funcEnrich(nabeel_node, cc_exp_phy, "bonferroni")
nabeel_exp_phy <- list(BP = nabeel_bp_exp_phy, MF = nabeel_mf_exp_phy, CC = nabeel_cc_exp_phy)
nabeel_exp_phy_df <- do.call(rbind, nabeel_exp_phy)[, c(1:13, 15:19)]
# write.csv(do.call(rbind, nabeel_exp_phy)[, c(1:13, 15:19)], file = "~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/GO/nabeel_exp_phy_bonferroni.csv", row.names = F)

## plot
nabeel_exp_phy_plot <- toDataFrame(nabeel_exp_phy_df)
metaPlot(nabeel_exp_phy_plot, "nabeel vs GO EXP PHY codes")
ggsave("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/GO/nabeel_exp_phy.pdf")

# 6. Laurent
# based on Experimental evidence codes and phylogenetically-inferred  annotations
laurent_bp_exp_phy <- funcEnrich(laurent_node, bp_exp_phy, "bonferroni")
laurent_mf_exp_phy <- funcEnrich(laurent_node, mf_exp_phy, "bonferroni")
laurent_cc_exp_phy <- funcEnrich(laurent_node, cc_exp_phy, "bonferroni")
laurent_exp_phy <- list(BP = laurent_bp_exp_phy, MF = laurent_mf_exp_phy, CC = laurent_cc_exp_phy)
laurent_exp_phy_df <- do.call(rbind, laurent_exp_phy)[, c(1:13, 15:19)]
# write.csv(do.call(rbind, laurent_exp_phy)[, c(1:13, 15:19)], file = "~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/GO/laurent_exp_phy_bonferroni.csv", row.names = F)

## plot
laurent_exp_phy_plot <- toDataFrame(laurent_exp_phy_df)
metaPlot(laurent_exp_phy_plot, "laurent vs GO EXP PHY codes")
ggsave("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/GO/laurent_exp_phy.pdf")

# 7. St-Germain
# based on Experimental evidence codes and phylogenetically-inferred  annotations
stgermain_bp_exp_phy <- funcEnrich(stgermain_node, bp_exp_phy, "bonferroni")
stgermain_mf_exp_phy <- funcEnrich(stgermain_node, mf_exp_phy, "bonferroni")
stgermain_cc_exp_phy <- funcEnrich(stgermain_node, cc_exp_phy, "bonferroni")
stgermain_exp_phy <- list(BP = stgermain_bp_exp_phy, MF = stgermain_mf_exp_phy, CC = stgermain_cc_exp_phy)
stgermain_exp_phy_df <- do.call(rbind, stgermain_exp_phy)[, c(1:13, 15:19)]
# write.csv(do.call(rbind, stgermain_exp_phy)[, c(1:13, 15:19)], file = "~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/GO/stgermain_exp_phy_bonferroni.csv", row.names = F)

## plot
stgermain_exp_phy_plot <- toDataFrame(stgermain_exp_phy_df)
metaPlot(stgermain_exp_phy_plot, "St-Germain vs GO EXP PHY codes")
ggsave("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/GO/stgermain_exp_phy.pdf")

# 8. Samavarchi
# based on Experimental evidence codes and phylogenetically-inferred  annotations
samavarchi_bp_exp_phy <- funcEnrich(samavarchi_node, bp_exp_phy, "bonferroni")
samavarchi_mf_exp_phy <- funcEnrich(samavarchi_node, mf_exp_phy, "bonferroni")
samavarchi_cc_exp_phy <- funcEnrich(samavarchi_node, cc_exp_phy, "bonferroni")
samavarchi_exp_phy <- list(BP = samavarchi_bp_exp_phy, MF = samavarchi_mf_exp_phy, CC = samavarchi_cc_exp_phy)
samavarchi_exp_phy_df <- do.call(rbind, samavarchi_exp_phy)[, c(1:13, 15:19)]
# write.csv(do.call(rbind, samavarchi_exp_phy)[, c(1:13, 15:19)], file = "~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/GO/samavarchi_exp_phy_bonferroni.csv", row.names = F)

## plot
samavarchi_exp_phy_plot <- toDataFrame(samavarchi_exp_phy_df)
metaPlot(samavarchi_exp_phy_plot, "Samavarchi vs GO EXP PHY codes")
ggsave("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/GO/samavarchi_exp_phy.pdf")
