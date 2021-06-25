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

# for all others
bp_all <- upload_GMT_file(gmtfile = "~/workplace/database/Homo_sapiens/GO/gobp.gmt")
mf_all <- upload_GMT_file(gmtfile = "~/workplace/database/Homo_sapiens/GO/gomf.gmt")
cc_all <- upload_GMT_file(gmtfile = "~/workplace/database/Homo_sapiens/GO/gocc.gmt")

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
write.csv(husci_all_df, file = "~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/GO/husci_go_fdr.csv", row.names = F)
## plot
husci_all_plot <- toDataFrame(husci_all_df)
metaPlot(husci_all_plot, "HuSCI")
ggsave("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/GO/husci_go_fdr.pdf")

# 2. Gordon
# based on Experimental evidence codes
gordon_bp <- funcEnrich(gordon_node, bp_all, "fdr")
gordon_mf <- funcEnrich(gordon_node, mf_all, "fdr")
gordon_cc <- funcEnrich(gordon_node, cc_all, "fdr")
gordon_go <- list(BP = gordon_bp, MF = gordon_mf, CC = gordon_cc)
gordon_go_df <- do.call(rbind, gordon_go)[, c(1:13, 15:19)]
write.csv(do.call(rbind, gordon_go)[, c(1:13, 15:19)], file = "~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/GO/gordon_go_fdr.csv", row.names = F)

## plot
gordon_go_plot <- toDataFrame(gordon_go_df)
metaPlot(gordon_go_plot, "Gordon et al")
ggsave("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/GO/gordon_go_fdr.pdf")

# 3. Stukalov
# based on Experimental evidence codes
stukalov_bp <- funcEnrich(stukalov_node, bp_all, "fdr")
stukalov_mf <- funcEnrich(stukalov_node, mf_all, "fdr")
stukalov_cc <- funcEnrich(stukalov_node, cc_all, "fdr")
stukalov_go <- list(BP = stukalov_bp, MF = stukalov_mf, CC = stukalov_cc)
stukalov_go_df <- do.call(rbind, stukalov_go)[, c(1:13, 15:19)]
write.csv(do.call(rbind, stukalov_go)[, c(1:13, 15:19)], file = "~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/GO/stukalov_go_fdr.csv", row.names = F)

## plot
stukalov_go_plot <- toDataFrame(stukalov_go_df)
metaPlot(stukalov_go_plot, "Stukalov et al")
ggsave("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/GO/stukalov_go_fdr.pdf")

# 4. Li
# based on Experimental evidence codes
li_bp <- funcEnrich(li_node, bp_all, "fdr")
li_mf <- funcEnrich(li_node, mf_all, "fdr")
li_cc <- funcEnrich(li_node, cc_all, "fdr")
li_go <- list(BP = li_bp, MF = li_mf, CC = li_cc)
li_go_df <- do.call(rbind, li_go)[, c(1:13, 15:19)]
write.csv(do.call(rbind, li_go)[, c(1:13, 15:19)], file = "~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/GO/li_go_fdr.csv", row.names = F)

## plot
li_go_plot <- toDataFrame(li_go_df)
metaPlot(li_go_plot, "Li et al")
ggsave("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/GO/li_go_fdr.pdf")

# 5. Nabeel
# based on Experimental evidence codes
nabeel_bp <- funcEnrich(nabeel_node, bp_all, "fdr")
nabeel_mf <- funcEnrich(nabeel_node, mf_all, "fdr")
nabeel_cc <- funcEnrich(nabeel_node, cc_all, "fdr")
nabeel_go <- list(BP = nabeel_bp, MF = nabeel_mf, CC = nabeel_cc)
nabeel_go_df <- do.call(rbind, nabeel_go)[, c(1:13, 15:19)]
write.csv(do.call(rbind, nabeel_go)[, c(1:13, 15:19)], file = "~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/GO/nabeel_go_fdr.csv", row.names = F)

## plot
nabeel_go_plot <- toDataFrame(nabeel_go_df)
metaPlot(nabeel_go_plot, "Nabeel et al")
ggsave("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/GO/nabeel_go_fdr.pdf")

# 6. Laurent
# based on Experimental evidence codes
laurent_bp <- funcEnrich(laurent_node, bp_all, "fdr")
laurent_mf <- funcEnrich(laurent_node, mf_all, "fdr")
laurent_cc <- funcEnrich(laurent_node, cc_all, "fdr")
laurent_go <- list(BP = laurent_bp, MF = laurent_mf, CC = laurent_cc)
laurent_go_df <- do.call(rbind, laurent_go)[, c(1:13, 15:19)]
write.csv(do.call(rbind, laurent_go)[, c(1:13, 15:19)], file = "~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/GO/laurent_go_fdr.csv", row.names = F)

## plot
laurent_go_plot <- toDataFrame(laurent_go_df)
metaPlot(laurent_go_plot, "Laurent et al")
ggsave("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/GO/laurent_go_fdr.pdf")

# 7. St-Germain
# based on Experimental evidence codes
stgermain_bp <- funcEnrich(stgermain_node, bp_all, "fdr")
stgermain_mf <- funcEnrich(stgermain_node, mf_all, "fdr")
stgermain_cc <- funcEnrich(stgermain_node, cc_all, "fdr")
stgermain_go <- list(BP = stgermain_bp, MF = stgermain_mf, CC = stgermain_cc)
stgermain_go_df <- do.call(rbind, stgermain_go)[, c(1:13, 15:19)]
write.csv(do.call(rbind, stgermain_go)[, c(1:13, 15:19)], file = "~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/GO/stgermain_go_fdr.csv", row.names = F)

## plot
stgermain_go_plot <- toDataFrame(stgermain_go_df)
metaPlot(stgermain_go_plot, "St-Germain et al")
ggsave("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/GO/stgermain_go_fdr.pdf")

# 8. Samavarchi
# based on Experimental evidence codes
samavarchi_bp <- funcEnrich(samavarchi_node, bp_all, "fdr")
samavarchi_mf <- funcEnrich(samavarchi_node, mf_all, "fdr")
samavarchi_cc <- funcEnrich(samavarchi_node, cc_all, "fdr")
samavarchi_go <- list(BP = samavarchi_bp, MF = samavarchi_mf, CC = samavarchi_cc)
samavarchi_go_df <- do.call(rbind, samavarchi_go)[, c(1:13, 15:19)]
write.csv(do.call(rbind, samavarchi_go)[, c(1:13, 15:19)], file = "~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/GO/samavarchi_go_fdr.csv", row.names = F)

## plot
samavarchi_go_plot <- toDataFrame(samavarchi_go_df)
metaPlot(samavarchi_go_plot, "Samavarchi et al")
ggsave("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/GO/samavarchi_go_fdr.pdf")
