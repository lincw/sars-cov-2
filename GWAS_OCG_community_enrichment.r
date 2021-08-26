# HuRI community viral target enrichment, HuSCI, Gordon et al and Stukalov et al, GWAS proteins
# Lin Chung-wen
# 19.08.2021

######
# load packages & functions
library(linkcomm)
library(dplyr)
library(openxlsx)
source("~/Documents/INET-work/virus_network/src/statCal.r")

toPlot <- function(enrich, interactor, file) {
    clusters <- enrich %>% filter(clustersig == "*") %>% select(cluster)
    pdf(file, width = 5, height = 5)
    par(xpd = FALSE)
    for (i in clusters$cluster) {
        el <- huri_ocg$edgelist[getEdgesIn(huri_ocg, clusterids = i), ]
        ig <- graph.edgelist(el, directed = FALSE)
        V(ig)$color <- ifelse(V(ig)$name %in% interactor, "red", "blue")
        plot(ig, vertex.label.dist = -3, vertex.label.color = "black", main = paste0("cluster = ", i), bty = "L")
        legend(0.6, 1.5, legend = c("viral target", "human protein"), pch = 21, pt.bg = c("red", "blue"), bty = "n")
    }
    dev.off()
}
# load date
load("~/Documents/INET-work/virus_network/statistic_results/community/HuRI_ocg.RData")
gordon <- read.xlsx("/Volumes/GoogleDrive/My\ Drive/VirHostome_CW/GitHub/data/extended_table/Extended_Table_2_PPIs.xlsx", sheet = "Gordon")
stukalov <- read.xlsx("/Volumes/GoogleDrive/My\ Drive/VirHostome_CW/GitHub/data/extended_table/Extended_Table_2_PPIs.xlsx", sheet = "Stukalov")
gwas_2021a <- read.csv("~/Documents/INET-work/virus_network/references/GWAS/Genetic\ mechanisms\ of\ critical\ illness\ in\ COVID-19/table1.csv", header = T)
gwas_2021b <- read.xlsx("~/Documents/INET-work/virus_network/statistic_results/GWAS/COVID_GWAS_hit_inHUSCI_v2.xlsx")

######
# GWAS list
gwas2021a <- c(gwas_2021a$Locus[c(1:4, 6:8)], "OAS1", "OAS2", "OAS3")
gwas2021b <- gwas_2021b$All.LD

######
# community viral target enrichment analyses
husci_enrich <- do.call(
    rbind.data.frame,
    statCal(huri_ocg, length(binary_node), length(unique(huri_ocg$nodeclusters$node)), binary_node)
)
husci_enrich$cluster <- row.names(husci_enrich)
husci_enrich$fdr <- p.adjust(husci_enrich$p, method = "fdr")
husci_enrich$clustersize <- huri_ocg$clustsizes[husci_enrich$cluster]
husci_enrich$clustersig <- ifelse(husci_enrich$p < 0.05 & husci_enrich$clustersize >= 4, "*", NA)

gwas2021a_enrich <- do.call(
    rbind.data.frame,
    statCal(huri_ocg, length(gwas2021a), length(unique(huri_ocg$nodeclusters$node)), gwas2021a)
    )
gwas2021a_enrich$cluster <- row.names(gwas2021a_enrich)
gwas2021a_enrich$fdr <- p.adjust(gwas2021a_enrich$p, method = "fdr")
gwas2021a_enrich$clustersize <- huri_ocg$clustsizes[gwas2021a_enrich$cluster]
gwas2021a_enrich$clustersig <- ifelse(gwas2021a_enrich$p < 0.05 & gwas2021a_enrich$clustersize >= 4, "*", NA)

gwas2021b_enrich <- do.call(
    rbind.data.frame,
    statCal(huri_ocg, length(gwas2021b), length(unique(huri_ocg$nodeclusters$node)), gwas2021b)
    )
gwas2021b_enrich$cluster <- row.names(gwas2021b_enrich)
gwas2021b_enrich$fdr <- p.adjust(gwas2021b_enrich$p, method = "fdr")
gwas2021b_enrich$clustersize <- huri_ocg$clustsizes[gwas2021b_enrich$cluster]
gwas2021b_enrich$clustersig <- ifelse(gwas2021b_enrich$p < 0.05 & gwas2021b_enrich$clustersize >= 4, "*", NA)

# redundant
gordon_list <- unique(gordon$PreyGene)
gordon_enrich <- do.call(
    rbind.data.frame,
    statCal(huri_ocg, length(gordon_list), length(unique(huri_ocg$nodeclusters$node)), gordon_list)
)
gordon_enrich$cluster <- row.names(gordon_enrich)
gordon_enrich$fdr <- p.adjust(gordon_enrich$p, method = "fdr")
gordon_enrich$clustersize <- huri_ocg$clustsizes[gordon_enrich$cluster]
gordon_enrich$clustersig <- ifelse(gordon_enrich$p < 0.05 & gordon_enrich$clustersize >= 4, "*", NA)

stukalov_list <- unique(stukalov$human)
stukalov_enrich <- do.call(
    rbind.data.frame,
    statCal(huri_ocg, length(stukalov_list), length(unique(huri_ocg$nodeclusters$node)), stukalov_list)
)
stukalov_enrich$cluster <- row.names(stukalov_enrich)
stukalov_enrich$fdr <- p.adjust(stukalov_enrich$p, method = "fdr")
stukalov_enrich$clustersize <- huri_ocg$clustsizes[stukalov_enrich$cluster]
stukalov_enrich$clustersig <- ifelse(stukalov_enrich$p < 0.05 & stukalov_enrich$clustersize >= 4, "*", NA)

######
# plotting
toPlot(husci_enrich, binary_node, "~/Documents/INET-work/virus_network/figure_results/community/husci_communities.pdf")

toPlot(gordon_enrich, gordon$PreyGene, "~/Documents/INET-work/virus_network/figure_results/community/gordon_communities.pdf")

toPlot(stukalov_enrich, stukalov$human, "~/Documents/INET-work/virus_network/figure_results/community/stukalov_communities.pdf")
toPlot(gwas2021a_enrich, gwas2021a, "~/Documents/INET-work/virus_network/figure_results/community/gwas2021a_communities.pdf")
toPlot(gwas2021b_enrich, gwas2021b, "~/Documents/INET-work/virus_network/figure_results/community/gwas2021b_communities.pdf")

######
# save results
toXLSX <- list(HuSCI = husci_enrich[, c(4, 2, 5:7)], Gordon = gordon_enrich[, c(4, 2, 5:7)], Stukalov = stukalov_enrich[, c(4, 2, 5:7)], GWAS_2021a = gwas2021a_enrich[, c(4, 2, 5:7)], GWAS_2021b = gwas2021b_enrich[, c(4, 2, 5:7)])
write.xlsx(toXLSX, file = "~/Documents/INET-work/virus_network/statistic_results/community/3dataset_enrichment.xlsx", overwrite = TRUE)

######
# save image
save.image("~/Documents/INET-work/virus_network/statistic_results/community/3dataset_enrichment.RData")
