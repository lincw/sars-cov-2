# HuRI community viral target enrichment, HuSCI, Gordon et al and Stukalov et al
# Lin Chung-wen
# 19.08.2021

######
# load packages & functions
library(linkcomm)
library(dplyr)
library(openxlsx)
source("~/Documents/INET-work/virus_network/src/statCal.r")

# load date
load("~/Documents/INET-work/virus_network/statistic_results/community/HuRI_ocg.RData")
gordon <- read.xlsx("/Volumes/GoogleDrive/My\ Drive/VirHostome_CW/GitHub/data/extended_table/Extended_Table_2_PPIs.xlsx", sheet = "Gordon")
stukalov <- read.xlsx("/Volumes/GoogleDrive/My\ Drive/VirHostome_CW/GitHub/data/extended_table/Extended_Table_2_PPIs.xlsx", sheet = "Stukalov")

######
# community viral target enrichment analyses
husci_enrich <- do.call(
    rbind.data.frame,
    statCal(huri_ocg, length(binary_node), length(unique(huri_ocg$nodeclusters$node)), binary_node)
)
husci_enrich$cluster <- row.names(husci_enrich)

gordon_list <- unique(gordon$PreyGene)
gordon_enrich <- do.call(
    rbind.data.frame,
    statCal(huri_ocg, length(gordon_list), length(unique(huri_ocg$nodeclusters$node)), gordon_list)
)
gordon_enrich$cluster <- row.names(gordon_enrich)

stukalov_list <- unique(stukalov$human)
stukalov_enrich <- do.call(
    rbind.data.frame,
    statCal(huri_ocg, length(stukalov_list), length(unique(huri_ocg$nodeclusters$node)), stukalov_list)
)
stukalov_enrich$cluster <- row.names(stukalov_enrich)

######
# save image
save.image("~/Documents/INET-work/virus_network/statistic_results/community/3dataset_enrichment.RData")
