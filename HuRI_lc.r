library(igraph)
library(linkcomm)
library(openxlsx)

# original
# huri <- read.csv("HuRI_Tong_withSymbol.csv", header = T)
# huri <- huri[, c(5:6)]
# huri_lc <- getLinkCommunities(huri, use.all.edges = TRUE, plot = FALSE)
# save.image(file = "HuRI_linkCommunities.RData")

# **2020.11.23 22:30** use individual assay, with LinkComm default parameters
huri_assay1 <- read.xlsx("~/Documents/INET-work/references/HuRI_binaryPPI/HuRI_individual_assay.xlsx", sheet = 1)
huri_assay2 <- read.xlsx("~/Documents/INET-work/references/HuRI_binaryPPI/HuRI_individual_assay.xlsx", sheet = 2)
huri_assay3 <- read.xlsx("~/Documents/INET-work/references/HuRI_binaryPPI/HuRI_individual_assay.xlsx", sheet = 3)

huri_lc1 <- getLinkCommunities(huri_assay1[, c(3, 4)], plot = FALSE)
huri_lc2 <- getLinkCommunities(huri_assay2[, c(3, 4)], plot = FALSE)
huri_lc3 <- getLinkCommunities(huri_assay3[, c(3, 4)], plot = FALSE)
save.image(file = "~/Documents/INET-work/virus_network/statistic_results/HuRI_linkcomm_simple_assays.RData")
