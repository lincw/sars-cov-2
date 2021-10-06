# GWAS 2021b HuRI community enrichment analysis 
# Lin Chung-wen
# 06.10.2021 **13:06**--

######
# environment
gwas_dic <- "~/Documents/INET-work/virus_network/references/GWAS/"

######
# packages and functions
library(linkcomm)
library(openxlsx)
library(dplyr)
library(rstatix)

source("~/Documents/INET-work/virus_network/src/statCal.r")

gwasEnrich <- function(gwas_protein) {
	huri_ocg_GWASenrich <- do.call(rbind.data.frame, statCal(huri_ocg, gwas_protein))
	names(huri_ocg_GWASenrich) <- c("n", "p (GWAS)", "p.signif (GWAS)")
	huri_ocg_GWASenrich$cluster <- row.names(huri_ocg_GWASenrich)
	huri_ocg_GWASenrich$size <- table(huri_ocg$nodeclusters$cluster)
	huri_ocg_GWASenrich$"viral_target" <- c()
	for (i in 1:nrow(huri_ocg_GWASenrich)) {
	    nodes <- huri_ocg$nodecluster[huri_ocg$nodeclusters$cluster == i, "node"]
	    huri_ocg_GWASenrich[i, "viral_target"] <- paste(nodes[nodes %in% binary_node], collapse = ",")
	}
	huri_ocg_GWASenrich$GWAScandidate <- c()
	for (i in 1:nrow(huri_ocg_GWASenrich)) {
	    nodes <- huri_ocg$nodecluster[huri_ocg$nodeclusters$cluster == i, "node"]
	    huri_ocg_GWASenrich[i, "GWAScandidate"] <- paste(nodes[nodes %in% gwas_protein], collapse = ",")
	}

	huri_ocg_GWASoutput <- list(
	  node = huri_ocg$nodeclusters,
	  stats = huri_ocg_GWASenrich[, c(4, 2, 3, 5, 7, 6)]
	)
}

######
# load data
huri_ocg <- readRDS("~/Documents/INET-work/HuRI/HuRI_ocg.RDS")
gwas2021a <- read.csv(file.path(gwas_dic, "Genetic\ mechanisms\ of\ critical\ illness\ in\ COVID-19/table1.csv"), header = T)
gwas2021b <- read.xlsx("~/Documents/INET-work/virus_network/statistic_results/GWAS/COVID_GWAS hits_v2.xlsx") # ref: Mapping the human genetic architecture of COVID-19 (url: https://doi.org/10.1038/s41586-021-03767-x)

gwas2021a_pro <- c(gwas2021a$Locus[c(1:4)], "OAS1", "OAS2", "OAS3", gwas2021a$Locus[c(6:8)])
gwas2021a_huri <- gwas2021a_pro[gwas2021a_pro %in% V(huri_ocg$igraph)$name]
gwas2021b_pro <- gwas2021b$All.LD
gwas2021b_huri <- gwas2021b_pro[gwas2021b_pro %in% V(huri_ocg$igraph)$name]

binary <- read.xlsx('/Volumes/GoogleDrive/My\ Drive/Paper_VirHostome_CoV2/04_Supplementary\ Information/Supplementary_Table_1.xlsx', sheet = '1b - HuSCI', startRow = 4)
binary_node <- unique(binary[, "Host.protein_symbol"])

######
# community enrichment analysis
## GWAS 2021a
gwas2021a_enrich <- gwasEnrich(gwas2021a_huri)

## GWAS 2021b
gwas2021b_enrich <- gwasEnrich(gwas2021b_huri)

stats_output <- do.call(cbind, c(gwas2021a_enrich$stats, gwas2021b_enrich$stats))
names(stats_output) <- c("community", "p_gwas2021a", "p_signif_gwas2021a", "community_size", "gwas2021a_candidate", "HuSCI_viral_target", "community", "p_gwas2021b", "p_signif_gwas2021b", "community_size", "gwas2021b_candidate", "HuSCI_viral_target")
gwas2021_output <- list(
	node = gwas2021a_enrich["node"],
	stats = stats_output)

write.xlsx(gwas2021_output, file = "~/Documents/INET-work/virus_network/statistic_results/community/HuRI_communities_withGWAS2021a-b.xlsx", overwrite = TRUE)
