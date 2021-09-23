# Matthias GWAS result transform
# It should be somewhere, but I can't find the code!!!!!!!!
# Lin Chung-wen
# 21.09.2021

######
# environment
gwas_path <- "~/Documents/INET-work/virus_network/statistic_results"

######
# load package
library(tidyr)
library(openxlsx)

######
# load data
# gwas <- read.delim(file.path(gwas_path, "GWAS_Mattias/GWAS_Mattias/magma_gtex_gwas_by_community_with_gwas_info.txt"), sep = "\t")
gwas <- read.delim(file.path(gwas_path, "GWAS_Mattias/GWAS_Mattias/magma_gtex_gwas_by_community_with_gwas_info_20210921.txt"), sep = "\t")
gwas$betaSE <- gwas$BETA/gwas$SE

gwas_beta <- gwas[, c(1, 2, 5)]
gwas_p <- gwas[, c(1, 2, 8)]
gwas_fdr <- gwas[, c(1, 2, 10)]
gwas_betaSE <- gwas[, c(1, 2, 21)]

gwas_beta_wider <- pivot_wider(gwas_beta, names_from = trait, values_from = BETA)
gwas_betaSE_wider <- pivot_wider(gwas_betaSE, names_from = trait, values_from = betaSE)
gwas_p_wider <- pivot_wider(gwas_p, names_from = trait, values_from = P)
gwas_fdr_wider <- pivot_wider(gwas_fdr, names_from = trait, values_from = FDR)

toList <- list(
	BETA = gwas_beta_wider,
	BETA_SE = gwas_betaSE_wider,
	P = gwas_p_wider,
	FDR = gwas_fdr_wider
	)

######
# save result
# write.xlsx(toList, file.path(gwas_path, "GWAS_Mattias/GWAS_Mattias/magma_gtex_gwas_by_community_with_gwas_info_20201221.xlsx"), over.write = TRUE)
write.xlsx(toList, file.path(gwas_path, "GWAS_Mattias/GWAS_Mattias/magma_gtex_gwas_by_community_with_gwas_info_20210921.xlsx"), over.write = TRUE)