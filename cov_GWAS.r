# profile of GWAS result from Mattias
## date: 2021.04.12

library(dplyr)

gwas <- read.delim("~/Documents/INET-work/virus_network/statistic_results/GWAS_Mattias/magma_gtex_gwas_by_community_with_gwas_info.txt", header = T, sep = "\t")
gwas %>% group_by(VARIABLE) %>% summarise(paste0(Phenotype, collapse = ",")) %>% write.table("/tmp/HPA.tsv", sep = "\t", row.names = F)
