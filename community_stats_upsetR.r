# upset plot for HuRI communities with either HuSCI or AP-MS data enrichment
# Lin Chung-wen
# 05.05.2022 **15:01** --
# 05.05.2022 **23:13* -- update the APMS_sign_community.xlsx with critical illness associated information

library(openxlsx)
library(UpSetR)
library(dplyr)
community <- "~/Documents/INET-work/virus_network/statistic_results/community/"

sign <- read.xlsx(file.path(community, "APMS_sign_community.xlsx"), sheet = 1)

sign_df <- sign %>% group_by(cluster, husci_enriched, gordon_enriched, stukalov_enriched, li_enriched, nabeel_enriched) %>% summarise(count = n())

sign_list <- list(
    'Husci' = sign_df[sign_df$husci_enriched == 1, 'cluster'],
    'Gordon et al.' = sign_df[sign_df$gordon_enriched == 1, 'cluster'],
    'Stukalov et al.' = sign_df[sign_df$stukalov_enriched == 1, 'cluster'],
    'Li et al.' = sign_df[sign_df$li_enriched == 1, 'cluster'],
    'Nabeel et al.' = sign_df[sign_df$nabeel_enriched == 1, 'cluster']
)
sign_list <- lapply(sign_list, function(x) as.data.frame(x)$cluster)

pdf("~/Documents/INET-work/virus_network/figure_results/community/sig_comm_upsetR.pdf", width = 7, height = 4)
upset(fromList(sign_list), nintersects = NA,
    point.size = 2, line.size = 0.5, text.scale = 1,
    empty.intersections = "on",
    queries = list(list(query = intersects, params = list(names(sign_list)), color = "blue", active = T)),
    sets.bar.color = c("darkred", rep("gray23", 4)))
dev.off()

######
# 05.05.2022 **23:13**
sign_critical <- sign[sign$GWAS_critical_illness_Nature2021a != 0, ]

sign_critical_df <- sign_critical %>% group_by(cluster, husci_enriched, gordon_enriched, stukalov_enriched, li_enriched, nabeel_enriched) %>% summarise(count = n())

sign_critical_list <- list(
    'HuSCI' = sign_critical_df[sign_critical_df$husci_enriched == 1, 'cluster'],
    'Gordon et al.' = sign_critical_df[sign_critical_df$gordon_enriched == 1, 'cluster'],
    'Stukalov et al.' = sign_critical_df[sign_critical_df$stukalov_enriched == 1, 'cluster'],
    'Li et al.' = sign_critical_df[sign_critical_df$li_enriched == 1, 'cluster'],
    'Nabeel et al.' = sign_critical_df[sign_critical_df$nabeel_enriched == 1, 'cluster']
)
sign_critical_list <- lapply(sign_critical_list, function(x) as.data.frame(x)$cluster)

pdf("~/Documents/INET-work/virus_network/figure_results/community/sig_comm_critical_upsetR.pdf", width = 7, height = 4)
upset(fromList(sign_critical_list), sets = c("HuSCI", "Gordon et al.", "Stukalov et al.", "Li et al.", "Nabeel et al."),
    point.size = 2, line.size = 0.5, text.scale = 1,
    empty.intersections = "on")
dev.off()
