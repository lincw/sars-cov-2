# cleanup enriched communities, only consider real communities
# Lin Chung-wen
# 30.09.2021, **15:11**-- # not help a lot for statistic

######
# environment
community_path <- "~/Documents/INET-work/virus_network/statistic_results/community"
figure <- "~/Documents/INET-work/virus_network/figure_results/community"

######
# packages
library(openxlsx)
library(plotrix)

######
# date
community <- read.xlsx(file.path(community_path, "HuRI_communities_withGWAS.xlsx"), sheet = 2)
true_community <- community[community$community.size >= 4, ]
gwas_community <- true_community[, c(1:3)]
husci_community <- true_community[, c(1,2, 6)]

gwas_community$p.adjust_bonferroni <- p.adjust(gwas_community$"p.(GWAS)", method = "bonferroni")
gwas_community$p.adjust_fdr <- p.adjust(gwas_community$"p.(GWAS)", method = "BH")

husci_community$p.adjust_bonferroni <- p.adjust(husci_community$"p.(vira.target)", method = "bonferroni")
husci_community$p.adjust_fdr <- p.adjust(husci_community$"p.(vira.target)", method = "BH")

######
# visualize of HuSCI
pdf(file.path(figure, "HuSCI_enrichedCommunity_stats.pdf"), height = 5)
par(xpd = TRUE, mar = c(14, 9, 3, 2))
boxplot(husci_community[, c(5, 4, 3)], col = terrain.colors(2), main = "HuSCI community stats", horizontal = TRUE, las = 1)
mtext("p-value", side = 1, line = 2)
stripchart(husci_community[, c(5, 4, 3)], vertical = FALSE, method = "jitter", add = TRUE, pch = 20, col = "grey", jitter = 0.2)
addtable2plot(-0.1, -4, summary(husci_community[, c(3, 4, 5)], digits = 2), display.rownames = TRUE, hlines = TRUE, title = "summary")
dev.off()

# visualize of GWAS
pdf(file.path(figure, "GWAS_enrichedCommunity_stats.pdf"), height = 5)
par(xpd = TRUE, mar = c(14, 9, 3, 2))
boxplot(gwas_community[, c(5, 4, 3)], col = terrain.colors(2), main = "GWAS community stats", horizontal = TRUE, las = 1)
mtext("p-value", side = 1, line = 2)
stripchart(gwas_community[, c(5, 4, 3)], vertical = FALSE, method = "jitter", add = TRUE, pch = 20, col = "grey", jitter = 0.2)
addtable2plot(-0.1, -4, summary(gwas_community[, c(3, 4, 5)], digits = 2), display.rownames = TRUE, hlines = TRUE, title = "summary")
dev.off()
