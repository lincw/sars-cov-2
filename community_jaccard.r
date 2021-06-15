# HuRI communities jaccard similarity analysis. The communities were used in GWAS trait analysis
# 10.06.2021
# Lin Chung-wen

#######
# load package
library(openxlsx)
library(linkcomm)
library(pheatmap)

######
# load OCG data
load("~/Documents/INET-work/virus_network/statistic_results/HuRI_ocg.RData")

######
# function, define Jaccard Similarity function
jaccard <- function(a, b) {
    intersection <- length(intersect(a, b))
    union <- length(a) + length(b) - intersection
    return (intersection/union)
}

######
# The communities used in GWAS trait analysis
community <- c(8, 145, 316, 364, 377, 415, 525, 571, 692, 731, 891, 926, 1353, 1515, 1627, 1652, 1833, 1882, 2398, 2545, 2769, 2831, 3035, 3205, 3564, 3682, 3941, 4160, 4227)

comm_member <- list()
for (i in 1:length(community)) {
    comm_member[[i]] <- getNodesIn(huri_ocg, clusterids = community[i])
}
write.xlsx(unlist(lapply(comm_member, function(x) paste0(x, collapse = ","))), file = "/tmp/gwas_trait.xlsx")

######
# Jaccard similarity calculation
jac_result <- c()
for (i in 1:length(community)) {
    for (j in 1:length(community)) {
        jac_result <- c(jac_result, jaccard(comm_member[[i]], comm_member[[j]]))
    }
}

jac_result_df <- matrix(jac_result, ncol = 29)
rownames(jac_result_df) <- community
colnames(jac_result_df) <- community
pdf("~/Documents/INET-work/virus_network/figure_results/gwas_trait_jaccard.pdf", width = 5, height = 5)
pheatmap(jac_result_df, cutree_rows = 8, cutree_cols = 8, main = "Jaccard similarity of community membership",
    legend_breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
    legend_labels = c("0", "0.2", "0.4", "0.6", "0.8", "score\n"),
    legend = TRUE)
dev.off()

df <- kmeans(jac_result_df, 8)

write.csv(jac_result_df, file = "/tmp/community_jaccard.csv")
