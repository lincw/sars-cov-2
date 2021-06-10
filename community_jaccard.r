# HuRI communities jaccard similarity analysis. The communities were used in GWAS trait analysis
# 10.06.2021
# Lin Chung-wen

#######
# load package
library(openxlsx)
library(linkcomm)

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

######
# Jaccard similarity calculation
jac_result <- c()
for (i in 1:length(community)) {
    for (j in 1:length(community)) {
        jac_result <- c(jac_result, jaccard(comm_member[[i]], comm_member[[j]]))
    }
}

jac_result_df <- matrix(jac_result, ncol = 29)
write.csv(jac_result_df, file = "/tmp/community_jaccard.csv")
