# random test of SARS-CoV-2 targeted specific communities.
# 22.04.2022 **17:03**--
# 25.04.2022 **09:30**-- show permutation test distribution of HuSCI community size
# Lin Chung-wen

######
# load library and data
library(openxlsx)
library(dplyr)
library(rethinking)
library(ggplot2)
library(ggpubr)

cluster <- read.xlsx("~/Documents/INET-work/virus_network/statistic_results/community/HuRI_communities_withHuSCI_20220419.xlsx", sheet = 1)
comm <- read.xlsx("~/Documents/INET-work/virus_network/statistic_results/community/HuRI_communities_withHuSCI_20220419.xlsx", sheet = 2)

######
# extract community size
comm_all <- comm[comm$size >= 4, ]
original <- comm[comm$size >= 4 & comm$p < 0.05, ]
targeted <- comm[comm$size >= 4 & !is.na(comm$viral_target), ]
non_targeted <- comm[comm$size >= 4 & is.na(comm$viral_target), ]

random <- sample_n(comm_all, nrow(original))

######
# proteins from community categories (comm_all, original, targeted, non_targeted)
commAll_protein <- unique(cluster[cluster$cluster %in% comm_all$cluster, "node"])
original_protein <- unique(cluster[cluster$cluster %in% original$cluster, "node"])
targeted_protein <- unique(cluster[cluster$cluster %in% targeted$cluster, "node"])
nonTar_protein <- unique(cluster[cluster$cluster %in% non_targeted$cluster, "node"])

# intersection ratio
## original vs commAll
ori_all <- sum(original_protein %in% commAll_protein) / length(unique(c(original_protein, commAll_protein)))

## original vs targeted
ori_nonTar <- sum(original_protein %in% nonTar_protein) / length(unique(c(original_protein, nonTar_protein)))

######
# permutation test of observed communities vs randomized communities
compareRand <- function(x, y) {
    rand <- sample(y, length(x))
    rand_mean <- mean(rand)
    return(rand_mean)
}

rand_mean <- mcreplicate(10000, compareRand(original$size, comm_all$size), mc.cores = parallel::detectCores())

rand_mean_df <- data.frame(rand_mean)

gp <- ggplot(rand_mean_df, aes(x = rand_mean, y = ..count..)) + 
    geom_histogram(color = 1, fill = "lightgrey", bins = 20) +
    geom_segment(aes(x = mean(original$size), y = 200, xend = mean(original$size), yend = 10), arrow = arrow(length = unit(0.3, "cm")), color = "red", size = 1) +
    # using `stat_bin` to show the value of bins
    # stat_bin(aes(y = ..count.., label = ..count..), geom = "text", vjust = -.5, bins = 20) +
    annotate("text", label = paste0("HuSCI: ", round(mean(original$size), 2), "\np < 0.0001"), x = mean(original$size) * 0.99, y = 350, cex = 2) +
    labs(x = "Random means", y = "Count", title = "Permutation test distribution") +
    theme_pubr()
ggsave(gp, file = "~/Documents/INET-work/virus_network/figure_results/community/community_size_permutation.pdf", width = 4, height = 4)