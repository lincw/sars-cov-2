# Statistic of SARS-CoV-2 infected organ
# The null hypothesis is SARS-CoV-2 associated proteins are diversely distributed in the reported infected organs.
# Lin Chung-wen
# 2022.04.29 **13:51** --

######
# set environment
dorwie <- "~/Documents/INET-work/virus_network/toSummarize/infected_organ/"

library(openxlsx)
library(dplyr)
library(tidyr)
library(rethinking)
library(ggplot2)
library(reshape2)
library(ggpubr)

######
# load data
dataset <- getSheetNames(file.path(dorwie, "data/dorard_raw.xlsx"))
infected_organ <- lapply(dataset, function(x) {
    read.xlsx(file.path(dorwie, "data/dorard_raw.xlsx"), sheet = x)
})
names(infected_organ) <- dataset

######
# summarize data
infected_organ_sum <- lapply(infected_organ, function(x) {
    x %>% group_by(dorwie) %>% summarise(count = n()) %>% drop_na()
})
for (i in dataset) {
    names(infected_organ_sum[[i]]) <- c("dorwie", i)
}

infected_organ_df <- Reduce(function(x, y) {merge(x = x, y = y, by = "dorwie", all = T)}, 
    infected_organ_sum
    )
infected_organ_df[is.na(infected_organ_df)] <- 0

######
# randomize HPA and compared with the other dataset
prot_name <- lapply(infected_organ, function(x) unique(x$from))

rand_analysis <- function(...) {
    rand_set <- list()
    for (rand in dataset[-1]) {
        rand_pro <- sample(prot_name[["HPA"]], length(prot_name[[rand]]))
        rand_dorwie <- infected_organ[["HPA"]][infected_organ[["HPA"]]$from %in% rand_pro, ]
        rand_stat <- rand_dorwie %>% 
            group_by(dorwie) %>% 
            summarise(count = n()) %>%
            drop_na()
        names(rand_stat) <- c("dorwie", rand)
        rand_set[[rand]] <- rand_stat
    }
    rand_set_df <- Reduce(function(x, y) {merge(x = x, y = y, by = "dorwie", all = T)}, rand_set)
    names(rand_set_df) <- c("dorwie", dataset[c(2:end(dataset)[[1]])])
    rand_set_df[is.na(rand_set_df)] <- 0

    ######
    # performed t-test analysis
    output <- lapply(dataset[-1], function(x) {
        # df <- matrix(c(infected_organ_df[, x], rand_set_df[, x]), ncol = 2)
        return(t.test(infected_organ_df[, x], rand_set_df[, x], paired = T)$p.value)
        })
    return(output)
}
permutate <- mcreplicate(10000, rand_analysis(), mc.cores = detectCores())

permutate_count <- sapply(dataset[-1], function(x) {
    df <- table(permutate[x, ] > 0.05)["TRUE"]
})
names(permutate_count) <- dataset[-1]

permutate_melt <- melt(permutate)
# permutate_count_melt <- melt(permutate_count)
# permutate_count_melt$Var1 <- dataset[-1]

ggplot(permutate_melt, aes(x = value, y = ..count..)) +
    geom_histogram(bins = 10, fill = "grey") +
    geom_vline(xintercept = 0.05, color = "red") +
    facet_wrap(~ Var1, scales = "free", ncol = 4) +
    labs(x = "paired t-test p-value", y = "count", main = "Permutation distributions") +
    theme_pubr() +
    theme(strip.background = element_blank(), strip.placement = "outside")
ggsave(file.path(dorwie, "figure/p_permutate.pdf"), width = 10, height = 5)

######
# save result
save.image(file.path(dorwie, "data/p_permutate.RData"))