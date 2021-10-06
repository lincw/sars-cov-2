# viral target enrichment analysis in HuRI OCG communities
library(rstatix)
library(dplyr)
statCal <- function(ocg, target) {
        clust <- ocg$nodeclusters
        clust2 <- clust %>% select(cluster) %>% group_by(cluster) %>% mutate(count = n()) %>% distinct()
        df <- ocg$nodeclusters[ocg$nodeclusters$node %in% target, ]
        df2 <- df %>% select(cluster) %>% group_by(cluster) %>% mutate(count = n()) %>% distinct()
        merged <- merge(clust2, df2, by = "cluster", all.x = TRUE)
        merged[is.na(merged)] <- 0
        merged$notTarget <- merged[, 2] - merged[, 3]
        row.names(merged) <- merged[, 1]
        merged <- merged[, c(2, 3)]
        names(merged) <- c("size", "target")
        data <- data.frame()
        data <- apply(merged, 1, function(x) {
                tab <- as.table(matrix(c(x[2], length(target) - x[2], x[1], vcount(ocg$igraph) - x[1]), ncol = 2, byrow = TRUE))
                dimnames(tab) <- list(
                        group = c("group1", "group2"),
                        in_community = c("yes", "no"))
                fisher_test(tab)
                })
        return(data)
}
