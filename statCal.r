# viral target enrichment analysis in HuRI OCG communities
library(rstatix)
statCal <- function(ocg, virTarget, totalN, viral_targets = binary_node) {
        clust <- ocg$nodeclusters
        clust2 <- clust %>% select(cluster) %>% group_by(cluster) %>% mutate(count = n()) %>% distinct()
        df <- ocg$nodeclusters[ocg$nodeclusters$node %in% viral_targets, ]
        df2 <- df %>% select(cluster) %>% group_by(cluster) %>% mutate(count = n()) %>% distinct()
        merged <- merge(clust2, df2, by = "cluster", all.x = TRUE)
        merged[is.na(merged)] <- 0
        merged$notTarget <- merged[, 2] - merged[, 3]
        # merged$cluster <- merged$cluster
        row.names(merged) <- merged[, 1]
        merged <- merged[, c(3, 4)]
        names(merged) <- c("Target", "notTarget")
        total <- c(virTarget, totalN - virTarget)
        names(total) <- c("Target", "notTarget")
        data <- data.frame()
        data <- apply(merged, 1, function(x) {
                tab <- rbind(x, total)
                fisher_test(tab)
                })
        return(data)
}
