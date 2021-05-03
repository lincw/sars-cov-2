ppi <- read.csv("~/Documents/INET-work/virus_network/pan_coronavirus/PPI_list_20210203.csv", header = T)
ppi[ppi$virus == "229E" & ppi$virus_protein == "NSP12", "human_target"]
virus <- unique(ppi$virus)[c(5, 2, 3, 4, 1, 7, 6, 8)]
human_t <- c("NSP1", "NSP3", "NSP7", "NSP13", "NSP14", "NSP16", "ORF4", "ORF4b")

jaccard_sim <- list()
jaccard_sim_name <- c()
for (v1 in 1:7) {
    for (v2 in 1:7) {
        for (i in human_t) {
            p1 <- ppi[ppi$virus == virus[v1] & ppi$virus_protein == i, "human_target"]
            p2 <- ppi[ppi$virus == virus[v2] & ppi$virus_protein == i, "human_target"]
            sim <- length(intersect(p1, p2))/length(unique(c(p1, p2)))
            jaccard_sim[paste0(i, "_", virus[v1], "-", virus[v2])] <- sim
            jaccard_sim_name <- c(jaccard_sim_name, paste0(virus[v1], "-", virus[v2]))
        }
    }
}

df <- matrix(unlist(jaccard_sim), ncol = 8, byrow = T)
colnames(df) <- human_t
rownames(df) <- unique(jaccard_sim_name)
write.csv(as.data.frame(df), file = "/tmp/df.csv")
