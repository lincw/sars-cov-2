# network profile calculation
# Lin Chung-wen
countRandomize <- function(node) {
    count <- c()
    gwas_rand1 <- randomGwas(node)
    count[1] <- gsize(gwas_rand1) # ecount
    count[2] <- vcount(gwas_rand1) # vcount
    gwas_rand1_strong <- components(gwas_rand1, mode = "strong")
    biggest_cluster_id <- which.max(gwas_rand1_strong$csize)
    vert_ids <- V(gwas_rand1)[gwas_rand1_strong$membership == biggest_cluster_id]
    gwas_rand1_strong_graph <- induced_subgraph(gwas_rand1, vert_ids)
    count[3] <- gsize(gwas_rand1_strong_graph) # str_ecount
    count[4] <- vcount(gwas_rand1_strong_graph) # str_vcount
    gwas_rand1_eb <- cluster_edge_betweenness(gwas_rand1_strong_graph, directed = FALSE)
    count[5] <- length(table(gwas_rand1_eb$membership)) # eb_size
    gwas_rand1_ocg <- getOCG.clusters(as_data_frame(gwas_rand1), init.class.sys = 1, cent.class.sys = 0, verbose = FALSE)
    count[6] <- length(gwas_rand1_ocg$clustsizes) # ocg_size
    return(count)
}
# count community size difference between largest one and second
countRandomize2 <- function(node) {
    count <- data.frame()
    gwas_rand1 <- randomGwas(node)
    gwas_rand1_strong <- components(gwas_rand1, mode = "strong")
    biggest_cluster_id <- which.max(gwas_rand1_strong$csize)
    vert_ids <- V(gwas_rand1)[gwas_rand1_strong$membership == biggest_cluster_id]
    gwas_rand1_strong_graph <- induced_subgraph(gwas_rand1, vert_ids)
    gwas_rand1_eb <- cluster_edge_betweenness(gwas_rand1_strong_graph, directed = FALSE)
    gwas_rand1_eb_diff <- sort(sizes(gwas_rand1_eb), decreasing = TRUE)[c(1, 2)]
    count[1, "eb"] <- gwas_rand1_eb_diff[1] - gwas_rand1_eb_diff[2]
    gwas_rand1_ocg <- getOCG.clusters(as_data_frame(gwas_rand1), init.class.sys = 1, cent.class.sys = 0, verbose = FALSE)
    gwas_rand1_ocg_diff <- sort(gwas_rand1_ocg$clustsizes, decreasing = TRUE)[c(1, 2)]
    count[1, "ocg"] <- gwas_rand1_ocg_diff[1] - gwas_rand1_ocg_diff[2]
    return(count)
}
