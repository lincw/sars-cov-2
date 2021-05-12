# Generation HuSCI randomized network, with degree preservation
# author: Lin Chung-wen
# date: 2021.05.11

library(igraph)
library(progress)

binary <- read.csv("../data/HuSCI_edge.csv", header = T)

graph <- graph_from_data_frame(binary, directed = FALSE)

pb <- progress_bar$new(total = 10000)

for (i in 1:10000) {
    pb$tick()
    random_g <- rewire(graph, keeping_degseq(niter = gsize(graph) * 10))
    write.csv(as_data_frame(random_g), file = paste0("../data/HuSCI_randomized_network/HuSCI_randomized_network_", i, ".csv"), row.names = FALSE)
}
