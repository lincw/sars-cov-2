# have network graph of SARS-CoV-2-host interactome
# date: 20201013
# author: LCW

library(igraph)
library(tidyr)
library(UpSetR)
###################################
# INET SARS-CoV-2 interactome
node_inet <- read.csv("~/Documents/INET-work/virus_network/Y2H_screening/20200923_verified/node_20200924.csv", header = T, as.is = T)
edge_inet <- read.csv("~/Documents/INET-work/virus_network/Y2H_screening/20200923_verified/edge_20200924.csv", header = T, as.is = T)

###################################
# Roth SARS-CoV-2 interactome
node_roth <- read.csv("~/Documents/INET-work/virus_network/Y2H_screening/source_data/dk_node.csv", header = T, as.is = T)
edge_roth <- read.csv("~/Documents/INET-work/virus_network/Y2H_screening/source_data/dk_BFG.csv", header = T, as.is = T)

###################################
# merge two interactome from INET and Roth
node_inet <- node_inet[, c(1, 5)]
node <- unique(rbind(node_inet, node_roth))
names(edge_inet) <- c("human", "virus")
edge_roth <- edge_roth[, c(1, 2)]
names(edge_roth) <- c("virus", "human")
edge <- rbind(edge_inet, edge_roth)
net <- graph_from_data_frame(d = edge, vertices = node, directed = FALSE)
V(net)$species <- node$species
V(net)$color <- ifelse(V(net)$species == "human", "blue", "pink")

############################################
# plotting
coords <- layout_(net, with_lgl())
pdf("merged_virhostome.pdf", width = 15, height = 15)
plot(net, vertex.size = 3, vertex.label.cex = .9, vertex.label.dist = 0.5, vertex.label.degree = 1, vertex.label.color = "black", margin = 0.1)
dev.off()
