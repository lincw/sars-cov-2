###################################
# have network graph and statistic of SARS-CoV-2-host interactome
# date: 20201008
# author: LCW

library(igraph)
library(tidyr)
library(UpSetR)
###################################
# INET SARS-CoV-2 interactome
node_inet <- read.csv("~/Documents/INET-work/virus_network/Y2H_screening/20200923_verified/node_20200924.csv", header = T, as.is = T)
edge_inet <- read.csv("~/Documents/INET-work/virus_network/Y2H_screening/20200923_verified/edge_20200924.csv", header = T, as.is = T)
net_inet <- graph_from_data_frame(d = edge_inet, vertices = node_inet, directed = F)
V(net_inet)$species <- node_inet$species
V(net_inet)$color <- ifelse(V(net_inet)$species == "human", "blue", "pink")
V(net_inet)[1]; E(net_inet)[1]; table(node_inet$species)[2];table(node_inet$species)[1];

###################################
# Roth SARS-CoV-2 interactome
node_roth <- read.csv("~/Documents/INET-work/virus_network/Y2H_screening/source_data/dk_node.csv", header = T, as.is = T)
edge_roth <- read.csv("~/Documents/INET-work/virus_network/Y2H_screening/source_data/dk_BFG.csv", header = T, as.is = T)
net_roth <- graph_from_data_frame(d = edge_roth, vertices = node_roth, directed = F)
V(net_roth)$species <- node_roth$species
V(net_roth)$color <- ifelse(V(net_roth)$species == "human", "blue", "pink")
V(net_roth)[1]; E(net_roth)[1]; table(node_roth$species)[2];table(node_roth$species)[1];

###################################
# Gordon SARS-CoV-2 interactome
node_gordon <- read.csv("~/Documents/INET-work/virus_network/Y2H_screening/source_data/gordon_node.csv", header = T, as.is = T)
edge_gordon <- read.csv("~/Documents/INET-work/virus_network/Y2H_screening/source_data/gordon_APMS.csv", header = T, as.is = T)
net_gordon <- graph_from_data_frame(d = edge_gordon, vertices = node_gordon, directed = F)
V(net_gordon)$species <- node_gordon$species
V(net_gordon)$color <- ifelse(V(net_gordon)$species == "human", "blue", "pink")
V(net_gordon)[1]; E(net_gordon)[1]; table(node_gordon$species)[2];table(node_gordon$species)[1];

###################################
# Li SARS-CoV-2 interactome
node_li <- read.csv("~/Documents/INET-work/virus_network/Y2H_screening/source_data/li_node.csv", header = T, as.is = T)
edge_li <- read.csv("~/Documents/INET-work/virus_network/Y2H_screening/source_data/li_APMS.csv", header = T, as.is = T)
net_li <- graph_from_data_frame(d = edge_li, vertices = node_li, directed = F)
V(net_li)$species <- node_li$species
V(net_li)$color <- ifelse(V(net_li)$species == "human", "blue", "pink")
V(net_li)[1]; E(net_li)[1]; table(node_li$species)[2];table(node_li$species)[1];

###################################
# Stukalov SARS-CoV-2 interactome
node_stukalov <- read.csv("~/Documents/INET-work/virus_network/Y2H_screening/source_data/stukalov_node.csv", header = T, as.is = T)
edge_stukalov <- read.csv("~/Documents/INET-work/virus_network/Y2H_screening/source_data/stukalov_APMS.csv", header = T, as.is = T)
# len <- dim(edge_stukalov)[1]
# node_stukalov <- unique(data.frame(symbol = c(edge_stukalov[, 2], edge_stukalov[, 1]), species = c(rep("human", len), rep("virus", len)))) # generate node file, 'cause excel doing shit!
net_stukalov <- graph_from_data_frame(d = edge_stukalov, vertices = node_stukalov, directed = F)
V(net_stukalov)$species <- node_stukalov$species
V(net_stukalov)$color <- ifelse(V(net_stukalov)$species == "human", "blue", "pink")
V(net_stukalov)[1]; E(net_stukalov)[1]; table(node_stukalov$species)[2];table(node_stukalov$species)[1];

###################################
# Laurent SARS-CoV-2 interactome
node_laurent <- read.csv("~/Documents/INET-work/virus_network/Y2H_screening/source_data/laurent_node.csv", header = T, as.is = T)
edge_laurent <- read.csv("~/Documents/INET-work/virus_network/Y2H_screening/source_data/laurent_BioID.csv", header = T, as.is = T)
net_laurent <- graph_from_data_frame(d = edge_laurent, vertices = node_laurent, directed = F)
V(net_laurent)$species <- node_laurent$species
V(net_laurent)$color <- ifelse(V(net_laurent)$species == "human", "blue", "pink")
V(net_laurent)[1]; E(net_laurent)[1]; table(node_laurent$species)[2];table(node_laurent$species)[1];

###################################
# merge edges into one dataframe
df_inet <- as_data_frame(net_inet)
df_inet <- unite(df_inet, "INET", c(from, to), remove = FALSE, sep = "-")
df_roth <- as_data_frame(net_roth)
df_roth <- unite(df_roth, "Roth", c(from, to), remove = FALSE, sep = "-")
df_gordon <- as_data_frame(net_gordon)
df_gordon <- unite(df_gordon, "Gordon", c(from, to), remove = FALSE, sep = "-")
df_li <- as_data_frame(net_li)
df_li <- unite(df_li, "Li", c(from, to), remove = FALSE, sep = "-")
df_stukalov <- as_data_frame(net_stukalov)
df_stukalov <- unite(df_stukalov, "Stukalov", c(from, to), remove = FALSE, sep = "-")
df_laurent <- as_data_frame(net_laurent)
df_laurent <- unite(df_laurent, "Laurent", c(from, to), remove = FALSE, sep = "-")

## manually merge these edges into one single file, the default data frame in R doesn't accept data frame with different row length
write.table(df_inet$INET, file = "/tmp/1.csv", sep = ",", row.names = F)
write.table(df_roth$Roth, file = "/tmp/2.csv", sep = ",", row.names = F)
write.table(df_gordon$Gordon, file = "/tmp/3.csv", sep = ",", row.names = F)
write.table(df_li$Li, file = "/tmp/4.csv", sep = ",", row.names = F)
write.table(df_stukalov$Stukalov, file = "/tmp/5.csv", sep = ",", row.names = F)
write.table(df_laurent$Laurent, file = "/tmp/6.csv", sep = ",", row.names = F)
