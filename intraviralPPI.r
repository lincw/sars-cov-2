# Intraviral
# Lin Chung-wen
# 2021.07.08

library(igraph)
sars2 <- read.xlsx("../../cloning/data/SARS2_seq_supplementary.xlsx")
sars2$"viral.protein"[c(11, 12)] <- "NSP12"
sars2 <- toupper(unique(sars2$"viral.protein"))

y2h_intra <- read.csv("../../cloning/data/intraviral_edge.csv", header = T)
y2h_intra <- data.frame(vAD = toupper(y2h_intra$vAD), vDB = toupper(y2h_intra$vDB))
y2h_graph <- graph_from_data_frame(y2h_intra)
y2h_intra <- unique(c(toupper(y2h_intra$vAD), toupper(y2h_intra$vDB)))

intra_viral <- read.csv("../../cloning/data/Y2H_list.csv", header = T)
intra_graph <- graph_from_data_frame(intra_viral)
intra_viral <- unique(c(intra_viral$FROM, intra_viral$TO))

binary_intra <- gsize(intersection(intra_graph, y2h_graph))

sars2_rep <- function(x) {
  sars2_sample <- data.frame(from = sample(x, 27, replace = TRUE), to = sample(x, 27, replace = TRUE))
  sars2_sample_graph <- simplify(graph_from_data_frame(sars2_sample))
  return(gsize(intersection(intra_graph, sars2_sample_graph)))
}

random_intra <- mcreplicate(10000, sars2_rep(sars2), mc.cores = detectCores())

paste("Observed in intra_viral identified human target within VirHostome:", binary_intra)
sign_intra <- round(1 - (abs(sig(as.numeric(random_intra), as.numeric(binary_intra))) / 10000), 4)
