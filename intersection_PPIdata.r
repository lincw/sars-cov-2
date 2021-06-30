# interaction and interactor overlap between SARS-CoV-2 interactomes
# Lin Chung-wen

######
# load package
library(openxlsx)
library(dplyr)
library(igraph)

######
# load date
husci <- read.csv("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/binary_edge_1126.csv", header = T)
husci <- mutate_all(husci, .funs = toupper)
others <- "/Users/chung-wen.lin/Documents/INET-work/virus_network/toSummarize/PPIs/Extended_Table_2_PPIs.xlsx"
gordon <- read.xlsx(others, sheet = "Gordon")
gordon <- mutate_all(gordon, .funs = toupper)

stukalov <- read.xlsx(others, sheet = "Stukalov")
stukalov <- mutate_all(stukalov, .funs = toupper)

li <- read.xlsx(others, sheet = "Li")
li <- mutate_all(li, .funs = toupper)

nabeel <- read.xlsx(others, sheet = "Nabeel")
nabeel <- mutate_all(nabeel, .funs = toupper)

laurent <- read.xlsx(others, sheet = "Laurent")
laurent <- mutate_all(laurent, .funs = toupper)

stgermain <- read.xlsx(others, sheet = "St_Germain")
stgermain <- mutate_all(stgermain, .funs = toupper)

samavarchi <- read.xlsx(others, sheet = "Samavarchi")
samavarchi <- mutate_all(samavarchi, .funs = toupper)

######
# interaction overlap
husci_g <- graph_from_data_frame(husci[, c(1:2)], direct = F)
gordon_g <- graph_from_data_frame(gordon[, c(2, 4)], direct = F)
stukalov_g <- graph_from_data_frame(stukalov[, c(1:2)], direct = F)
li_g <- graph_from_data_frame(li[, c(1:2)], direct = F)
nabeel_g <- graph_from_data_frame(nabeel[, c(1, 3)], direct = F)
laurent_g <- graph_from_data_frame(laurent[, c(1:2)], direct = F)
stgermain_g <- graph_from_data_frame(stgermain[, c(1:2)], direct = F)
samavarchi_g <- graph_from_data_frame(samavarchi[, c(1:2)], direct = F)

husci_gordon_edge <- intersection(husci_g, gordon_g)
husci_stukalov_edge <- intersection(husci_g, stukalov_g)
husci_li_edge <- intersection(husci_g, li_g)
husci_nabeel_edge <- intersection(husci_g, nabeel_g)
husci_laurent_edge <- intersection(husci_g, laurent_g)
husci_stgermain_edge <- intersection(husci_g, stgermain_g)
husci_samavarchi_edge <- intersection(husci_g, samavarchi_g)

######
# interactor overlap
husci_gordon_node <- unique(as_edgelist(husci_g)[, 2])[unique(as_edgelist(husci_g)[, 2]) %in% unique(as_edgelist(gordon_g)[, 2])]
husci_stukalov_node <- unique(as_edgelist(husci_g)[, 2])[unique(as_edgelist(husci_g)[, 2]) %in% unique(as_edgelist(stukalov_g)[, 2])]
husci_li_node <- unique(as_edgelist(husci_g)[, 2])[unique(as_edgelist(husci_g)[, 2]) %in% unique(as_edgelist(li_g)[, 2])]
husci_nabeel_node <- unique(as_edgelist(husci_g)[, 2])[unique(as_edgelist(husci_g)[, 2]) %in% unique(as_edgelist(nabeel_g)[, 2])]
husci_laurent_node <- unique(as_edgelist(husci_g)[, 2])[unique(as_edgelist(husci_g)[, 2]) %in% unique(as_edgelist(laurent_g)[, 2])]
husci_stgermain_node <- unique(as_edgelist(husci_g)[, 2])[unique(as_edgelist(husci_g)[, 2]) %in% unique(as_edgelist(stgermain_g)[, 2])]
husci_samavarchi_node <- unique(as_edgelist(husci_g)[, 2])[unique(as_edgelist(husci_g)[, 2]) %in% unique(as_edgelist(samavarchi_g)[, 2])]

all4 <- unique(c(husci_gordon_node, husci_stukalov_node, husci_li_node, husci_nabeel_node))
all7 <- unique(c(husci_gordon_node, husci_stukalov_node, husci_li_node, husci_nabeel_node, husci_laurent_node, husci_stgermain_node, husci_samavarchi_node))
