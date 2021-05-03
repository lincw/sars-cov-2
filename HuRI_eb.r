library(igraph)
library(openxlsx)

# original
# huri <- read.csv("HuRI_Tong_withSymbol.csv", header = T)
# huri <- huri[, c(5:6)]
# huri_graph <- graph_from_data_frame(huri, directed = FALSE)
# huri_eb <- cluster_edge_betweenness(huri_graph, directed = FALSE)
# save.image(file = "HuRI_edgebetweenness.RData")

# **2020.11.23 22:26** use individual assay
huri_assay1 <- read.xlsx("~/Documents/INET-work/references/HuRI_binaryPPI/HuRI_individual_assay.xlsx", sheet = 1)
huri_assay2 <- read.xlsx("~/Documents/INET-work/references/HuRI_binaryPPI/HuRI_individual_assay.xlsx", sheet = 2)
huri_assay3 <- read.xlsx("~/Documents/INET-work/references/HuRI_binaryPPI/HuRI_individual_assay.xlsx", sheet = 3)
huri <- list()
huri[[1]] <- graph_from_data_frame(huri_assay1[, c(3, 4)], directed = FALSE)
huri[[2]] <- graph_from_data_frame(huri_assay2[, c(3, 4)], directed = FALSE)
huri[[3]] <- graph_from_data_frame(huri_assay3[, c(3, 4)], directed = FALSE)

huri_eb1 <- cluster_edge_betweenness(huri[[1]], directed = FALSE)
huri_eb2 <- cluster_edge_betweenness(huri[[2]], directed = FALSE)
huri_eb3 <- cluster_edge_betweenness(huri[[3]], directed = FALSE)
save.image(file = "~/Documents/INET-work/virus_network/statistic_results/HuRI_edgebetweenness_assays.RData")

# output membership into excel worksheet
library(openxlsx)
ceb_tb1 <- membership(huri_eb1)
ceb_df1 <- as.data.frame(matrix(ceb_tb1))
names(ceb_df1) <- "community"
ceb_df1$gene <- names(ceb_tb1)
ceb_df1 <- ceb_df1[order(ceb_df1$community), ]

ceb_tb2 <- membership(huri_eb2)
ceb_df2 <- as.data.frame(matrix(ceb_tb2))
names(ceb_df2) <- "community"
ceb_df2$gene <- names(ceb_tb2)
ceb_df2 <- ceb_df2[order(ceb_df2$community), ]

ceb_tb3 <- membership(huri_eb3)
ceb_df3 <- as.data.frame(matrix(ceb_tb3))
names(ceb_df3) <- "community"
ceb_df3$gene <- names(ceb_tb3)
ceb_df3 <- ceb_df3[order(ceb_df3$community), ]

wb <- createWorkbook()
addWorksheet(wb, "assay1")
addWorksheet(wb, "assay2")
addWorksheet(wb, "assay3")
writeData(wb, sheet = "assay1", ceb_df1)
writeData(wb, sheet = "assay2", ceb_df2)
writeData(wb, sheet = "assay3", ceb_df3)
saveWorkbook(wb, "/tmp/HuRI_assays_edgebetweenness.xlsx", overwrite = TRUE)