library(openxlsx)
library(linkcomm)
library(pheatmap)
# library(ComplexHeatmap)
library(circlize) # color used for ComplexHeatmap
library(seriation)
library(dendextend)
# library(ctc) # output hclust data as .cdt for java treeview
library(foreach)
library(doParallel)

jaccard <- function(a, b) {
    intersection <- length(intersect(a, b))
    union <- length(a) + length(b) - intersection
    return (intersection/union)
}

cores <- detectCores()
cl <- makeCluster(cores)
registerDoParallel(cl)

######
# load OCG data
huri_ocg <- readRDS("~/Documents/INET-work/HuRI/HuRI_ocg.RDS")

# id mapping table
id_table <- read.xlsx("~/Documents/INET-work/virus_network/statistic_results/community/HuRI_GO_annotation.xlsx", sheet = "id")

tar_en <- read.xlsx("~/Documents/INET-work/virus_network/statistic_results/community/HuRI_communities_withHuSCI.xlsx", sheet = 2)
husci_comm <- tar_en[tar_en$p < 0.05 & tar_en$size >= 4, "cluster"] # n = 204
non_comm <- tar_en[is.na(tar_en$viral_target) & tar_en$size >= 4, "cluster"] # n = 2576

comm_member_husci <- list()
for (i in 1:length(husci_comm)) {
    comm_member_husci[[i]] <- getNodesIn(huri_ocg, clusterids = husci_comm[i])
}
names(comm_member_husci) <- husci_comm

jac_husci <- foreach(i = 1:length(husci_comm)) %dopar% {
    sapply(comm_member_husci, function(x) {jaccard(comm_member_husci[[i]], x)})
}
jac_husci_mt <- as.matrix(do.call(rbind, jac_husci))
rownames(jac_husci_mt) <- husci_comm
colnames(jac_husci_mt) <- husci_comm

jac_husci_2 <- c()
for (i in 1:length(husci_comm)) {
    for (j in 1:length(husci_comm)) {
        jac_husci_2 <- c(jac_husci_2, jaccard(comm_member_husci[[i]], comm_member_husci[[j]]))
    }
}
jac_husci2_mt <- matrix(jac_husci_2, nrow = length(husci_comm), byrow = T)
rownames(jac_husci2_mt) <- husci_comm
colnames(jac_husci2_mt) <- husci_comm

comm_member_noTar <- list()
for (i in 1:length(non_comm)) {
    comm_member_noTar[[i]] <- getNodesIn(huri_ocg, clusterids = non_comm[i])
}

jac_noTar <- foreach(i = 1:length(non_comm)) %dopar% {
    sapply(comm_member_noTar, function(x) {jaccard(comm_member_noTar, x)})
    # setTxtProgressBar(txtProgressBar(min = 0, max = length(non_comm), style = 3), i)
}
stopCluster(cl)

save.image(file = "overlap.RData")