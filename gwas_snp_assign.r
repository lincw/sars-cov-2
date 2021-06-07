# GWAS hit candidates plus 1st interactors from HuRI
# Lin Chung-wen

######
# load package
library(openxlsx)
library(igraph)
library(rethinking)
library(linkcomm)
source("~/Documents/INET-work/virus_network/src/plotOCGGraph.r")

######
# functions
# 1. degree preserving randomized HuRI network
# 2. count final node list appeared in HuSCI
# 3. extract the largest connected component
# 4. community detection with OCG,
random <- function(node) {
    df <- c()
    huri_deg_g <- rewire(huri_g, keeping_degseq(niter = gsize(huri_g) * 10))
    gwas_random_final <- combineNetwork(huri_deg_g, node)
    df[1] <- vcount(gwas_random_final)
    df[2] <- ecount(gwas_random_final)
    strong <- components(gwas_random_final, mode = "strong")
    strong_graph <- induced_subgraph(gwas_random_final, names(strong$membership[strong[1]$membership == 1]))
    df[3] <- vcount(strong_graph)
    df[4] <- ecount(strong_graph)
    return(df)
}
source("~/Documents/INET-work/virus_network/src/combineNetwork.r")

######
# load dataset
huri <- read.csv("~/Documents/INET-work/references/HuRI_binaryPPI/HuRI_Tong_withSymbol.csv", header = T)
gwas <- read.xlsx("~/Documents/INET-work/virus_network/statistic_results/GWAS/COVID_GWAS hits_v2.xlsx")
husci <- read.csv("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/binary_node_1126.csv", header = TRUE)
gwas_2 <- read.xlsx("~/Documents/INET-work/virus_network/statistic_results/GWAS/COVID_GWAS_hit_inHUSCI_v2.xlsx") # to filter hits with high degree
    # the gwas_2 was generated after 1st analysis
gwas_file <- "~/Documents/INET-work/virus_network/statistic_results/GWAS/COVID_GWAS hits_source_v2.xlsx"
candidate <- read.xlsx(gwas_file, sheet = "all")
ctcl <- read.xlsx(gwas_file, sheet = "critical")
hosp <- read.xlsx(gwas_file, sheet = "hospitalization")
infct <- read.xlsx(gwas_file, sheet = "infection")

######
# HuRI graph generation
huri_symbol <- huri[, c(5:6)]
huri_g <- graph_from_data_frame(huri_symbol, directed = FALSE) # V:8274, E:52573
huri_g <- simplify(huri_g, remove.loops = FALSE) # V:8274, E:52558

######
# protein list filter
husci_sym <- husci[husci$group == "human", "node"]
husci_huri <- V(huri_g)$name[V(huri_g)$name %in% husci_sym] # HuSCI in HuRI whole
gwas_huri <- gwas$All.LD[gwas$All.LD %in% V(huri_g)$name] # GWAS hit in HuRI
gwas_all_candidate_huri <- gwas_huri[gwas_huri %in% V(huri_g)$name]
######
# interactor of GWAS hit
gwas_hit_1st <- make_ego_graph(huri_g, nodes = sort(gwas$All.LD[gwas$All.LD %in% V(huri_g)$name]), order = 1, mode = "all") #17 of 42 GWAS hits in HuRI
# subnetwork of GWAS hit from HuRI
all_graph_merge <- combineNetwork(huri_g, gwas_all_candidate_huri)

######
# community detection with OCG
all_strong <- components(all_graph_merge, mode = "strong")
all_largest <- induced_subgraph(all_graph_merge, names(all_strong$membership[all_strong[1]$membership == 1]))
all_ocg <- getOCG.clusters(as_data_frame(all_largest), init.class.sys = 1, cent.class.sys = 0)
all_ocg2 <- getOCG.clusters(as_data_frame(all_largest), cent.class.sys = 0)

plotOCGraph(all_ocg, main = "V2 with maximal cliques")
plotOCGraph(all_ocg2, main = "V2 with Centered cliques")

######
# randomize HuRI and extract the subnetwork from GWAS hits to their 1st interactors
random_gwas_all <- c()
random_gwas_all <- c(random_gwas_all, mcreplicate(10000, random(gwas_all_candidate_huri), mc.cores = detectCores()))
random_gwas_all <- as.data.frame(matrix(random_gwas_all, ncol = 4, byrow = T))
names(random_gwas_all) <- c("ncount_subnetwork", "ecount_subnetwork", "ncount_largest", "ecount_largest")

pdf("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/random_GWAS_hit_v2.pdf", width = 3, height = 3)
dens_gwas <- hist(random_gwas_all[, 2], breaks = 15, plot = FALSE)
par(mgp = c(2, 0.7, 0), ps = 8)
plot(dens_gwas, col = rgb(0.75, 0.75, 0.75, 1/2), border = NA, las = 1, xlab = "Number of edges", ylab = "Frequency density", main = "", yaxt = "n", cex.sub = 0.5)
mytitle <- "COVID19 GWAS hit"
mysubtitle <- "version 2"
mtext(side = 3, line = 1, cex = 1, mytitle)
mtext(side = 3, line = 0.5, cex = .7, mysubtitle)
axis(side = 2, at = seq(0, 2200, 400), labels = seq(0, 0.2, 0.04), las = 2)
arrows(370, 300, 370, 10, col = "#922687", lwd = 2, length = 0.1)
text(370 - 30, 400, paste0("p = ", 4/10000), cex = 0.4, pos = 4)
dev.off()

######
# get one gene for individual SNP
## 1. all GWAS hits
all_candidate <- expand.grid(strsplit(candidate[c(5:13), 5], split = ","))
all_candidate_uniq <- unname(unique(as.matrix(unlist(all_candidate)))[, 1])
all_candidate_huri <- all_candidate_uniq[all_candidate_uniq %in% V(huri_g)$name]
## 2. Critical illness
ctcl_candidate <- expand.grid(strsplit(ctcl[c(3:6), 5], split = ","))
ctcl_candidate_uniq <- unname(unique(as.matrix(unlist(ctcl_candidate)))[, 1])
ctcl_candidate_huri <- ctcl_candidate_uniq[ctcl_candidate_uniq %in% V(huri_g)$name]
## 3. Hospitalization
hosp_candidate <- expand.grid(strsplit(hosp[c(5:9), 5], split = ","))
hosp_candidate_uniq <- unname(unique(as.matrix(unlist(hosp_candidate)))[, 1])
hosp_candidate_huri <- hosp_candidate_uniq[hosp_candidate_uniq %in% V(huri_g)$name]
## 4. Reported infection
infct_candidate <- expand.grid(strsplit(infct[c(3:7), 5], split = ","))
infct_candidate_uniq <- unname(unique(as.matrix(unlist(infct_candidate)))[, 1])
infct_candidate_huri <- infct_candidate_uniq[infct_candidate_uniq %in% V(huri_g)$name]
## check subnetwork from GWAS candidates, one for each SNP
to_check <- list(All_GWAS = all_candidate, CTCL_GWAS = ctcl_candidate, HOSP_GWAS = hosp_candidate, INFCT_GWAS = infct_candidate)
for (i in 1:length(to_check)) {
    pdf(paste0("/tmp/", names(to_check[i]), ".pdf"))
    count_check <- list()
    for (j in 1:dim(to_check[[i]])[1]) {
        genes <- as.vector(as.matrix(unname(to_check[[i]][j, ])))
        if(table(genes %in% V(huri_g)$name)["TRUE"] == dim(to_check[[i]])[2]) {
            check <- c()
            network <- combineNetwork(huri_g, genes)
            plot(network, vertex.size = 4, vertex.label.cex = .5, vertex.label.dist = 1)
            title(paste0(genes, collapse = "-"), cex.main = 0.8)
            mtext(paste0("vcount: ", vcount(network)), side = 3, line = 0.5, cex = 0.5)
            mtext(paste0("ecount: ", ecount(network)), side = 3, line = -0.5, cex = 0.5)
            count_check[[paste0(genes, collapse = "-")]] <- c(vcount(network), ecount(network))
        }
    }
    dev.off()
    check_df <- data.frame(matrix(unlist(count_check), ncol = 2, byrow = T))
    row.names(check_df) <- names(count_check)
    colnames(check_df) <- c("vcount", "ecount")
    assign(names(to_check[i]), check_df)
}
par(mfrow = c(1, 4))
boxplot(All_GWAS, las = 2)
boxplot(CTCL_GWAS, las = 2)
boxplot(HOSP_GWAS, las = 2)
boxplot(INFCT_GWAS, las = 2)
