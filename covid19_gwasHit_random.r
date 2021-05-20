# COVID-19 GWAS hit analysis
# Lin Chung-win

## load package ----
library(openxlsx)
library(igraph)
library(rethinking)

## load dataset ----
huri <- read.csv("~/Documents/INET-work/references/HuRI_binaryPPI/HuRI_Tong_withSymbol.csv", header = T)
gwas <- read.xlsx("~/Documents/INET-work/virus_network/statistic_results/GWAS/COVID_GWAS hits.xlsx")
husci <- read.csv("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/binary_node_1126.csv", header = TRUE)
gwas_2 <- read.xlsx("~/Documents/INET-work/virus_network/statistic_results/GWAS/COVID_GWAS_hit_inHUSCI.xlsx") # to filter hits with high degree
    # the gwas_2 was generated after 1st analysis

## HuRI graph generation ----
huri_symbol <- huri[, c(5:6)]
huri_g <- graph_from_data_frame(huri_symbol, directed = FALSE) # V:8274, E:52573
huri_g <- simplify(huri_g, remove.loops = FALSE) # V:8274, E:52558

## protein list filter ----
husci_sym <- husci[husci$group == "human", "node"]
husci_huri <- V(huri_g)$name[V(huri_g)$name %in% husci_sym] # HuSCI in HuRI whole
gwas_huri <- gwas$All.LD[gwas$All.LD %in% V(huri_g)$name] # GWAS hit in HuRI

## function generation----
check <- function(x) {
    value <- table(V(x)$name %in% husci_sym)[2]
    value <- ifelse(is.na(value), "0", value)
    return(as.numeric(value))
}

check_update <- function(x) {
    V(x)$name[V(x)$name %in% husci_sym]
}

## interactor of GWAS hit
gwas_hit_1st <- make_ego_graph(huri_g, nodes = sort(gwas$All.LD[gwas$All.LD %in% V(huri_g)$name]), order = 1, mode = "all") #16 of 48 in HuRI
gwas_check <- length(unique(unlist(lapply
    (gwas_hit_1st, check_update)))) #11 interactors

## filter interactor of GWAS hit
gwas_hit_1st_sel <- make_ego_graph(huri_g, nodes = sort(gwas_2[gwas_2[, 2] < 30 & gwas_2[, 2] > 0, 1]), order = 1, mode = "all") #100 for exclude MUC1; 30 for exclude MUC1&TMEM65
gwas_check_sel <- length(unique(unlist(lapply
    (gwas_hit_1st_sel, check_update)))) #11 from exclude MUC1; 10 from exclude MUC1&TMEM65

## protein list filter 2 ----
gwas_husci <- unique(unlist(lapply
        (gwas_hit_1st_sel, check_update))
        ) # HuSCI in GWAS hit interactor,
gwas_1st_hit_l <- unique(unlist(lapply(gwas_hit_1st_sel, function(x)
    V(x)$name))) # GWAS 1st hit list,

gwas_n_husci <- unique(gwas_1st_hit_l[
        !gwas_1st_hit_l %in% gwas_husci]) # GWAS hit interactor not in HuSCI
huri_nhusci <- unique(V(huri_g)$name[!V(huri_g)$name %in% husci_sym])
huri_nhusci_n1st <- unique(huri_nhusci[!huri_nhusci %in% gwas_1st_hit_l])
                #1. HuRI not in HuSCI: 8274 - 146
                #2. HuRI not in gwas hit 1st list (GWAS hit involved in 1st list): 8274 - 146 - 209 + 11(ruplicate)

## degree preserving randomized network analysis ----
random_check <- function(sel = 0) {
    huri_deg_g <- rewire(huri_g, keeping_degseq(niter = gsize(huri_g) * 10))
    if (sel == 0) {
        gwas_hit_deg_1st <- make_ego_graph(huri_deg_g, gwas$All.LD[gwas$All.LD %in% V(huri_g)$name], order = 1, mode = "all")
    } else if (sel == 100) {
        gwas_hit_deg_1st <- make_ego_graph(huri_deg_g, gwas_2[gwas_2[, 2] < 100 & gwas_2[, 2] > 0, 1], order = 1, mode = "all")
    } else if (sel == 30) {
        gwas_hit_deg_1st <- make_ego_graph(huri_deg_g, gwas_2[gwas_2[, 2] < 30 & gwas_2[, 2] > 0, 1], order = 1, mode = "all")
    }
    gwas_deg_check <- length(unique(unlist(lapply(gwas_hit_deg_1st, check_update))))
    return(gwas_deg_check)
}

random_out <- c()
random_out <- c(random_out, mcreplicate(10000, random_check(0), mc.cores = detectCores()))

random_out_sel <- c()
random_out_sel <- c(random_out_sel, mcreplicate(10000, random_check(30), mc.cores = detectCores()))

## fisher test ----

df <- matrix(c(
    length(gwas_husci), # gwas hit interactor in HuSCI
    length(gwas_n_husci), # gwas hit interactor not in HuSCI
    length(husci_huri[!husci_huri %in% gwas_husci]), # HuSCI in HuRI, not in gwas hit interactor
    length(huri_nhusci_n1st)) # not_gwas_not_husci
    , ncol = 2, byrow = T)

fisher.test(df)
## all
# p-value: 0.001139
# odds ratio: 3.26253

## exclude MUC1
# p-value: 5.236e-06
# odds ratio: 6.283085

## exclude MUC1 and TMEM65
# p-value: 1.739e-06
# odds ratio: 8.107661

# plot permutation analysis ----
pdf("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/random_GWAS_hit_2.pdf", width = 3, height = 3)
dens_gwas <- hist(random_out_sel, breaks = 15, plot = FALSE)
par(mgp = c(2, 0.7, 0), ps = 8)
plot(dens_gwas, xlim = c(0, 20), ylim = c(0, 0.22), col = rgb(0.75, 0.75, 0.75, 1/2), freq = FALSE, border = NA, las = 1, xlab = "Number of interactor in HuSCI", ylab = "Frequency density", main = "COVID19 GWAS hit\nMUC1&TMEM65 excluded", xaxt = "n")
axis(side = 1, at = seq(1, 20, 4) - 0.5, labels = seq(0, 19, 4))
box(col = "black")
arrows(length(gwas_husci) + 0.5, 400/10000, length(gwas_husci) + 0.5, 50/10000, col = "red", lwd = 2, length = 0.1)
text(length(gwas_husci) - 2, 550/10000, paste0("p = ", 0.0042), cex = 0.4, pos = 4)
dev.off()