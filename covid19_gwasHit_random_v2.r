# COVID-19 GWAS hit analysis
# Lin Chung-win
## v2: this is very wired! One typo from Pascal makes me unable to replicate the analysis?!

######
# load package
library(openxlsx)
library(igraph)
library(rethinking)
library(gprofiler2)
library(linkcomm)
library(gplots)
# functions
source("~/Documents/INET-work/virus_network/src/plotOCGGraph.r")

go2check <- function(x, organism) {
    goquery <- gost(query = x, organism = organism, correction_method = "bonferroni", evcodes = T)
    goquery$result$inCommunity <- paste0( goquery$meta$query_metadata$queries$query_1, collapse = ",")
    goquery$result$annotatedInCommunity <- paste0(goquery$meta$genes_metadata$query$query_1$ensgs, collapse = ",")
    return(goquery$result)
}

check <- function(x) {
    value <- table(V(x)$name %in% husci_sym)[2]
    value <- ifelse(is.na(value), "0", value)
    return(as.numeric(value))
}

checkUpdate <- function(x) {
    V(x)$name[V(x)$name %in% husci_sym]
}

source("~/Documents/INET-work/virus_network/src/combineNetwork.r")

huriRewire <- function(remove.loops = FALSE, ...) {
    huri_re <- rewire(huri_g, keeping_degseq(niter = gsize(huri_g) * 10))
    huri_sim <- simplify(huri_re, remove.loops = remove.loops)
    return(huri_sim)
}

huriRewireHusci <- function(node, remove.loops) {
    huri_re <- huriRewire(remove.loops)
    merged <- combineNetwork(huri_re, node)
    merged_inHuSCI <- table(V(merged)$name %in% husci_sym)["TRUE"]
    return(merged_inHuSCI)
}

huriRewireHusciMulti <- function(node, node1, node2, node3, remove.loops) {
    df <- c()
    huri_re <- huriRewire(remove.loops)
    merged <- combineNetwork(huri_re, node)
    df <- c(df, table(V(merged)$name %in% husci_sym)["TRUE"])
    merged <- combineNetwork(huri_re, node1)
    df <- c(df, table(V(merged)$name %in% husci_sym)["TRUE"])
    merged <- combineNetwork(huri_re, node2)
    df <- c(df, table(V(merged)$name %in% husci_sym)["TRUE"])
    merged <- combineNetwork(huri_re, node3)
    df <- c(df, table(V(merged)$name %in% husci_sym)["TRUE"])
    return(df)
    # 1. viral targets from whole GWAS hits
    # 2. viral targets from critical illness
    # 3. viral targets from hospitalization
    # 4. viral targets from reported infection
}

ocgRewire <- function(node) {
    rewire <- huriRewire()
    combine <- combineNetwork(rewire, node)
    strong <- components(combine, mode = "strong")
    strong1 <- induced_subgraph(combine, names(strong$membership[strong$membership == 1]))
    strong_ocg <- getOCG.clusters(as_data_frame(strong1), init.class.sys = 3, cent.class.sys = 0)
    return(strong_ocg)
}
ocgRewireRatio <- function(node) {
    df <- c()
    strong_ocg <- ocgRewire(node)
    data <- as.numeric(sort(strong_ocg$clustsizes, decreasing = TRUE))
    df <- c(df, data)
    if(data[2] > 4) {
        df <- c((data[1] / data[2]), df)
    } else {
        df <- c("NA", df)
    }
    return(df)
}
# edit plot parameters
trace("plot.igraph", edit = T)
######
# load dataset
huri <- read.csv("~/Documents/INET-work/references/HuRI_binaryPPI/HuRI_Tong_withSymbol.csv", header = T)
gwas <- read.xlsx("~/Documents/INET-work/virus_network/statistic_results/GWAS/COVID_GWAS hits_v2.xlsx")
husci <- read.csv("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/binary_node_1126.csv", header = TRUE)
gwas_2 <- read.xlsx("~/Documents/INET-work/virus_network/statistic_results/GWAS/COVID_GWAS_hit_inHUSCI_v2.xlsx") # to filter hits with high degree
    # the gwas_2 was generated after 1st analysis
gwas_file <- "~/Documents/INET-work/virus_network/statistic_results/GWAS/COVID_GWAS hits_source_v2.xlsx"
gwas_candidate <- read.xlsx(gwas_file, sheet = "all")
gwas_ortholog <- read.xlsx(gwas_file, sheet = "Ortholog")
gwas_ctcl <- read.xlsx(gwas_file, sheet = "critical")
gwas_hosp <- read.xlsx(gwas_file, sheet = "hospitalization")
gwas_infct <- read.xlsx(gwas_file, sheet = "infection")
######
# 1. HuRI graph generation
huri_symbol <- huri[, c(5:6)]
huri_g_ori <- graph_from_data_frame(huri_symbol, directed = FALSE) # V:8274, E:52573
huri_g <- simplify(huri_g_ori, remove.loops = FALSE) # V:8274, E:52558
huri_g_noloop <- simplify(huri_g_ori) # V:8274 E:52078 without self-loops
# protein list filter
husci_sym <- husci[husci$group == "human", "node"]
husci_huri <- V(huri_g)$name[V(huri_g)$name %in% husci_sym] # HuSCI in HuRI whole
gwas_huri <- gwas$All.LD[gwas$All.LD %in% V(huri_g)$name] # GWAS hit in HuRI

######
# 2. interactor of GWAS hit
gwas_hit_1st <- make_ego_graph(huri_g, nodes = sort(gwas$All.LD[gwas$All.LD %in% V(huri_g)$name]), order = 1, mode = "all") #17 of 42 in HuRI
gwas_hit_1st_noLoop <- make_ego_graph(huri_g_noloop, nodes = sort(gwas$All.LD[gwas$All.LD %in% V(huri_g)$name]), order = 1, mode = "all") #17 of 42 in HuRI

# 2.1 interactors of GWAS hits with critical illness phenotypes
ctcl <- gwas[, 2][gwas[, 5] == 1]
ctcl <- unique(ctcl[!is.na(ctcl)])
ctcl_huri <- ctcl[ctcl %in% V(huri_g)$name]
ctcl_1st <- combineNetwork(huri_g, ctcl_huri)
ctcl_1st_noLoop <- simplify(ctcl_1st)
gwas_ctcl_husci <- V(ctcl_1st)$name[V(ctcl_1st)$name %in% husci_sym]
gwas_ctcl_husci_length <- length(gwas_ctcl_husci)

hosp <- gwas[, 2][gwas[, 6] == 1]
hosp <- unique(hosp[!is.na(hosp)])
hosp_huri <- hosp[hosp %in% V(huri_g)$name]
hosp_1st <- combineNetwork(huri_g, hosp_huri)
hosp_1st_noLoop <- simplify(hosp_1st)
gwas_hosp_husci <- V(hosp_1st)$name[V(hosp_1st)$name %in% husci_sym]
gwas_hosp_husci_length <- length(gwas_hosp_husci)

infct <- gwas[, 2][gwas[, 7] == 1]
infct <- unique(infct[!is.na(infct)])
infct_huri <- infct[infct %in% V(huri_g)$name]
infct_1st <- combineNetwork(huri_g, infct_huri)
infct_1st_noLoop <- simplify(infct_1st)
gwas_infct_husci <- V(infct_1st)$name[V(infct_1st)$name %in% husci_sym]
gwas_infct_husci_length <- length(gwas_infct_husci)

overlap_gwas <- list(
    ctcl = ctcl,
    hosp = hosp,
    infct = infct
)
venn(overlap_gwas); title("Overlap between all LD")
overlap_huri <- list(
    ctcl = ctcl_huri,
    hosp = hosp_huri,
    infct = infct_huri
)
venn(overlap_huri); title("Overlap between all LD in HuRI")
overlap_husci <- list(
    ctcl = V(ctcl_1st)$name[V(ctcl_1st)$name %in% husci_sym],
    hosp = V(hosp_1st)$name[V(hosp_1st)$name %in% husci_sym],
    infct = V(infct_1st)$name[V(infct_1st)$name %in% husci_sym]
)
venn(overlap_husci); title("Overlap between 1st node in HuSCI")
######
# 3. **rewiring analysis of HuRI**, to see if the HuSCI viral target is significant.
# load gwas loci info, with 3 phenotype (critical illness, hospitalization and infection)
# subnetwork of GWAS hit from HuRI
# inherite from above code
gwas_all_list_df <- lapply(gwas_hit_1st, as_data_frame)
gwas_all_df <- do.call(rbind, gwas_all_list_df)
gwas_all_g_merge <- graph_from_data_frame(gwas_all_df, directed = FALSE)
# to have interaction between 1st interactors
gwas_all_final <- simplify(induced_subgraph(huri_g, names(V(gwas_all_g_merge))), remove.loops = F) # V:118, E:370; E:351, removing self-loops
gwas_all_final_noLoop <- simplify(gwas_all_final)
gwas_all_husci <- V(gwas_all_final)$name[V(gwas_all_final)$name %in% husci_sym] # 11 viral targets
gwas_all_husci_length <- length(gwas_all_husci)
# visualization of gwas subnetwork
bd <- ifelse(V(gwas_all_final)$name %in% husci_sym, "orange", NA)
label <- ifelse(V(gwas_all_final)$name %in% c(husci_sym, gwas_huri), V(gwas_all_final)$name, NA)
color <- ifelse(V(gwas_all_final)$name %in% gwas_huri, "red", "grey")
pdf("/tmp/GWAS_subnetwork.pdf")
plot(gwas_all_final, vertex.frame.color = bd, vertex.size = 3, vertex.label.dist = 1, vertex.label.color = "black", vertex.label = label, vertex.color = color, vertex.frame.width = 2)
write_graph(gwas_all_final, "/tmp/GWAS_subnetwork.gml", format = "gml")
dev.off()
# 3.1. do statistical analysis for all 17 GWAS hit candidate genes
all_re_husci <- c()
all_re_husci <- c(all_re_husci, mcreplicate(10000, huriRewireHusciMulti(gwas_huri, ctcl_huri, hosp_huri, infct_huri, FALSE), mc.cores = detectCores()))
all_re_husci[is.na(all_re_husci)] <- 0
all_re_husci_df <- as.data.frame(matrix(all_re_husci, ncol = 4, byrow = T))
all_re_husci_noLoop <- c()
all_re_husci_noLoop <- c(all_re_husci_noLoop, mcreplicate(10000, huriRewireHusciMulti(gwas_huri, ctcl_huri, hosp_huri, infct_huri, TRUE), mc.cores = detectCores()))
all_re_husci_noLoop[is.na(all_re_husci_noLoop)] <- 0
all_re_husci_noLoop_df <- as.data.frame(matrix(all_re_husci_noLoop, ncol = 4, byrow = T))

inHuSCI_length <- list(
    allGWAS = gwas_all_husci_length,
    allGWSAS_noLoop = gwas_all_husci_length,
    critical = gwas_ctcl_husci_length,
    critical_noLoop = gwas_ctcl_husci_length,
    hospitalization = gwas_hosp_husci_length,
    hospitalization_noLoop = gwas_hosp_husci_length,
    infection = gwas_infct_husci_length,
    infection_noLoop = gwas_infct_husci_length
)

# 3.2 for critical illness, hospitalization and infection phenotype associated genes
ctcl_re_husci <- c()
ctcl_re_husci <- c(ctcl_re_husci, mcreplicate(10000, huriRewireHusci(ctcl_huri, FALSE), mc.cores = detectCores()))
ctcl_re_husci[is.na(ctcl_re_husci)] <- 0
ctcl_re_husci_noLoop <- c()
ctcl_re_husci_noLoop <- c(ctcl_re_husci_noLoop, mcreplicate(10000, huriRewireHusci(ctcl_huri, TRUE), mc.cores = detectCores()))
ctcl_re_husci_noLoop[is.na(ctcl_re_husci_noLoop)] <- 0

hosp_re_husci <- c()
hosp_re_husci <- c(hosp_re_husci, mcreplicate(10000, huriRewireHusci(hosp_huri, FALSE)), mc.cores = detectCores())
hosp_re_husci[is.na(hosp_re_husci)] <- 0
hosp_re_husci_noLoop <- c()
hosp_re_husci_noLoop <- c(hosp_re_husci_noLoop, mcreplicate(10000, huriRewireHusci(hosp_huri, TRUE), mc.cores = detectCores()))
hosp_re_husci_noLoop[is.na(hosp_re_husci_noLoop)] <- 0

infct_re_husci <- c()
infct_re_husci <- c(infct_re_husci, mcreplicate(10000, huriRewireHusci(infct_huri, FALSE)), mc.cores = detectCores())
infct_re_husci[is.na(infct_re_husci)] <- 0
infct_re_husci_noLoop <- c()
infct_re_husci_noLoop <- c(infct_re_husci_noLoop, mcreplicate(10000, huriRewireHusci(infct_huri, TRUE), mc.cores = detectCores()))
infct_re_husci_noLoop[is.na(infct_re_husci_noLoop)] <- 0

inHuSCI_summary <- list(
    allGWAS = as.numeric(all_re_husci),
    allGWSAS_noLoop = as.numeric(all_re_husci_noLoop),
    ctcl = as.numeric(ctcl_re_husci),
    ctcl_noLoop = as.numeric(ctcl_re_husci_noLoop),
    hosp = as.numeric(hosp_re_husci),
    hosp_noLoop = as.numeric(hosp_re_husci_noLoop),
    infct = as.numeric(infct_re_husci),
    infct_noLoop = as.numeric(infct_re_husci_noLoop)
)
write.xlsx(t(do.call(rbind, inHuSCI_summary)), file = "~/Documents/INET-work/virus_network/statistic_results/GWAS/10000_randomHuRI.xlsx")

# plotting
pdf("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/random_GWAS_viral_target_v2.pdf", width = 3, height = 3)
par(mgp = c(2, 0.7, 0), ps = 8)
dens_gwas <- hist(all_re_husci, breaks = 15, plot = FALSE)
plot(dens_gwas, col = rgb(0.75, 0.75, 0.75, 1/2), freq = FALSE, border = NA, las = 1, xlim = c(0, 20), xlab = "Number of viral targets", ylab = "Frequency density", main = "", cex.sub = 0.5)
mytitle <- "COVID19 GWAS loci candidate genes\n(all candidate genes)"
mtext(side = 3, line = 1, cex = 1, mytitle)
arrows(gwas_all_husci_length + 0.5, 0.04, gwas_all_husci_length + 0.5, 0.02, col = "#922687", lwd = 2, length = 0.1)
text(gwas_all_husci_length - 2, 0.06, paste0("observed = ", gwas_all_husci_length, "\np = ", table(all_re_husci >= gwas_all_husci_length)["TRUE"]/10000), cex = 0.4, pos = 4)

dens_gwas <- hist(ctcl_re_husci, breaks = 15, plot = FALSE)
plot(dens_gwas, col = rgb(0.75, 0.75, 0.75, 1/2), freq = FALSE, border = NA, las = 1, xlim = c(0, 20), xlab = "Number of viral targets", ylab = "Frequency density", main = "", cex.sub = 0.5)
mytitle <- "COVID19 GWAS loci candidate genes\n(critical illness)"
mtext(side = 3, line = 1, cex = 1, mytitle)
arrows(gwas_ctcl_husci_length + 0.5, 0.06, gwas_ctcl_husci_length + 0.5, 0.02, col = "#922687", lwd = 2, length = 0.1)
text(gwas_ctcl_husci_length - 1, 0.08, paste0("observed = ", gwas_ctcl_husci_length, "\np = ", table(ctcl_re_husci >= gwas_ctcl_husci_length)["TRUE"]/10000), cex = 0.4, pos = 4)

dens_gwas <- hist(hosp_re_husci, breaks = 15, plot = FALSE)
plot(dens_gwas, col = rgb(0.75, 0.75, 0.75, 1/2), freq = FALSE, border = NA, las = 1, xlim = c(0, 20), xlab = "Number of viral targets", ylab = "Frequency density", main = "", cex.sub = 0.5)
mytitle <- "COVID19 GWAS loci candidate genes\n(hospitalization)"
mtext(side = 3, line = 1, cex = 1, mytitle)
arrows(gwas_hosp_husci_length + 0.5, 0.04, gwas_hosp_husci_length + 0.5, 0.02, col = "#922687", lwd = 2, length = 0.1)
text(gwas_hosp_husci_length - 2, 0.06, paste0("observed = ", gwas_hosp_husci_length, "\np = ", table(hosp_re_husci >= gwas_hosp_husci_length)["TRUE"]/10000), cex = 0.4, pos = 4)

dens_gwas <- hist(infct_re_husci, breaks = 10, plot = FALSE)
plot(dens_gwas, col = rgb(0.75, 0.75, 0.75, 1/2), freq = FALSE, border = NA, las = 1, xlim = c(0, 10), xlab = "Number of viral targets", ylab = "Frequency density", main = "", cex.sub = 0.5)
mytitle <- "COVID19 GWAS loci candidate genes\n(infection)"
mtext(side = 3, line = 1, cex = 1, mytitle)
arrows(gwas_infct_husci_length + 0.5, 0.1, gwas_infct_husci_length + 0.5, 0.02, col = "#922687", lwd = 2, length = 0.1)
text(gwas_infct_husci_length - 1, 0.2, paste0("observed = ", gwas_infct_husci_length, "\np = ", table(infct_re_husci >= gwas_infct_husci_length)["TRUE"]/10000), cex = 0.4, pos = 4)

dev.off()

boxplot(inHuSCI_summary, las = 1)
pdf("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/random_GWAS_viral_target_v2.1.pdf", width = 4, height = 4)
for (i in 1:length(inHuSCI_summary)) {
    dens_gwas <- hist(inHuSCI_summary[[i]], breaks = seq(0, max(inHuSCI_summary[[i]]), by = 1), plot = FALSE, right = FALSE)
    plot(dens_gwas, col = rgb(0.75, 0.75, 0.75, 1/2), freq = FALSE, border = NA, las = 1, xlim = c(0, max(inHuSCI_summary[[i]])), xlab = "Number of viral targets", ylab = "Frequency density", main = "", cex.sub = 0.5)
    mytitle <- paste0("COVID19 GWAS loci candidate genes\n", names(inHuSCI_summary[i]))
    mtext(side = 3, line = 1, cex = 1, mytitle)
    arrows(inHuSCI_length[[i]] + 0.5, 0.04, inHuSCI_length[[i]] + 0.5, 0.01, col = "#922687", lwd = 2, length = 0.1)
    text(inHuSCI_length[[i]], 0.05, paste0("observed = ", inHuSCI_length[[i]], "\n p = ", table(inHuSCI_summary[[i]] >= inHuSCI_length[[i]])["TRUE"]/10000), cex = 0.4, pos = 4)
}
dev.off()
#######
# 4. keep one ortholog at a time
ortholog_count <- list()
for (i in 1:dim(all_candidate)[1]) {
    ortholog_gra <- combineNetwork(huri_g, c(as.matrix(all_candidate[i, ])))
    V(ortholog_gra)$color <- ifelse(V(ortholog_gra)$name %in% husci_sym, "red", "grey")
    name <- paste0(c(as.matrix(all_candidate[i, ])), collapse = "-")
    ortholog_count[[name]] <- data.frame(Vcount = vcount(ortholog_gra), Ecount = ecount(ortholog_gra))
    pdf(paste0("/tmp/", name, ".pdf"))
    plot(ortholog_gra, vertex.size = 4, vertex.label.cex = .5, vertex.label.dist = 1, vertex.label.color = "black")
    title(name, cex.main = .5)
    dev.off()
}

# subgraph for critical illness, hospitalization and infection
gwas_ctcl_graph <- make_ego_graph(huri_g, nodes = gwas_ctcl_candidate_huri, mode = "all")
gwas_hosp_graph <- make_ego_graph(huri_g, nodes = gwas_hosp_candidate_huri, mode = "all")
gwas_infct_graph <- make_ego_graph(huri_g, nodes = gwas_infct_candidate_huri, mode = "all")
### HuSCI in 1st interactor
gwas_1stN_inHuSCI <- table(names(V(gwas_all_final)) %in% husci_sym)['TRUE']

gwas_ctcl_inHuSCI <- table( unique( unlist( lapply(gwas_ctcl_graph, function(x) names(V(x))) ) ) %in% husci_sym)["TRUE"]
gwas_hosp_inHuSCI <- table( unique( unlist( lapply(gwas_hosp_graph, function(x) names(V(x))) ) ) %in% husci_sym)["TRUE"]
gwas_infct_inHuSCI <- table( unique( unlist( lapply(gwas_infct_graph, function(x) names(V(x))) ) ) %in% husci_sym)["TRUE"]
######
# degree preserving randomized HuRI network and count final node list appeared in HuSCI
randomGwas <- function(network, node) {
    gwas_random_final <- combineNetwork(network, node)
    gwas_1stN_inHuSCI <- table(names(V(gwas_random_final)) %in% husci_sym)['TRUE']
    return(gwas_1stN_inHuSCI)
}
random_gwas_all <- c()
random_gwas_all <- c(random_gwas_all, mcreplicate(10000, randomGwas(huri_g, gwas_all_candidate_huri)), mc.cores = detectCores())
random_gwas_all_final <- as.numeric(random_gwas_all)
random_gwas_all_final[is.na(random_gwas_all_final)] <- 0

random_gwas_ctcl <- c()
random_gwas_ctcl <- c(random_gwas_ctcl, mcreplicate(10000, randomGwas(huri_g, ctcl), mc.cores = detectCores()))
random_gwas_ctcl_final <- as.numeric(random_gwas_ctcl)
random_gwas_ctcl_final[is.na(random_gwas_ctcl_final)] <- 0

random_gwas_hosp <- c()
random_gwas_hosp <- c(random_gwas_hosp, mcreplicate(10000, randomGwas(huri_g, hosp), mc.cores = detectCores()))
random_gwas_hosp_final <- as.numeric(random_gwas_hosp)
random_gwas_hosp_final[is.na(random_gwas_hosp_final)] <- 0

random_gwas_infct <- c()
random_gwas_infct <- c(random_gwas_infct, mcreplicate(10000, randomGwas(huri_g, infct), mc.cores = detectCores()))
random_gwas_infct_final <- as.numeric(random_gwas_infct)
random_gwas_infct_final[is.na(random_gwas_infct_final)] <- 0

######
# plot GWAS hits of 3 phenotypes permutation
source("gwas_phenotype_plot.r")

######
# 5. community detected
gwas_all_strong <- components(gwas_all_final, mode = "strong")
gwas_all_strong1 <- induced_subgraph(gwas_all_final, names(gwas_all_strong$membership[gwas_all_strong$membership == 1]))
gwas_all_strong1 <- simplify(gwas_all_strong1) # looks like doesn't hurt
gwas_all_ocg <- getOCG.clusters(as_data_frame(gwas_all_strong1), init.class.sys = 3, cent.class.sys = 0)
## visualized with whole network
pdf('/tmp/gwas_ocg.pdf')
plotOCGraph(gwas_all_ocg)
dev.off()

gwas_all_ocg_dif <- sort(gwas_all_ocg$clustsizes, decreasing = T)[1] / sort(gwas_all_ocg$clustsizes, decreasing = T)[2]
## randomized network, gwas all
all_ocg_rand_count <- list()
for (i in 1:1350) {
    tryCatch({
        all_ocg_rand_count[[i]] <- ocgRewireRatio(gwas_huri)
    }, error = function(e){
        print("Wired error!")
    })
}
all_ocg_sum <- as.numeric(unlist(lapply(all_ocg_rand_count, function(x) x[1])))[c(1:1000)]
all_ocg_sum3 <- as.numeric(unlist(lapply(all_ocg_rand_count, function(x) x[3])))[c(1:1000)]
## visualized with boxplot
boxplot(all_ocg_sum, horizontal = TRUE, cex.axis = 2)
points(gwas_all_ocg_dif, 1, pch = 19, col = "red", cex = 1.5)
## visualized with histogram
pdf("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/random_GWAS_commuRatio_v2.pdf", width = 3, height = 3)
dens <- hist(all_ocg_sum, breaks = 15, plot = FALSE)
par(mgp = c(2, 0.7, 0), ps = 8)
plot(dens, col = rgb(0.75, 0.75, 0.75, 1/2), freq = FALSE, border = NA, las = 1, xlim = c(0, 20), xlab = "Ratio of communities", ylab = "Frequency density", main = "", cex.sub = 0.5)
arrows(gwas_all_ocg_dif + 0.5, 0.1, gwas_all_ocg_dif + 0.5, 0.02, col = "#922687", lwd = 2, length = 0.1)
text(gwas_all_ocg_dif - 2, 0.12, paste0("observed = ", gwas_all_ocg_dif, "\np = ", table(all_ocg_sum >= gwas_all_ocg_dif)["TRUE"]/1000, "/", table(all_ocg_sum >= gwas_all_ocg_dif)["FALSE"]/1000), cex = 0.4, pos = 4)
dev.off()
## randomized network to reveal viral target significance


# 5.1.  based on 3 different phenotype
ctcl_network <- combineNetwork(huri_g, gwas_ctcl_candidate_huri)
ctcl_strong <- components(ctcl_network, mode = "strong")
ctcl_network_sub <- induced_subgraph(ctcl_network, names(ctcl_strong$membership[ctcl_strong[1]$membership == 1]))
hosp_network <- combineNetwork(huri_g, hosp)
hosp_strong <- components(hosp_network, mode = "strong")
hosp_network_sub <- induced_subgraph(hosp_network, names(hosp_strong$membership[hosp_strong[1]$membership == 1]))
infct_network <- combineNetwork(huri_g, infct)
infct_strong <- components(infct_network, mode = "strong")
infct_network_sub <- induced_subgraph(infct_network, names(infct_strong$membership[infct_strong[1]$membership == 1]))

ctcl_ocg <- getOCG.clusters(as_data_frame(ctcl_network_sub), init.class.sys = 1, cent.class.sys = 0)
hosp_ocg <- getOCG.clusters(as_data_frame(hosp_network_sub), init.class.sys = 1, cent.class.sys = 0)
infct_ocg <- getOCG.clusters(as_data_frame(infct_network_sub), init.class.sys = 1, cent.class.sys = 0)
######
# plotting graph
plotOCGraph()
plotOCGraph(ctcl_ocg)
plotOCGraph(hosp_ocg)
plotOCGraph(infct_ocg)
plot(gwas_all_final, vertex.size = 3, vertex.label = NA)
plot(ctcl_network, vertex.size = 3, vertex.label = NA)
plot(hosp_network, vertex.size = 3, vertex.label = NA)
plot(infct_network, vertex.size = 3, vertex.label = NA)
