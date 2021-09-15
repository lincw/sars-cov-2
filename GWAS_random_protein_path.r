# randomly choose proteins and have the average shortest path from either HuRI or BioPlex3
# Lin Chung-wen
# 26.08.2021
# 31.08.2021
# 06.09.2021, permute function not used now!
######
# load packages
library(igraph)
library(openxlsx)
library(rethinking)
library(data.table)
mean_dist <- function(network, node, verbose = FALSE) {
	distance <- distances(network, node, node)
	if (verbose == TRUE) {
		print(distance)
	}
	return(mean(distance[lower.tri(distance)]))
}
source("~/Documents/INET-work/virus_network/src/combineNetwork.r")
sqrt_deg_lo <- function(network, node) {
	deg <- degree(networks[[network]], node)
	# deg_out <- deg - round(deg ^ 0.5, 0)
	deg_out <- deg - round(deg * 0.9, 0)
}
sqrt_deg_hi <- function(network, node) {
	deg <- degree(networks[[network]], node)
	# deg_out <- deg + round(deg ^ 0.5, 0)
	deg_out <- deg + round(deg * 1.1, 0)
}
viralTarget <- function(network, node, viral_list) {
    subnet <- combineNetwork(network, node)
    ls <- V(subnet)$name[V(subnet)$name %in% viral_list]
    return(ls)
}
dens <- function(network, node) {
    net <- combineNetwork(network, node)
    edge_dens <- edge_density(net)
    return(edge_dens)
}

avg_path_rand <- function(network, node, verbose = FALSE) {
	sample_node <- sample(V(network)$name, length(node))
	sample_avg <- mean_dist(network, sample_node)
	return(sample_avg)
}

plotDistance <- function(value, observe, phenotype, xlabel) {
    dens_gwas <- hist(value, plot = FALSE, right = FALSE)
    ymax <- round(max(dens_gwas$count)/1000, 1) * 1000
    plot(dens_gwas, col = rgb(0.75, 0.75, 0.75, 1/2), border = NA, las = 1, yaxt = "n", xlab = xlabel, main = "", cex.sub = 0.5)
    mtext(side = 3, line = 1, cex = 1, phenotype)
    axis(side = 2, at = seq(0, ymax, by = 500), labels = seq(0, ymax/10000, by = 0.05), las = 1)
    arrows(observe, 500, observe, 0, col = "#922687", lwd = 2, length = 0.1)
    text(median(value), max(dens_gwas$counts), paste0("median = ", round(median(value), 4)), col = "grey", cex = 0.5)
    text(observe, 700, paste0("observed = ", round(observe, 3), "\np = ", round(table(value >= observe)["FALSE"]/length(value[!is.infinite(value)]), 4)), cex = 0.4, pos = 4)
}

permut <- function(network, node_list, dist_list, times, pdf_extension) {
	results <- list()
	pdf(paste0(names(networks[network]), "_", pdf_extension, ".pdf"), width = 3, height = 3)
	for (no in 1:length(node_list)) {
		right <- "FALSE"
		result <- c()
		result <- c(result, mcreplicate(times, avg_path_rand(networks[[network]], node_list[[no]]), mc.cores = detectCores()))

		sample_notInf <- result[!is.infinite(result)]

		sig <- table(sample_notInf >= dist_list[[names(node_list[no])]])[right] / length(sample_notInf)
		if (is.na(sig)) {
			sig_value <- paste0(" < ", 1/length(sample_notInf))
		} else {
			sig_value <- round(as.numeric(sig), 4)
		}
		cat(paste0("Node: ", names(node_list[no]), "; Network: ", names(networks[network]), "; avg shortest path: ", round(dist_list[[names(node_list[no])]], 3), "\n"))
		print(paste0("right tail: ", right, "; median = ", round(median(sample_notInf), 3), "; sig = ", sig_value))

		plotDistance(result, 3000, dist_list[[names(node_list[no])]], names(node_list[no]), pdf_extension)

		results[[names(node_list[no])]] <- result
	}
	dev.off()

	return(results)
}

permut_onlyPlot <- function(result_list, network, node_list, dist_list, pdf_extension) {
	pdf(paste0(names(networks[network]), "_", pdf_extension, ".pdf"), width = 4, height = 4)
	for (no in 1:length(node_list)) {
		right <- "FALSE"
		result <- result_list[[names(node_list[no])]]

		sample_notInf <- result[!is.infinite(result)]

		sig <- table(sample_notInf >= dist_list[[names(node_list[no])]])[right] / length(sample_notInf)
		if (is.na(sig)) {
			sig_value <- paste0(" < ", 1/length(sample_notInf))
		} else {
			sig_value <- round(as.numeric(sig), 4)
		}

		plotDistance(sample_notInf, dist_list[[names(node_list[no])]], names(node_list[no]), pdf_extension)
	}
	dev.off()
}

rand_distance <- function(times, network, nodes_list, ...) {
	dist_result <- list()
	dist_name <- names(nodes_list)
	pb <- txtProgressBar(0, times, style = 3)
	for (i in 1:times) {
	    re <- rewire(networks[[network]], keeping_degseq(niter = gsize(networks[[network]]) * 10))
		for (n in dist_name) {
			dist_result[[n]] <- c(dist_result[[n]], mean_dist(re, nodes_list[[n]]$sample))
		}
    	setTxtProgressBar(pb, i)
    	Sys.sleep(time = 0)
		}
	close(pb)
	return(dist_result)
}

######
# load data
huri <- read.csv("~/Documents/INET-work/references/HuRI_binaryPPI/HuRI_Tong_withSymbol.csv", header = T)
bioplex <- read.delim("/Volumes/GoogleDrive/My Drive/VirHostome_CW/GitHub/data/extended_table/BioPlex.3.0_edge.tsv", header = T)

husci <- read.csv("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/binary_node_1126.csv", header = T)
husci_node <- unique(husci[husci$group == "human", "node"])
# as documentary
gwas2021a <- read.csv("~/Documents/INET-work/virus_network/references/GWAS/Genetic\ mechanisms\ of\ critical\ illness\ in\ COVID-19/table1.csv", header = T)
gwas2021b <- read.xlsx("~/Documents/INET-work/virus_network/statistic_results/GWAS/COVID_GWAS_hit_inHUSCI_v2.xlsx")

gwas2021a <- c(gwas2021a$Locus[c(1:4, 6:8)], "OAS1", "OAS2", "OAS3")
gwas2021b_ctcl <- gwas2021b[gwas2021b$COVID_ctcl == 1 & !is.na(gwas2021b$COVID_ctcl), "All.LD"]
gwas2021b_hosp <- gwas2021b[gwas2021b$COV_hosp == 1 & !is.na(gwas2021b$COV_hosp), "All.LD"]
gwas2021b_infct <- gwas2021b[gwas2021b$COVID_infct == 1 & !is.na(gwas2021b$COVID_infct), "All.LD"]

######
# have network
networks <- list(
	HuRI = igraph::simplify(graph_from_data_frame(huri[, c(5:6)], directed = FALSE), remove.loops = FALSE),
	BioPlex = igraph::simplify(graph_from_data_frame(bioplex[, c(5:6)], directed = FALSE), remove.loops = FALSE),
	HuRI_noLoop = igraph::simplify(graph_from_data_frame(huri[, c(5:6)], directed = FALSE), remove.loops = TRUE),
	BioPlex_noLoop = igraph::simplify(graph_from_data_frame(bioplex[, c(5:6)], directed = FALSE), remove.loops = TRUE)
	)

######
# have protein list
nodes_huri <- list(
	gwas_2021a = gwas2021a[gwas2021a %in% V(networks[["HuRI"]])$name],
	gwas_2021a_noParalog = c("LZTFL1", "CCHCR1", "TYK2", "OAS1"),
	gwas_2021b = gwas2021b$All.LD[gwas2021b$All.LD %in% V(networks[["HuRI"]])$name],
	gwas_2021b_noParalog = c("FDX2", "ICAM3", "KAT7", "LZTFL1", "NSF", "NXPE3", "OAS1", "PLEKHA4", "RAVER1", "SLC6A20", "SPPL2C", "STH", "TMEM65", "TYK2"),
	ctcl = gwas2021b_ctcl[gwas2021b_ctcl %in% V(networks[["HuRI"]])$name],
	ctcl_noParalog = c("FDX2", "ICAM3", "KAT7", "LZTFL1", "OAS1", "RAVER1", "TYK2"),
	hosp = gwas2021b_hosp[gwas2021b_hosp %in% V(networks[["HuRI"]])$name],
	hosp_noParalog = c("FDX2", "ICAM3", "LZTFL1", "NSF", "OAS1", "RAVER1", "SPPL2C", "STH", "TMEM65", "TYK2"),
	infct = gwas2021b_infct[gwas2021b_infct %in% V(networks[["HuRI"]])$name],
	infct_noParalog = c("LZTFL1", "NXPE3", "OAS1", "PLEKHA4", "SLC6A20"),
	gwas_union = unique(union(gwas2021a, gwas2021b$All.LD))[unique(union(gwas2021a, gwas2021b$All.LD)) %in% V(networks[["HuRI"]])$name]
)

nodes_bioplex <- list(
	gwas_2021a = gwas2021a[gwas2021a %in% V(networks[["BioPlex"]])$name],
	gwas_2021b = gwas2021b$All.LD[gwas2021b$All.LD %in% V(networks[["BioPlex"]])$name],
	gwas_2021b_noParalog = c("ARHGAP27", "ARL17A", "CEP97", "DPP9", "FOXP4", "ICAM3", "IFNAR2", "LRRC37A2", "LZTFL1", "MAPT", "NUCB1", "OAS2", "PPP1R15A", "RAVER1", "RPL24", "SLC6A20", "TMEM65", "TULP2", "TYK2", "WNT3", "ZBTB11"),
	ctcl = gwas2021b_ctcl[gwas2021b_ctcl %in% V(networks[["BioPlex"]])$name],
	ctcl_noParalog = c("DPP9", "ICAM3", "IFNAR2", "LZTFL1", "OAS2", "RAVER1", "TYK2"),
	hosp = gwas2021b_hosp[gwas2021b_hosp %in% V(networks[["BioPlex"]])$name],
	hosp_noParalog = c("ARHGAP27", "ARL17A", "DPP9", "FOXP4", "ICAM3", "IFNAR2", "LRRC37A2", "LZTFL1", "MAPT", "OAS2", "RAVER1", "TMEM65", "TYK2", "WNT3"),
	infct = gwas2021b_infct[gwas2021b_infct %in% V(networks[["BioPlex"]])$name]
)

######
# degree calculation
deg_huri <- degree(networks[["HuRI"]])
deg_nodes_huri <- list(
	gwas_2021a = data.frame(
		gene = nodes_huri[["gwas_2021a"]],
		degree = degree(networks[["HuRI"]], nodes_huri[["gwas_2021a"]]),
		low = sqrt_deg_lo("HuRI", nodes_huri[["gwas_2021a"]]),
		high = sqrt_deg_hi("HuRI", nodes_huri[["gwas_2021a"]])),
	gwas_2021a_noParalog = data.frame(
		gene = nodes_huri[["gwas_2021a_noParalog"]],
		degree = degree(networks[["HuRI"]], nodes_huri[["gwas_2021a_noParalog"]]),
		low = sqrt_deg_lo("HuRI", nodes_huri[["gwas_2021a_noParalog"]]),
		high = sqrt_deg_hi("HuRI", nodes_huri[["gwas_2021a_noParalog"]])
		),
	gwas_2021b = data.frame(
		gene = nodes_huri[["gwas_2021b"]],
		degree = degree(networks[["HuRI"]], nodes_huri[["gwas_2021b"]]),
		low = sqrt_deg_lo("HuRI", nodes_huri[["gwas_2021b"]]),
		high = sqrt_deg_hi("HuRI", nodes_huri[["gwas_2021b"]])
		),
	gwas_2021b_noParalog = data.frame(
		gene = nodes_huri[["gwas_2021b_noParalog"]],
		degree = degree(networks[["HuRI"]], nodes_huri[["gwas_2021b_noParalog"]]),
		low = sqrt_deg_lo("HuRI", nodes_huri[["gwas_2021b_noParalog"]]),
		high = sqrt_deg_hi("HuRI", nodes_huri[["gwas_2021b_noParalog"]])
		),
	ctcl = data.frame(
		gene = nodes_huri[["ctcl"]],
		degree = degree(networks[["HuRI"]], nodes_huri[["ctcl"]]),
		low = sqrt_deg_lo("HuRI", nodes_huri[["ctcl"]]),
		high = sqrt_deg_hi("HuRI", nodes_huri[["ctcl"]])
		),
	ctcl_noParalog = data.frame(
		gene = nodes_huri[["ctcl_noParalog"]],
		degree = degree(networks[["HuRI"]], nodes_huri[["ctcl_noParalog"]]),
		low = sqrt_deg_lo("HuRI", nodes_huri[["ctcl_noParalog"]]),
		high = sqrt_deg_hi("HuRI", nodes_huri[["ctcl_noParalog"]])
		),
	hosp = data.frame(
		gene = nodes_huri[["hosp"]],
		degree = degree(networks[["HuRI"]], nodes_huri[["hosp"]]),
		low = sqrt_deg_lo("HuRI", nodes_huri[["hosp"]]),
		high = sqrt_deg_hi("HuRI", nodes_huri[["hosp"]])
		),
	hosp_noParalog = data.frame(
		gene = nodes_huri[["hosp_noParalog"]],
		degree = degree(networks[["HuRI"]], nodes_huri[["hosp_noParalog"]]),
		low = sqrt_deg_lo("HuRI", nodes_huri[["hosp_noParalog"]]),
		high = sqrt_deg_hi("HuRI", nodes_huri[["hosp_noParalog"]])
		),
	infct = data.frame(
		gene = nodes_huri[["infct"]],
		degree = degree(networks[["HuRI"]], nodes_huri[["infct"]]),
		low = sqrt_deg_lo("HuRI", nodes_huri[["infct"]]),
		high = sqrt_deg_hi("HuRI", nodes_huri[["infct"]])
		),
	infct_noParalog = data.frame(
		gene = nodes_huri[["infct_noParalog"]],
		degree = degree(networks[["HuRI"]], nodes_huri[["infct_noParalog"]]),
		low = sqrt_deg_lo("HuRI", nodes_huri[["infct_noParalog"]]),
		high = sqrt_deg_hi("HuRI", nodes_huri[["infct_noParalog"]])
		)
	)

# average shortest path of individual GWAS protein list
dist_huri <- list(
	gwas_2021a = mean_dist(networks[["HuRI"]], nodes_huri[["gwas_2021a"]]),
	gwas_2021a_noParalog = mean_dist(networks[["HuRI"]], nodes_huri[["gwas_2021a_noParalog"]]),
	gwas_2021b = mean_dist(networks[["HuRI"]], nodes_huri[["gwas_2021b"]]),
	gwas_2021b_noParalog = mean_dist(networks[["HuRI"]], nodes_huri[["gwas_2021b_noParalog"]]),
	ctcl = mean_dist(networks[["HuRI"]], nodes_huri[["ctcl"]]),
	ctcl_noParalog = mean_dist(networks[["HuRI"]], nodes_huri[["ctcl_noParalog"]]),
	hosp = mean_dist(networks[["HuRI"]], nodes_huri[["hosp"]]),
	hosp_noParalog = mean_dist(networks[["HuRI"]], nodes_huri[["hosp_noParalog"]]),
	infct = mean_dist(networks[["HuRI"]], nodes_huri[["infct"]]),
	infct_noParalog = mean_dist(networks[["HuRI"]], nodes_huri[["infct_noParalog"]])
	)

dist_bioplex <- list(
	gwas_2021a = mean_dist(networks[["BioPlex"]], nodes_bioplex[["gwas_2021a"]]),
	gwas_2021b = mean_dist(networks[["BioPlex"]], nodes_bioplex[["gwas_2021b"]]),
	gwas_2021b_noParalog = mean_dist(networks[["BioPlex"]], nodes_bioplex[["gwas_2021b_noParalog"]]),
	ctcl = mean_dist(networks[["BioPlex"]], nodes_bioplex[["ctcl"]]),
	ctcl_noParalog = mean_dist(networks[["BioPlex"]], nodes_bioplex[["ctcl_noParalog"]]),
	hosp = mean_dist(networks[["BioPlex"]], nodes_bioplex[["hosp"]]),
	hosp_noParalog = mean_dist(networks[["BioPlex"]], nodes_bioplex[["hosp_noParalog"]]),
	infct = mean_dist(networks[["BioPlex"]], nodes_bioplex[["infct"]])
	)

######
# viral target number in individual GWAS protein list generated GWAS subnetwork
viral_huri <- list(
    gwas_2021a = viralTarget(networks[["HuRI"]], nodes_huri[["gwas_2021a"]], husci_node),
    gwas_2021a_noParalog = viralTarget(networks[["HuRI"]], nodes_huri[["gwas_2021a_noParalog"]], husci_node),
    gwas_2021b = viralTarget(networks[["HuRI"]], nodes_huri[["gwas_2021b"]], husci_node),
    gwas_2021b_noParalog = viralTarget(networks[["HuRI"]], nodes_huri[["gwas_2021b_noParalog"]], husci_node),
    ctcl = viralTarget(networks[["HuRI"]], nodes_huri[["ctcl"]], husci_node),
    ctcl_noParalog = viralTarget(networks[["HuRI"]], nodes_huri[["ctcl_noParalog"]], husci_node),
    hosp = viralTarget(networks[["HuRI"]], nodes_huri[["hosp"]], husci_node),
    hosp_noParalog = viralTarget(networks[["HuRI"]], nodes_huri[["hosp_noParalog"]], husci_node),
    infct = viralTarget(networks[["HuRI"]], nodes_huri[["infct"]], husci_node),
    infct_noParalog = viralTarget(networks[["HuRI"]], nodes_huri[["infct_noParalog"]], husci_node)
)
viral_huri_length <- lapply(viral_huri, length)
viral_bioplex <- list(
    gwas_2021a = viralTarget(networks[["BioPlex"]], nodes_bioplex[["gwas_2021a"]], husci_node),
    gwas_2021b = viralTarget(networks[["BioPlex"]], nodes_bioplex[["gwas_2021b"]], husci_node),
    ctcl = viralTarget(networks[["BioPlex"]], nodes_bioplex[["ctcl"]], husci_node),
    hosp = viralTarget(networks[["BioPlex"]], nodes_bioplex[["hosp"]], husci_node),
    infct = viralTarget(networks[["BioPlex"]], nodes_bioplex[["infct"]], husci_node)
)
viral_bioplex_length <- lapply(viral_bioplex, length)

######
# subnetwork density of individual GWAS protein list
dens_huri <- list(
    gwas_2021a = dens(networks[["HuRI"]], nodes_huri[["gwas_2021a"]]),
    gwas_2021a_noParalog = dens(networks[["HuRI"]], nodes_huri[["gwas_2021a_noParalog"]]),
    gwas_2021b = dens(networks[["HuRI"]], nodes_huri[["gwas_2021b"]]),
    gwas_2021b_noParalog = dens(networks[["HuRI"]], nodes_huri[["gwas_2021b_noParalog"]]),
    ctcl = dens(networks[["HuRI"]], nodes_huri[["ctcl"]]),
    ctcl_noParalog = dens(networks[["HuRI"]], nodes_huri[["ctcl_noParalog"]]),
    hosp = dens(networks[["HuRI"]], nodes_huri[["hosp"]]),
    hosp_noParalog = dens(networks[["HuRI"]], nodes_huri[["hosp_noParalog"]]),
    infct = dens(networks[["HuRI"]], nodes_huri[["infct"]]),
    infct_noParalog = dens(networks[["HuRI"]], nodes_huri[["infct_noParalog"]])
)
dens_bioplex <- list(
    gwas_2021a = dens(networks[["BioPlex"]], nodes_bioplex[["gwas_2021a"]]),
    gwas_2021b = dens(networks[["BioPlex"]], nodes_bioplex[["gwas_2021b"]]),
    ctcl = dens(networks[["BioPlex"]], nodes_bioplex[["ctcl"]]),
    hosp = dens(networks[["BioPlex"]], nodes_bioplex[["hosp"]]),
    infct = dens(networks[["BioPlex"]], nodes_bioplex[["infct"]])
)
######
# subnetwork interactions (gsize) of individual GWAS protein list
interact_huri <- list(
    gwas_2021a = gsize(combineNetwork(networks[["HuRI"]], nodes_huri[["gwas_2021a"]])),
    gwas_2021a_noParalog = gsize(combineNetwork(networks[["HuRI"]], nodes_huri[["gwas_2021a_noParalog"]])),
    gwas_2021b = gsize(combineNetwork(networks[["HuRI"]], nodes_huri[["gwas_2021b"]])),
    gwas_2021b_noParalog = gsize(combineNetwork(networks[["HuRI"]], nodes_huri[["gwas_2021b_noParalog"]])),
    ctcl = gsize(combineNetwork(networks[["HuRI"]], nodes_huri[["ctcl"]])),
    ctcl_noParalog = gsize(combineNetwork(networks[["HuRI"]], nodes_huri[["ctcl_noParalog"]])),
    hosp = gsize(combineNetwork(networks[["HuRI"]], nodes_huri[["hosp"]])),
    hosp_noParalog = gsize(combineNetwork(networks[["HuRI"]], nodes_huri[["hosp_noParalog"]])),
    infct = gsize(combineNetwork(networks[["HuRI"]], nodes_huri[["infct"]])),
    infct_noParalog = gsize(combineNetwork(networks[["HuRI"]], nodes_huri[["infct_noParalog"]]))
)
interact_bioplex <- list(
    gwas_2021a = gsize(combineNetwork(networks[["BioPlex"]], nodes_bioplex[["gwas_2021a"]])),
    gwas_2021b = gsize(combineNetwork(networks[["BioPlex"]], nodes_bioplex[["gwas_2021b"]])),
    ctcl = gsize(combineNetwork(networks[["BioPlex"]], nodes_bioplex[["ctcl"]])),
    hosp = gsize(combineNetwork(networks[["BioPlex"]], nodes_bioplex[["hosp"]])),
    infct = gsize(combineNetwork(networks[["BioPlex"]], nodes_bioplex[["infct"]]))
)

######
# a. randomly select the same number proteins as node list from 1) HuRI and 2) BioPlex, calculating the average shortest path is significant than the observation (as node list) or not

# HuRI, using function
huri_result <- permut("HuRI", nodes_huri, dist_huri, 10000)

############################
# Node: gwas_2021a; Network: HuRI; avg shortest path: 3.1
# [1] "right tail: FALSE; median = 3.8; sig = 0.0455"
# Node: gwas_2021a_noParalog; Network: HuRI; avg shortest path: 2.833
# [1] "right tail: FALSE; median = 3.833; sig = 0.0182"
# Node: gwas_2021b; Network: HuRI; avg shortest path: 3.515
# [1] "right tail: FALSE; median = 3.838; sig = 0.0855"
# Node: gwas_2021b_noParalog; Network: HuRI; avg shortest path: 3.571
# [1] "right tail: FALSE; median = 3.824; sig = 0.1595"
# Node: ctcl; Network: HuRI; avg shortest path: 3.578
# [1] "right tail: FALSE; median = 3.844; sig = 0.2052"
# Node: ctcl_noParalog; Network: HuRI; avg shortest path: 3.762
# [1] "right tail: FALSE; median = 3.81; sig = 0.4199"
# Node: hosp; Network: HuRI; avg shortest path: 3.41
# [1] "right tail: FALSE; median = 3.833; sig = 0.0588"
# Node: hosp_noParalog; Network: HuRI; avg shortest path: 3.467
# [1] "right tail: FALSE; median = 3.822; sig = 0.1205"
# Node: infct; Network: HuRI; avg shortest path: 3.867
# [1] "right tail: FALSE; median = 3.8; sig = 0.5178"
# Node: infct_noParalog; Network: HuRI; avg shortest path: 3.8
# [1] "right tail: FALSE; median = 3.8; sig = 0.4499"

## replot, if necessary
# permut_onlyPlot(huri_result, "HuRI", nodes_huri, dist_huri)

# HuRI, using for loop
times <- 10000
pb <- txtProgressBar(0, times, style = 3)
protein_rand_list <- list()
meandist_rand_result <- list()
viral_target_rand_result <- list()
density_rand_result <- list()
interaction_rand_result <- list()
for (a in 1:times) {
	dn_huri_sam <- deg_nodes_huri
	for (i in 1:length(dn_huri_sam)) {
		protein_list <- sample(names(deg_huri), nrow(dn_huri_sam[[i]]), replace = FALSE)
		meandist_rand_result[[names(dn_huri_sam[i])]] <- c(meandist_rand_result[[names(dn_huri_sam[i])]], mean_dist(networks[["HuRI"]], protein_list))

        rand_network <- combineNetwork(networks[["HuRI"]], protein_list)
        protein_rand_list[[names(dn_huri_sam[i])]][[a]] <- protein_list

        viral_target_rand_result[[names(dn_huri_sam[i])]] <- c(viral_target_rand_result[[names(dn_huri_sam[i])]], as.numeric(table(V(rand_network)$name %in% husci_node)["TRUE"]))

        density_rand_result[[names(dn_huri_sam[i])]] <- c(density_rand_result[[names(dn_huri_sam[i])]], edge_density(rand_network))

        interaction_rand_result[[names(dn_huri_sam[i])]] <- c(interaction_rand_result[[names(dn_huri_sam[i])]], gsize(rand_network))
	}
    setTxtProgressBar(pb, a)
    Sys.sleep(time = 0)
}
close(pb)

viral_target_rand_result <- lapply(viral_target_rand_result, function(x) {
	x[is.na(x)] <- 0
	x
	})

permut_onlyPlot(meandist_rand_result, "HuRI", nodes_huri, dist_huri, "avg shortest path, random")
permut_onlyPlot(viral_target_rand_result, "HuRI", nodes_huri, viral_huri_length, "viral target, random")
permut_onlyPlot(density_rand_result, "HuRI", nodes_huri, dens_huri, "network density, random")
permut_onlyPlot(interaction_rand_result, "HuRI", nodes_huri, interact_huri, "interactions, random")

# BioPlex3
bioplex_result <- permut("BioPlex", nodes_bioplex, dist_bioplex, 10000)
# Node: gwas_2021a; Network: BioPlex; avg shortest path: 3.952
# [1] "right tail: FALSE; median = 3.81; sig = 0.6187"
# Node: gwas_2021b; Network: BioPlex; avg shortest path: 3.732
# [1] "right tail: FALSE; median = 3.837; sig = 0.2671"
# Node: gwas_2021b_noParalog; Network: BioPlex; avg shortest path: 3.757
# [1] "right tail: FALSE; median = 3.838; sig = 0.3334"
# Node: ctcl; Network: BioPlex; avg shortest path: 3.689
# [1] "right tail: FALSE; median = 3.822; sig = 0.2772"
# Node: ctcl_noParalog; Network: BioPlex; avg shortest path: 3.857
# [1] "right tail: FALSE; median = 3.81; sig = 0.5022"
# Node: hosp; Network: BioPlex; avg shortest path: 3.882
# [1] "right tail: FALSE; median = 3.838; sig = 0.5705"
# Node: hosp_noParalog; Network: BioPlex; avg shortest path: 3.967
# [1] "right tail: FALSE; median = 3.835; sig = 0.6982"
# Node: infct; Network: BioPlex; avg shortest path: 3.511
# [1] "right tail: FALSE; median = 3.844; sig = 0.1007"

## replot, if necessary
# permut_onlyPlot(bioplex_result, "BioPlex", nodes_bioplex, dist_bioplex)

# b. random select proteins with similar degree profile (degree of proteins plus/minus round of square of degree proteins by 1)

# b.1. proteins with similar degree sharing similar average shortest path in HuRI?
# random proteins selection, degree +- round(degree ^ 0.5, 1)
times <- 10000
pb <- txtProgressBar(0, times, style = 3)
meandist_sqrt_result <- list()
protein_result <- list()
viral_target_result <- list()
density_result <- list()
interaction_result <- list()
for (a in 1:times) {
	dn_huri_sam <- deg_nodes_huri
	for (i in 1:length(dn_huri_sam)) {
		for (j in 1:nrow(dn_huri_sam[[i]])) {
			dn_huri_sam[[i]][j, "sample"] <- sample(names(deg_huri[between(deg_huri, dn_huri_sam[[i]][j, 3], dn_huri_sam[[i]][j, 4])]), 1)
			# dn_huri_sam[[i]][j, "sample_degree"] <- degree(networks[["HuRI"]], dn_huri_sam[[i]][j, "sample"])
		}
		toCheck <- all(table(dn_huri_sam[[i]]$sample) == 1)
		protein_result[[names(dn_huri_sam[i])]][[a]] <- dn_huri_sam[[i]]$sample
		while (toCheck != TRUE) { # code error with `if`
			for (j in 1:nrow(dn_huri_sam[[i]])) {
				dn_huri_sam[[i]][j, "sample"] <- sample(names(deg_huri[between(deg_huri, dn_huri_sam[[i]][j, 3], dn_huri_sam[[i]][j, 4])]), 1)
			}
			toCheck <- all(table(dn_huri_sam[[i]]$sample) == 1)
			protein_result[[names(dn_huri_sam[i])]][[a]] <- dn_huri_sam[[i]]$sample
		}
		meandist_sqrt_result[[names(dn_huri_sam[i])]] <- c(meandist_sqrt_result[[names(dn_huri_sam[i])]], mean_dist(networks[["HuRI"]], dn_huri_sam[[i]]$sample))
        
        rand_network <- combineNetwork(networks[["HuRI"]], dn_huri_sam[[i]]$sample)
        
        viral_target_result[[names(dn_huri_sam[i])]] <- c(viral_target_result[[names(dn_huri_sam[i])]], as.numeric(table(V(rand_network)$name %in% husci_node)["TRUE"]))
        density_result[[names(dn_huri_sam[i])]] <- c(density_result[[names(dn_huri_sam[i])]], edge_density(rand_network))
        interaction_result[[names(dn_huri_sam[i])]] <- c(interaction_result[[names(dn_huri_sam[i])]], gsize(rand_network))
	}
    setTxtProgressBar(pb, a)
    Sys.sleep(time = 0)
}
close(pb)
viral_target_result <- lapply(viral_target_result, function(x) {
	x[is.na(x)] <- 0
	x
	})

permut_onlyPlot(meandist_sqrt_result, "HuRI", nodes_huri, dist_huri, "avg shortest path")
permut_onlyPlot(viral_target_result, "HuRI", nodes_huri, viral_huri_length, "viral target")
permut_onlyPlot(density_result, "HuRI", nodes_huri, dens_huri, "network density")
permut_onlyPlot(interaction_result, "HuRI", nodes_huri, interact_huri, "interactions")

######
# Topological twins,
# Pascal: repeat the network randomized shortest path with 4 or 5 different proteins that have the same number of interactions

# observed and random protein shortest path calculation, in degree preserving randomized network
rand_result_huri <- rand_distance(10000, "HuRI", deg_nodes_huri_sample)

unlist(obs_rand_huri)
unlist(lapply(mean_result_huri, function(x) mean(x[!is.infinite(x)])))

mean_result_huri <- list()
pb <- txtProgressBar(0, 1000, style = 3)
for (i in 1:1000) {
    re <- rewire(networks[["HuRI"]], keeping_degseq(niter = gsize(networks[["HuRI"]]) * 10))
	for (j in 1:length(deg_nodes_huri)) {
	    mean_dis <- mean_dist(re, deg_nodes_huri[[j]]$sample)
	    mean_result_huri[[names(deg_nodes_huri)[j]]] <- c(mean_result_huri[[names(deg_nodes_huri)[j]]], mean_dis)
	}
    setTxtProgressBar(pb, i, title = "mean distances of randomized network")
    Sys.sleep(time = 0)
}
close(pb)

means_value <- lapply(mean_result_huri, function(x) {
	x_new <- x[!is.infinite(x)]
	mean(x_new)
	})

permut_onlyPlot(means_value, "HuRI", nodes_huri, dist_huri)

######
# save result
save.image("~/Documents/INET-work/virus_network/statistic_results/GWAS/random_protein_distance.RData")
