# randomly choose proteins and have the average shortest path from either HuRI or BioPlex3
# Lin Chung-wen
# 26.08.2021

######
# load packages
library(igraph)
library(openxlsx)
library(rethinking)

mean_dist <- function(network, node, verbose = FALSE) {
	distance <- distances(network, node, node)
	if (verbose == TRUE) {
		print(distance)
	}
	return(mean(distance[lower.tri(distance)]))
}

avg_path_rand <- function(network, node, verbose = FALSE) {
	sample_node <- sample(V(network)$name, length(node))
	sample_avg <- mean_dist(network, sample_node)
	return(sample_avg)
}

to_Plot <- function(x, title, xlabel, ymax, ...) {
	dens_path <- hist(x, plot = FALSE, right = FALSE)
	plot(dens_path, col = rgb(0.75, 0.75, 0.75, 1/2), border = NA, las = 1, freq = FALSE, xlab = "", ylab = "Frequency", yaxt = "n", main = "", cex.sub = 0.5)
    mytitle <- paste0("Average shortest path\nviral targets in ", title)
    mtext(side = 3, line = 0.5, cex = 0.5, mytitle)
    mtext(side = 1, line = 1, cex = 0.3, xlabel)
    axis(side = 1, at = seq(0, xmax, by = 5) + 0.5, labels = seq(0, xmax, by = 5))
    arrows(length + 0.5, y1, length + 0.5, 0, col = "#922687", lwd = 2, length = 0.1)
    text(median(value) + 4, max(dens_gwas$counts / 10000), paste0("median = ", median(value)), col = "grey", cex = 0.5)
    text(length - 2, y2, paste0("observed = ", length, "\np = ", table(value >= length)["TRUE"]/10000), cex = 0.4, pos = 4)

}

plotDistance <- function(value, ymax, observe, phenotype) {
    dens_gwas <- hist(value, plot = FALSE, right = FALSE)
    plot(dens_gwas, col = rgb(0.75, 0.75, 0.75, 1/2), border = NA, las = 1, yaxt = "n", xlab = "Average shortest path", main = "", cex.sub = 0.5)
    mtext(side = 3, line = 1, cex = 1, phenotype)
    axis(side = 2, at = seq(0, ymax, by = 500), labels = seq(0, ymax/10000, by = 0.05), las = 1)
    arrows(observe, 500, observe, 0, col = "#922687", lwd = 2, length = 0.1)
    text(median(value), max(dens_gwas$counts), paste0("median = ", round(median(value), 4)), col = "grey", cex = 0.5)
    text(observe - 0.3, 700, paste0("observed = ", round(observe, 3), "\np = ", round(table(value >= observe)["FALSE"]/length(value[!is.infinite(value)]), 4)), cex = 0.4, pos = 4)
}

permut <- function(network, node_list,dist_list, times) {
	results <- list()
	pdf(paste0(names(networks[network]), "_avg_path.pdf"), width = 3, height = 3)	
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

		plotDistance(result, 3000, dist_list[[names(node_list[no])]], names(node_list[no]))

		results[[names(node_list[no])]] <- result
	}	
	dev.off()

	return(results)
}

permut_onlyPlot <- function(result_list, network, node_list,dist_list) {
	pdf(paste0(names(networks[network]), "_avg_path.pdf"), width = 4, height = 4)	
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

		plotDistance(sample_notInf, 3000, dist_list[[names(node_list[no])]], names(node_list[no]))
	}	
	dev.off()
}

######
# load data
huri <- read.csv("~/Documents/INET-work/references/HuRI_binaryPPI/HuRI_Tong_withSymbol.csv", header = T)
bioplex <- read.delim("/Volumes/GoogleDrive/My Drive/VirHostome_CW/GitHub/data/extended_table/BioPlex.3.0_edge.tsv", header = T)

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
	HuRI = simplify(graph_from_data_frame(huri[, c(5:6)], directed = FALSE), remove.loops = TRUE),
	BioPlex = simplify(graph_from_data_frame(bioplex[, c(5:6)], directed = FALSE), remove.loops = TRUE)
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
	infct_noParalog = c("LZTFL1", "NXPE3", "OAS1", "PLEKHA4", "SLC6A20")
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
# randomly select the same number proteins as node list from 1) HuRI and 2) BioPlex, calculating the average shortest path is significant than the observation (as node list) or not

# HuRI
huri_result <- permut("HuRI", nodes_huri, dist_huri, 10000)

# Node: gwas_2021a; Network: HuRI; avg shortest path: 3.1
# [1] "right tail: FALSE; median = 3.8; sig = 0.0456"
# Node: gwas_2021a_noParalog; Network: HuRI; avg shortest path: 2.833
# [1] "right tail: FALSE; median = 3.833; sig = 0.0194"
# Node: gwas_2021b; Network: HuRI; avg shortest path: 3.515
# [1] "right tail: FALSE; median = 3.838; sig = 0.083"
# Node: gwas_2021b_noParalog; Network: HuRI; avg shortest path: 3.571
# [1] "right tail: FALSE; median = 3.824; sig = 0.1634"
# Node: ctcl; Network: HuRI; avg shortest path: 3.578
# [1] "right tail: FALSE; median = 3.844; sig = 0.2119"
# Node: ctcl_noParalog; Network: HuRI; avg shortest path: 3.762
# [1] "right tail: FALSE; median = 3.81; sig = 0.4172"
# Node: hosp; Network: HuRI; avg shortest path: 3.41
# [1] "right tail: FALSE; median = 3.833; sig = 0.0629"
# Node: hosp_noParalog; Network: HuRI; avg shortest path: 3.467
# [1] "right tail: FALSE; median = 3.822; sig = 0.1235"
# Node: infct; Network: HuRI; avg shortest path: 3.867
# [1] "right tail: FALSE; median = 3.8; sig = 0.5183"
# Node: infct_noParalog; Network: HuRI; avg shortest path: 3.8
# [1] "right tail: FALSE; median = 3.8; sig = 0.4469"

# replot, if necessary
# permut_onlyPlot(huri_result, "HuRI", nodes_huri, dist_huri)

# BioPlex3
bioplex_result <- permut("BioPlex", nodes_bioplex, dist_bioplex, 10000)
# Node: gwas_2021a; Network: BioPlex; avg shortest path: 3.952
# [1] "right tail: FALSE; median = 3.81; sig = 0.6252"
# Node: gwas_2021b; Network: BioPlex; avg shortest path: 3.732
# [1] "right tail: FALSE; median = 3.837; sig = 0.2681"
# Node: gwas_2021b_noParalog; Network: BioPlex; avg shortest path: 3.757
# [1] "right tail: FALSE; median = 3.838; sig = 0.3296"
# Node: ctcl; Network: BioPlex; avg shortest path: 3.689
# [1] "right tail: FALSE; median = 3.822; sig = 0.2889"
# Node: ctcl_noParalog; Network: BioPlex; avg shortest path: 3.857
# [1] "right tail: FALSE; median = 3.81; sig = 0.5084"
# Node: hosp; Network: BioPlex; avg shortest path: 3.882
# [1] "right tail: FALSE; median = 3.838; sig = 0.5796"
# Node: hosp_noParalog; Network: BioPlex; avg shortest path: 3.967
# [1] "right tail: FALSE; median = 3.835; sig = 0.6902"
# Node: infct; Network: BioPlex; avg shortest path: 3.511
# [1] "right tail: FALSE; median = 3.844; sig = 0.1007"

# replot, if necessary
# permut_onlyPlot(bioplex_result, "BioPlex", nodes_bioplex, dist_bioplex)

######
# Topological twins, 
# Pascal: repeat the network randomized shortest path with 4 or 5 different proteins that have the same number of interactions
# 

# 1. degree calculation
deg_nodes_huri <- list(
	gwas_2021a = data.frame(degree = degree(huri_g, nodes_huri[[1]])),
	gwas_2021a_noParalog = data.frame(degree = degree(huri_g, nodes_huri[[2]])),
	gwas_2021b = data.frame(degree = degree(huri_g, nodes_huri[[3]])),
	gwas_2021b_noParalog = data.frame(degree = degree(huri_g, nodes_huri[[4]])),
	ctcl = data.frame(degree = degree(huri_g, nodes_huri[[5]])),
	ctcl_noParalog = data.frame(degree = degree(huri_g, nodes_huri[[6]]))
	)

deg_nodes_bioplex <- list(
	gwas_2021a = data.frame(degree = degree(bioplex_g, nodes_bioplex[[1]])),
	gwas_2021b = data.frame(degree = degree(bioplex_g, nodes_bioplex[[2]])),
	gwas_2021b_noParalog = data.frame(degree = degree(bioplex_g, nodes_bioplex[[3]])),
	ctcl = data.frame(degree = degree(bioplex_g, nodes_bioplex[[4]])),
	ctcl_noParalog = data.frame(degree = degree(bioplex_g, nodes_bioplex[[5]]))
	)

######
# save result
save.image("~/Documents/INET-work/virus_network/statistic_results/GWAS/random_protein_distance.RData")