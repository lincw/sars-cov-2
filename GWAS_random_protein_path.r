# randomly choose proteins and have the average shortest path from either HuRI or BioPlex3
# Lin Chung-wen
# 26.08.2021

######
# load packages
library(igraph)
library(openxlsx)
library(rethinking)

avg_path_rand <- function(network, node, verbose = FALSE) {
	sample_node <- sample(V(network)$name, length(node))
	sample_distance <- distances(network, v = sample_node, to = sample_node)
	sample_avg <- mean(sample_distance[lower.tri(sample_distance)])
	if (verbose == TRUE) {
		print(sample_distance)
	}
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

######
# load data
huri <- read.csv("~/Documents/INET-work/references/HuRI_binaryPPI/HuRI_Tong_withSymbol.csv", header = T)
bioplex <- read.delim("/Volumes/GoogleDrive/My Drive/VirHostome_CW/GitHub/data/extended_table/BioPlex.3.0_edge.tsv", header = T)

######
# have network
networks <- list(
	HuRI = simplify(graph_from_data_frame(huri[, c(5:6)], directed = FALSE), remove.loops = TRUE),
	BioPlex = simplify(graph_from_data_frame(bioplex[, c(5:6)], directed = FALSE), remove.loops = TRUE)
	)

######
# have protein list
nodes_huri <- list(
	gwas_2021a = c("CCHCR1", "OAS2", "LZTFL1", "OAS1", "TYK2"),
	gwas_2021a_noParalog = c("CCHCR1", "LZTFL1", "OAS1", "TYK2"),
	gwas_2021b = c("FDX2", "ICAM1", "ICAM3", "ICAM4", "KAT7", "LZTFL1", "NSF", "NXPE3", "OAS1", "OAS2", "PLEKHA4", "RAVER1", "SLC6A20", "SPPL2C", "STH", "TMEM65", "TYK2"),
	gwas_2021b_noParalog = c("FDX2", "ICAM3", "KAT7", "LZTFL1", "NSF", "NXPE3", "OAS1", "PLEKHA4", "RAVER1", "SLC6A20", "SPPL2C", "STH", "TMEM65", "TYK2"),
	ctcl = c("FDX2", "ICAM1", "ICAM3", "ICAM4", "KAT7", "LZTFL1", "OAS1", "OAS2", "RAVER1", "TYK2"),
	ctcl_noParalog = c("FDX2", "ICAM3", "KAT7", "LZTFL1", "OAS1", "RAVER1", "TYK2")
)

nodes_bioplex <- list(
	gwas_2021a = c("LZTFL1", "HLA-G", "CCHCR1", "OAS2", "DPP9", "TYK2", "IFNAR2"),
	gwas_2021b = c("ARHGAP27", "ARL17A", "CEP97", "DPP9", "FOXP4", "ICAM1", "ICAM3", "ICAM4", "ICAM5", "IFNAR2", "LRRC37A2", "LZTFL1", "MAPT", "NUCB1", "OAS2", "PPP1R15A", "RAVER1", "RPL24", "SLC6A20", "TMEM65", "TULP2", "TYK2", "WNT3", "ZBTB11"),
	gwas_2021b_noParalog = c("ARHGAP27", "ARL17A", "CEP97", "DPP9", "FOXP4", "ICAM3", "IFNAR2", "LRRC37A2", "LZTFL1", "MAPT", "NUCB1", "OAS2", "PPP1R15A", "RAVER1", "RPL24", "SLC6A20", "TMEM65", "TULP2", "TYK2", "WNT3", "ZBTB11"),
	ctcl = c("DPP9", "ICAM1", "ICAM3", "ICAM4", "ICAM5", "IFNAR2", "LZTFL1", "OAS2", "RAVER1", "TYK2"),
	ctcl_noParalog = c("DPP9", "ICAM3", "ICAM5", "IFNAR2", "LZTFL1", "OAS2", "RAVER1", "TYK2")
)

######
# randomly select the same number proteins as node list from 1) HuRI and 2) BioPlex, calculating the average shortest path is significant than the observation (as node list) or not

# HuRI
huri_result <- list()
for (no in 1:length(nodes_huri)) {
	right <- "FALSE"
	huri_rand_avg <- c()
	gwas_distance <- distances(networks[["HuRI"]], v = nodes_huri[[no]], to = nodes_huri[[no]])
	gwas_avg_path <- mean(gwas_distance[lower.tri(gwas_distance)])
	huri_rand_avg <- c(huri_rand_avg, mcreplicate(100, avg_path_rand(networks[["HuRI"]], nodes_huri[[no]]), mc.cores = detectCores()))

	sample_notInf <- huri_rand_avg[!is.infinite(huri_rand_avg)]

	sig <- table(sample_notInf >= gwas_avg_path)[right] / length(sample_notInf)
	if (is.na(sig_value)) {
		sig_value <- paste0(" < ", 1/length(sample_notInf))
	} else {
		sig_value <- as.numeric(sig)
	}
	cat(paste0("Node: ", names(nodes_huri[no]), "; Network: ", names(networks["HuRI"]), "; avg shortest path: ", round(gwas_avg_path, 3), "\n"))
	print(paste0("right tail: ", right, "; median = ", round(median(sample_notInf), 3), "; sig = ", round(sig_value, 3)))

}	


# Node: gwas_2021a; Network: huri; avg shortest path: 3.1
# [1] "right tail: FALSE; median = 3.8; sig = 0.044"
# Node: gwas_2021a_noParalog; Network: huri; avg shortest path: 2.833
# [1] "right tail: FALSE; median = 3.833; sig = 0.011"
# Node: gwas_2021b; Network: huri; avg shortest path: 3.515
# [1] "right tail: FALSE; median = 3.897; sig = 0.037"
# Node: gwas_2021b_noParalog; Network: huri; avg shortest path: 3.571
# [1] "right tail: FALSE; median = 3.813; sig = 0.247"
# Node: ctcl; Network: huri; avg shortest path: 3.578
# [1] "right tail: FALSE; median = 3.967; sig = 0.202"
# Node: ctcl_noParalog; Network: huri; avg shortest path: 3.762
# [1] "right tail: FALSE; median = 3.857; sig = 0.388"

# BioPlex3
for (no in 1:length(nodes_bioplex)) {
	right <- "FALSE"
	gwas_distance <- distances(networks[[2]], v = nodes_bioplex[[no]], to = nodes_bioplex[[no]])
	gwas_avg_path <- mean(gwas_distance[lower.tri(gwas_distance)])

	cat(paste0("Node: ", names(nodes_bioplex[no]), "; Network: ", names(networks[2]), "; avg shortest path: ", round(gwas_avg_path, 3), "\n"))
	sig_value <- sig_rand_sample(networks[[2]], nodes_huri[[no]], 100, right, "FALSE")

	print(paste0("right tail: ", right, "; median = ", round(sig_value[1], 3), "; sig = ", round(as.numeric(sig_value[2]), 3)))
}
# Node: gwas_2021a; Network: bioplex; avg shortest path: 3.952
# [1] "right tail: FALSE; median = 3.9; sig = 0.598"
# Node: gwas_2021b; Network: bioplex; avg shortest path: 3.732
# [1] "right tail: FALSE; median = 3.833; sig = 0.47"
# Node: gwas_2021b_noParalog; Network: bioplex; avg shortest path: 3.757
# [1] "right tail: FALSE; median = 3.831; sig = 0.379"
# Node: ctcl; Network: bioplex; avg shortest path: 3.689
# [1] "right tail: FALSE; median = 3.901; sig = 0.224"
# Node: ctcl_noParalog; Network: bioplex; avg shortest path: 3.893
# [1] "right tail: FALSE; median = 3.856; sig = 0.58"

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
