# make HuRI community graph
# Lin Chung-wen

#######
# load packages and function
library(linkcomm)
# trace("plot.igraph", edit = TRUE) # set as 0.08
source("~/Documents/INET-work/virus_network/src/plotOCGGraph.r")
######
# load R Data
load("~/Documents/INET-work/virus_network/statistic_results/HuRI_ocg.RData")

######
# plotting
pdf("/tmp/HuRI_403_926.pdf", width = 20, height = 20)
plotOCGraph(huri_ocg, clusterids = c(8, 145, 316, 364, 377, 415, 525, 571, 692, 731, 891, 926, 1353, 1515, 1627, 1652, 1833, 1882, 2398, 2545, 2769, 2831, 3035, 3205, 3564, 3682, 3941, 4160, 4227), vertex.radius = 0.03, vertex.label.cex = 1, layout = layout.sphere)
dev.off()