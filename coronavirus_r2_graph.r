library(igraph)
library(circular)

#################################################
# 1. inet virhost network
node_inet <- read.csv("~/Documents/INET-work/virus_network/Y2H_screening/20200923_verified/node_20200924.csv", header = T, as.is = T)
edge_inet <- read.csv("~/Documents/INET-work/virus_network/Y2H_screening/20200923_verified/edge_20200924.csv", header = T, as.is = T)
net_inet <- graph_from_data_frame(d = edge_inet, vertices = node_inet, directed = F)
V(net_inet)$species <- node_inet$species
V(net_inet)$color <- ifelse(V(net_inet)$species == "human", "blue", "pink")
# coords <- layout.reingold.tilford(net, circular = T)
coords <- layout_(net_inet, with_kk())
node <- read.csv("~/Documents/INET-work/virus_network/Y2H_screening/20200923_verified/node_20200924.csv", header = T, as.is = T)
edge <- read.csv("~/Documents/INET-work/virus_network/Y2H_screening/20200923_verified/edge_20200924.csv", header = T, as.is = T)
net <- graph_from_data_frame(d = edge, vertices = node, directed = F)
V(net)$species <- node$species
V(net)$color <- ifelse(V(net)$species == "human", "blue", "pink")
# coords <- layout.reingold.tilford(net, circular = T)
coords <- layout_(net, with_kk())
coords <- layout_(net, with_lgl())
# radian.rescale <- function(x, start = 0, direction = 1) {
#   c.rotate <- function(x) (x + start) %% (2 * pi) * direction
#   c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
# }
# lab.locs <- radian.rescale(x = 1:dim(coords)[1], direction = -1, start = 0)
pdf("inet2.pdf", width = 15, height = 15)
plot(net_inet, layout = coords, vertex.size = 3, vertex.label.cex = .9, vertex.label.dist = 0.5, vertex.label.degree = 1, vertex.label.color = "black", margin = 0.1)
dev.off()

# 2. Roth lab fBFG virhost network
node_dk <- read.csv("~/Documents/INET-work/virus_network/references/Roth_fBFG/node.csv", header = T, as.is = T)
edge_dk <- read.csv("~/Documents/INET-work/virus_network/references/Roth_fBFG/edge.csv", header = T, as.is = T)
net_dk <- graph_from_data_frame(d = edge_dk, vertices = node_dk, directed = F)
V(net_dk)$species <- node_dk$species
V(net_dk)$color <- ifelse(V(net_dk)$species == "human", "darkgreen", "pink")
coords <- layout_(net_dk, with_kk())
pdf("Roth.pdf", width = 15, height = 15)
plot(net_dk, layout = coords, vertex.size = 3, vertex.label.cex = .9, vertex.label.dist = 0.5, vertex.label.degree = 1, vertex.label.color = "black", margin = 0.1)
dev.off()

## intersection of net_inet and net_dk
net_bi <- graph.union(net_inet, net_dk, byname = T)
V(net_bi)$species <- c(rep("human", 124), rep("virus", 154 - 124), rep("human", 376 - 154))
V(net_bi)$color <- ifelse(V(net_bi)$species == "human", "blue", "pink")
coords <- layout_(net_bi, with_kk())
pdf("INETandRoth_intersection.pdf", width = 15, height = 15)
plot(net_bi, layout = coords, vertex.size = 3, vertex.label.cex = .9, vertex.label.dist = 0.5, vertex.label.degree = 1, vertex.label.color = "black", margin = 0.1)
dev.off()


=======
plot(net_inet, layout = coords, vertex.size = 3, vertex.label.cex = .9, vertex.label.dist = 0.5, vertex.label.degree = 1, vertex.label.color = "black", margin = 0.1)
dev.off()

# 2. Roth lab fBFG virhost network
node_dk <- read.csv("~/Documents/INET-work/virus_network/references/Roth_fBFG/node.csv", header = T, as.is = T)
edge_dk <- read.csv("~/Documents/INET-work/virus_network/references/Roth_fBFG/edge.csv", header = T, as.is = T)
net_dk <- graph_from_data_frame(d = edge_dk, vertices = node_dk, directed = F)
V(net_dk)$species <- node_dk$species
V(net_dk)$color <- ifelse(V(net_dk)$species == "human", "darkgreen", "pink")
coords <- layout_(net_dk, with_kk())
pdf("Roth.pdf", width = 15, height = 15)
plot(net_dk, layout = coords, vertex.size = 3, vertex.label.cex = .9, vertex.label.dist = 0.5, vertex.label.degree = 1, vertex.label.color = "black", margin = 0.1)
dev.off()

## intersection of net_inet and net_dk
net_bi <- graph.union(net_inet, net_dk, byname = T)
V(net_bi)$species <- c(rep("human", 124), rep("virus", 154 - 124), rep("human", 376 - 154))
V(net_bi)$color <- ifelse(V(net_bi)$species == "human", "blue", "pink")
coords <- layout_(net_bi, with_kk())
pdf("INETandRoth_intersection.pdf", width = 15, height = 15)
plot(net_bi, layout = coords, vertex.size = 3, vertex.label.cex = .9, vertex.label.dist = 0.5, vertex.label.degree = 1, vertex.label.color = "black", margin = 0.1)
dev.off()


plot(net, layout = coords, vertex.size = 3, vertex.label.cex = .9, vertex.label.dist = 0.5, vertex.label.degree = 1, vertex.label.color = "black", margin = 0.1)
dev.off()

layouts <- grep("^layout_", ls("package:igraph"), value=TRUE)[-1]
# Remove layouts that do not apply to our graph.
layouts <- layouts[!grepl("bipartite|merge|norm|sugiyama|tree", layouts)]

par(mfrow=c(3,3), mar=c(1,1,1,1))
pdf("inet2.pdf", width = 15, height = 15)
for (layout in layouts) {
    print(layout)
    l <- do.call(layout, list(net))
    plot(net, layout = l, vertex.size = 3, vertex.label.cex = .9, vertex.label.dist = 0.5, vertex.label.degree = 1, vertex.label.color = "black", margin = 0.1, main=layout)
}
dev.off()
#################################################################
# change vertex label
# ref: https://stackoverflow.com/questions/23209802/placing-vertex-label-outside-a-circular-layout-in-igraph

### Here's one way to do it.

# library(igraph)
# library(ggplot2)
# library(scales)

# ## The igraph docs say that vertex.label.degree controls the position
# ## of the labels with respect to the vertices. It's interpreted as a
# ## radian, like this:
# ##
# ## Value is : Label appears ... the node
# ## -pi/2: above
# ## 0: to the right of
# ## pi/2: below
# ## pi: to the left of
# ##
# ## We can generalize this. vertex.label.degree can take a vector as
# ## well as a scalar for its argument. So we write a function to
# ## calculate the right position for a label based on its vertex's location
# ## on the circle.

# ## Get the labels aligned consistently around the edge of the circle
# ## for any n of nodes.
# ## This code borrows bits of ggplot2's polar_coord function
# ## start = offset from 12 o'clock in radians
# ## direction = 1 for clockwise; -1 for anti-clockwise.

# radian.rescale <- function(x, start=0, direction=1) {
#   c.rotate <- function(x) (x + start) %% (2 * pi) * direction
#   c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
# }

# ### Example
# ## Generate some fake data
# n <- 15
# g <- erdos.renyi.game(n, 0.5)
# ## Obviously labeling in this way this only makes sense for graphs
# ## laid out as a circle to begin with
# la <- layout.circle(g)

# lab.locs <- radian.rescale(x=1:n, direction=-1, start=0)
# plot(g, layout=la, vertex.size=2, vertex.label.dist=1,
#      vertex.label.degree=lab.locs)
