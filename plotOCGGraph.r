# ref: https://github.com/alextkalinka/linkcomm/blob/master/R/OCG_main.R
plotOCGraph <- function(x, clusterids = 1:x$numbers[3], nodes = NULL, pie.local = TRUE, incident = TRUE, layout = layout.fruchterman.reingold, vertex.radius = 0.03, scale.vertices = 0.05, edge.color = "grey", vertex.label.color = "black", vertex.label.cex = 0.5, pal = brewer.pal(7,"Set2"), shownodesin = 0, vlabel = TRUE, random = TRUE, ...)
	# x is an "OCG" object.
	{
	# Make an edgelist based on clusterids or nodes.
	if(is.null(nodes)){
		el <- x$edgelist[getEdgesIn(x, clusterids = clusterids),]
		ig <- graph.edgelist(el,directed = FALSE)
		nodes <- V(ig)$name
	}else{
		if(incident){
			# Only clusters to which named nodes belong.
			clusterids = which.communities(x, nodes = nodes)
			el <- x$edgelist[getEdgesIn(x, clusterids = clusterids),]
			ig <- graph.edgelist(el,directed = FALSE)
			nodes <- V(ig)$name
		}else{
			# Get all clusters to which named nodes and their first-order neighbourhood belong.
			clusterids = which.communities(x, nodes = nodes)
			el <- x$edgelist[getEdgesIn(x, clusterids = clusterids, all = TRUE),]
			nodes <- unique(c(as.character(el[,1]),as.character(el[,2])))
			ig <- graph.edgelist(el,directed = FALSE)
			nodes <- V(ig)$name
			}
		}
	crf <- colorRampPalette(pal,bias=1)
	cols <- crf(length(clusterids))
	if(random){
		cols <- sample(cols,length(clusterids),replace=FALSE)
		}
	names(cols) <- clusterids

	if(shownodesin == 0){
		vnames <- V(ig)$name
	}else{ # Show nodes that belong to more than x number of communities.
		vnames <- V(ig)$name
		inds <- NULL
		for(i in 1:length(vnames)){
			if(x$numclusters[which(names(x$numclusters)==vnames[i])] < shownodesin){
				inds <- append(inds,i)
				}
			}
		vnames[inds] <- ""
		}
	if(vlabel==FALSE){
		vnames = NA
		}

	# Get node community membership by edges.
	if(pie.local){
		edge.memb <- numberEdgesIn(x, clusterids = clusterids, nodes = nodes)
	}else{
		edge.memb <- numberEdgesIn(x, nodes = nodes)
		}

	cat("   Getting node layout...")
	lay <- layout(ig) # vertex x-y positions; will serve as centre points for node pies.
	lay <- layout.norm(lay, xmin=-1, xmax=1, ymin=-1, ymax=1)
	rownames(lay) <- V(ig)$name
	cat("\n")

	node.pies <- .nodePie(edge.memb=edge.memb, layout=lay, nodes=nodes, edges=100, radius=vertex.radius, scale=scale.vertices)
	cat("\n")

	dev.hold(); on.exit(dev.flush())
	# Plot graph.
	plot(ig, layout=lay, vertex.shape="none", vertex.label=NA, vertex.label.dist=0, vertex.size = 2, vertex.label.color=vertex.label.color, edge.color=edge.color, ...)
	labels <- list()
	# Plot node pies and node names.
	for(i in 1:length(node.pies)){
		yp <- NULL
		for(j in 1:length(node.pies[[i]])){
			seg.col <- cols[which(names(cols)==names(edge.memb[[i]])[j])]
			polygon(node.pies[[i]][[j]][,1], node.pies[[i]][[j]][,2], col = seg.col)
			yp <- append(yp, node.pies[[i]][[j]][,2])
			}
		lx <- lay[which(rownames(lay)==names(node.pies[i])),1] + 0.1
		ly <- max(yp) + 0.02 # Highest point of node pie.
		labels[[i]] <- c(lx, ly)
		}
	# Plot node names after nodes so they overlay them.
	# for(i in 1:length(labels)){
	# 	text(labels[[i]][1], labels[[i]][2], labels = vnames[which(nodes==names(node.pies[i]))], cex = vertex.label.cex, col = vertex.label.color)
	# 	}
	}

.nodePie <- function(edge.memb, layout, nodes, edges, radius, scale = NULL)
	{
	# Returns x-y coordinates for node pies using their layout coordinates as centre points.
	base_radius <- radius
	node.pies <- list()

	t2xy <- function(t){
		t2p <- 2*pi*t + 90 * pi/180
		list(x = radius * cos(t2p), y = radius * sin(t2p))
		}
	for(i in 1:length(edge.memb)){
		prog <- paste(c("   Constructing node pies...",floor((i/length(nodes))*100),"%\r"),collapse="")
		cat(prog)
		flush.console()

		x <- edge.memb[[i]]
		x <- c(0, cumsum(x)/sum(x))
		dx <- diff(x)
		if(length(dx)==1){
			if(is.nan(dx)){
				x <- c(0,1)
				dx <- 1
				}
			}

		if(!is.null(scale)){
			if(length(dx) > 1){
				radius <- base_radius + base_radius*scale*length(dx)
			}else{
				radius <- base_radius
				}
			}

		for (j in 1:length(dx)) {
			n <- max(2, floor(edges * dx[j]))
			P <- t2xy(seq.int(x[j], x[j + 1], length.out = n))
			if(n != 100){
				P$x <- c(P$x, 0)
				P$y <- c(P$y, 0)
				}
			P$x <- P$x + layout[i,1]
			P$y <- P$y + layout[i,2]
			mat <- cbind(P$x,P$y)
			colnames(mat) <- c("x","y")
			node.pies[[as.character(nodes[i])]][[j]] <- mat
			}

		}

	return(node.pies)

	}
