## Function to prune tree using 'richness' information, assumed to have minimally two columns, "taxon" and "n.taxa"
##   Perhaps relax on these column names, may cause too many problems
## May also include 'exemplar' column; in that case, rename relevant tip.label before pruning.
prepareData <- function (phy, richness, verbose)
{
	if (is.null(richness)) { # Assume tree represents single species tips and is completely sampled
		richness <- data.frame(taxon=phy$tip.label, n.taxa=1);
	}
# Rename exemplar taxa with taxon name in richness file
	if (!is.null(richness$exemplar))
	{
# Change relevant tip.labels in phy; individual 'exemplar' may be NA, use original tip.label.
# Ordering in richness file should NOT be assumed to match order of tip.labels
		i.na <- is.na(richness$exemplar);
		phy$tip.label[match(richness$exemplar[!i.na], phy$tip.label)] <- as.character(richness$taxon[!i.na]);
	}
	
# make sure things are in the correct order and of correct format
	if (length(richness[1,]) == 2)
	{
		if (colnames(richness)[1] != "taxon" || colnames(richness)[2] != "n.taxa")
		{
			if (class(richness[,1]) == "factor" & class(richness[,2]) == "integer")
			{
				colnames(richness) = c("taxon", "n.taxa");
			} else if (class(richness[,1]) == "integer" & class(richness[,2]) == "factor")
			{
				colnames(richness) = c("n.taxa", "taxon");
			} else {
				cat("turboMEDUSA thinks your richness data is in an incorrect format. See ?runTurboMEDUSA.\n")
				stop;
			}
		}
	}
	
# Prune tree down to lineages with assigned richnesses
	temp <- richness[, "n.taxa"];
	names(temp) <- richness[, "taxon"];
	pruned <- treedata(phy, temp, warnings=verbose)  # geiger function calling ape (namecheck)
	
# checking for typos
	if (length(pruned$phy$tip.label) == 0) {
		cat("\n\nWARNING: turboMEDUSA encountered a serious error. Tree has no tips after processing richness information! \nIt is likely that an incorrect richness file is being used. \nAnalysis cannot proceed. Please examine the information below to identify the error.\\nn");
		cat("\nTree tip labels:\n");
		print(phy$tip.label);
		cat("\nRichness taxon labels:\n");
		print(as.character(richness[, "taxon"]));
		stop;
	} else if (length(phy$tip.label) != length(pruned$phy$tip.label)) {
		cat("turboMEDUSA thinks there is a typo in either the tree or richness files, as one or more tips were dropped from the tree.\n");
		stop;
	}

	phy <- pruned$phy;
# Check the tree
	#	plotNN(phy)					# Node numbers (ape-style) plotted
	
	return(list(phy=phy, richness=richness));
}


## Original default was to fit 20 models (or less if the tree was small).
## Changing to a stop-criterion (stop="modelLimit") e.g. when k = n-1 (i.e. when denominator of aicc correction is undefined).
## k <- (3*i-1) # when both birth and death are estimated, where i is the number of piecewise models
  ## This occurs when i = n/3
  ## If Yule, max i = n/2
## n <- (2*num.taxa - 1) == (2*length(richness[,1]) - 1) # i.e. total number of nodes in tree (internal + pendant)
## Alternatively use aicc threshold itself as a stopping criterion (stop="threshold").
# AICc = AIC + 2*k*(k+1)/(n-k-1);
getMaxModelLimit <- function (richness, modelLimit, model, stop, verbose)
{
	samp.size <- (2*length(richness[,1]) - 1)
	if (model == "bd" || model == "mixed")
	{
		max.modelLimit <- as.integer(samp.size/3) - ((!(samp.size %% 3)) * 1);
	} else {
		max.modelLimit <- as.integer(samp.size/2) - ((!(samp.size %% 2)) * 1);
	}
	
	if (stop == "modelLimit")
	{
		if (modelLimit > max.modelLimit) {modelLimit <- max.modelLimit;}
	} else {
		modelLimit <- max.modelLimit;
	}
	
	if (verbose)
	{
		cat("\nLimiting consideration to a maximum of ", modelLimit, " piecewise", sep="");
		if (model == "bd")
		{
			cat(" birth-death models")
		} else if (model == "mixed")
		{
			cat(" mixed models")
		} else {
			cat(" pure-birth (Yule) models");
		}
		if (stop == "threshold") {cat(" (or until threshold is not satisfied)");}
		cat(".\n\n")
	}
	return(modelLimit);
}


## Fitted curve from random b-d simulations
## Value corresponds to 95th percentile of AICc(split) - AICc(no-split) for no-split simulations
## x-shifted power function
getThreshold <- function (x)
{
	a = -3.5941052380332650E+01;
	b =  6.7372587299747000E+00;
	c = -1.0061508340754866E-01;
	Offset =  2.7516678664333408E+01;
	y <- a * (x-b)^c + Offset;
	if (y < 0 || is.nan(y)) y <- 0;
	return(y);
}


## The makeCacheMedusa function is like the first half of the original splitEdgeMatrix().
## It works through and reorders the edges, then works out start and end times of these
## based on the phylogeny's branching times.
##
## In addition, every node's descendants are also calculated.  The element 'desc' is a list.
## $desc[i] contains the indices within $edge, $t.start, etc., of all descendants of node 'i'
## (in ape node numbering format).
makeCacheMedusa <- function (phy, richness, all.nodes, mc, numCores)
{
	n.tips <- length(phy$tip.label);
	n.int <- nrow(phy$edge) - n.tips;
	
## Ape numbers the tips first
	i.int <- seq_len(n.int);
	interior <- phy$edge[,2] %in% phy$edge[,1];
	bt <- branching.times(phy);
	
# Consider only internal edges first
	edges.int <- phy$edge[interior,];
	colnames(edges.int) <- c("anc", "dec");
	
	t.0 <- bt[match(edges.int[,1], (n.tips+1):max(edges.int))];
	t.1 <- c(t.0[i.int] - phy$edge.length[interior]);
	
	z.internal <- cbind(edges.int, t.0, t.1, t.len=t.0 - t.1,
		n.0=rep(1, n.int), n.t=rep(NA, n.int));
	
# Now, pendant edges; 
	edges.pendant <- phy$edge[match(seq_len(n.tips), phy$edge[,2]),];
	colnames(edges.pendant) <- c("anc", "dec");
	
	t.0 <- bt[match(edges.pendant[,1], (n.tips+1):max(edges.pendant))];
	t.1 <- rep(0, n.tips);
# cannot assume richness ordering necessarily matches that of tip labels
	ext.richness <- richness$n.taxa[match(phy$tip.label, richness$taxon)];
	
	z.pendant <- cbind(edges.pendant, t.0, t.1, t.len=t.0 - t.1,
		n.0=rep(1, n.tips), n.t=ext.richness);
	
	z <- rbind(z.internal, z.pendant);
	z <- cbind(z,partition=rep(1, length(z[,1]))); # Stores piecewise model structure
	rownames(z) <- NULL;
	
# Used for identifying descendant nodes below i.e. tracking breakpoints
	all.edges <- as.matrix(z[,c("anc","dec")]);
	desc.stem <- list();
	desc.node <- list();
	
	if (mc)
	{
		desc.stem <- mclapply(seq_len(max(all.edges)), descendantsCutAtStem.idx, all.edges=all.edges, mc.cores=numCores);
		desc.node <- mclapply(seq_len(max(all.edges)), descendantsCutAtNode.idx, all.edges=all.edges, mc.cores=numCores);
	} else {
		desc.stem <- lapply(seq_len(max(all.edges)), descendantsCutAtStem.idx, all.edges=all.edges);
		desc.node <- lapply(seq_len(max(all.edges)), descendantsCutAtNode.idx, all.edges=all.edges);
	}
	
	# for (i in 1:length(desc.node)) #  *** NEED TO FIX THIS SHIT!!! ***
	# {
		# if (length(desc.node[[i]]) == 0) # tips
		# {
			# desc.node[[i]] = descendantsCutAtStem.idx(node.list=i, all.edges=all.edges)
		# }
	# }
	
## Needed downstream; don't recalculate
 ## Gives the number of tips associated with an internal node; determines whether a node is 'virgin' or not
	num.tips <- list()
	if (mc)
	{
		num.tips <- mclapply(all.nodes, getNumTips, phy=phy, mc.cores=numCores);
	} else {
		num.tips <- lapply(all.nodes, getNumTips, phy=phy);
	}
	
	res <- list(z=z, desc.stem=desc.stem, desc.node=desc.node, num.tips=num.tips);
	return(res);
}


## Needed for determining whther nodes are virgin nodes
getNumTips <- function (node, phy)
{
	n <- length(node.leaves(phy,node));
	return(n);
}


## This generates the indices of all descendants of a node, using ape's edge matrix.
## Deals with row numbers of the edge matrix rather than node numbers of the tree.
descendantsCutAtStem <- function (node, all.edges)
{
	ans <- numeric();
	ans <- node;
	repeat {
		node <- all.edges[all.edges[,1] %in% node,2];
		if (length(node) > 0)
		{
			ans <- c(ans, node);
		} else {break;}
	}
	return(unlist(ans));
}


## The function 'descendants' returns the indices of all descendants within the edge matrix.
descendantsCutAtStem.idx <- function (node.list, all.edges)
{
	which(all.edges[,1] == node.list | all.edges[,2] %in% descendantsCutAtStem(node.list, all.edges));
}


## This generates the indices of all descendants of a node, using ape's edge matrix.
## Deals with row numbers of the edge matrix rather than node numbers of the tree.
descendantsCutAtNode <- function (node, all.edges)
{
	ans <- numeric();
	repeat {
		node <- all.edges[all.edges[,1] %in% node,2];
		if (length(node) > 0)
		{
			ans <- c(ans, node);
		} else {break;}
	}
	return(unlist(ans));
}


## The function 'descendants' returns the indices of all descendants within the edge matrix.
descendantsCutAtNode.idx <- function (node.list, all.edges)
{
	which(all.edges[,1] == node.list | all.edges[,2] %in% descendantsCutAtNode(node.list, all.edges));
}


#Identify user-desired node. May involve pruned data, and/or mrca
getNodeNumber <- function (phy, richness, mrca, verbose)
{
	nodeNumber <- integer();
	if (!is.null(richness))
	{
		phy <- prepareData(phy=phy, richness=richness, verbose=verbose)$phy;
	}
	if (!is.null(mrca))
	{
		
		
		
		
	}
	
	
	return(nodeNumber);
}


findMrca = function(phy, tips)
{
	tt <- match(tips, phy$tip.label);
	getMRCA(phy, tt); #ape function
}

extractCladeFromTree <- function (phy, tips)
{
	node <- findMrca(phy, tips)
	cat(node, "\n")
	phy <- extract.clade(phy, node); # ape function

	return(phy)
}
