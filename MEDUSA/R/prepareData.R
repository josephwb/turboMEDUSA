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
				cat("MEDUSA thinks your richness data is in an incorrect format. See ?MEDUSA.\n")
				stop;
			}
		}
	}
	
	if (class(phy) != "phylo") {cat("\n\nWARNING: tree is not of class \"phylo\". Stopping.\n"); stop;}
	
# Prune tree down to lineages with assigned richnesses
	temp <- richness[, "n.taxa"];
	names(temp) <- richness[, "taxon"];
	pruned <- treedata(phy, temp, warnings=verbose)  # geiger function calling ape (namecheck)
	
# checking for typos
	if (length(pruned$phy$tip.label) == 0) {
		cat("\n\nWARNING: MEDUSA encountered a serious error. Tree has no tips after processing richness information! \nIt is likely that an incorrect richness file is being used. \nAnalysis cannot proceed. Please examine the information below to identify the error.\\nn");
		cat("\nTree tip labels:\n");
		print(phy$tip.label);
		cat("\nRichness taxon labels:\n");
		print(as.character(richness[, "taxon"]));
		stop;
	} else if (length(phy$tip.label) != length(pruned$phy$tip.label)) {
		cat("MEDUSA thinks there is a typo in either the tree or richness files, as one or more tips were dropped from the tree.\n");
		stop;
	}

	phy <- pruned$phy;
# Check the tree
	#	plotNN(phy)					# Node numbers (ape-style) plotted
	
	return(list(phy=phy, richness=richness));
}

# Determine desired model from passed in fixed parameters (if present)
configureModel <- function (model, epsilon, r, b, d, initialR, initialE)
{
	sp <- NULL;
	fixPar <- NULL;
	if (!is.null(epsilon)) # user-defined epsilon
	{
		if (epsilon <= 0 | epsilon >= 1) {cat("\n\nWARNING: value of epsilon (", epsilon, ") is invalid; must be > 0 and < 1. Stopping analysis.\n", sep=""); stop;}
		sp <- c(initialR, epsilon);
		fixPar <- epsilon;
		model <- "fixedEpsilon";
	} else if (!is.null(r)) # user-defined net diversification rate
	{
		if (r <= 0) {cat("\n\nWARNING: value of r (", r, ") is invalid; must be > 0. Stopping analysis.\n", sep=""); stop;}
		sp <- c(r, initialE);
		fixPar <- r;
		model <- "fixedR";
	} else if (!is.null(d)) # user-defined extiction rate
	{
		if (d <= 0) {cat("\n\nWARNING: value of d (", d, ") is invalid; must be > 0. Stopping analysis.\n", sep=""); stop;}
		sp <- c(initialR, d);
		fixPar <- d;
		model <- "fixedD";
	} else if (!is.null(b)) # user-defined speciation rate
	{
		if (b <= 0) {cat("\n\nWARNING: value of b (", b, ") is invalid; must be > 0. Stopping analysis.\n", sep=""); stop;}
		sp <- c((b/2), initialE);
		fixPar <- b;
		model <- "fixedB";
	} else {
		sp <- c(initialR, initialE);
	}
	return(list(model=model, sp=sp, fixPar=fixPar))
}


## Original default was to fit 20 models (or less if the tree was small).
## Changing to a stop-criterion (stop="modelLimit") e.g. when k = n-1 (i.e. when denominator of aicc correction is undefined).
## k <- (3*i-1) # when both birth and death are estimated, where i is the number of piecewise models
  ## This occurs when i = n/3
  ## If Yule, max i = n/2
## n <- (2*num.taxa - 1) == (2*length(richness[,1]) - 1) # i.e. total number of nodes in tree (internal + pendant)
## Alternatively use aicc threshold itself as a stopping criterion (stop="threshold").
# AICc = AIC + 2*k*(k+1)/(n-k-1);
getMaxModelLimit <- function (richness, modelLimit, model, stop)
{
	samp.size <- (2*length(richness[,1]) - 1)
	if (model == "bd" || model == "mixed")
	{
		max.modelLimit <- as.integer(samp.size/3) - ((!(samp.size %% 3)) * 1);
	} else {
		max.modelLimit <- as.integer(samp.size/2) - ((!(samp.size %% 2)) * 1); # models estimating only one diversification parameter
	}
	
	if (stop == "modelLimit")
	{
		if (modelLimit > max.modelLimit) {modelLimit <- max.modelLimit;}
		cat("Limiting consideration to ", modelLimit, " piecewise", sep="");
		if (model == "bd") {
			cat(" birth-death models")
		} else if (model == "mixed") {
			cat(" mixed models")
		} else if (model == "yule") {
			cat(" pure-birth (Yule) models");
		} else {
			cat(" diversification models")
		}
		cat(".\n\n");
	}
	
	return(modelLimit);
}


## Fitted curve from random b-d simulations
## Value corresponds to 95th percentile of AICc(split) - AICc(no-split) for no-split simulations
## x-shifted power function
getThreshold <- function (treeSize, fixThreshold, stop)
{
	a = -3.5941052380332650E+01;
	b =  6.7372587299747000E+00;
	c = -1.0061508340754866E-01;
	Offset =  2.7516678664333408E+01;
	threshold <- a * (treeSize-b)^c + Offset;
	if (threshold < 0 || is.nan(threshold)) threshold <- 0;
	
	if (!is.null(fixThreshold))
	{
		if (is.nan(fixThreshold) || class(fixThreshold) != "numeric" || fixThreshold < 0)
		{
			cat("Provided threshold value of '", fixThreshold, "' is invalid.\n", sep="")
			cat("Will proceed with value determined from simulations.\n")
			cat("Appropriate AICc threshold for tree of ", treeSize, " tips is: ", threshold, ".\n", sep="");
			return(threshold);
		} else {
			cat("Using provided threshold value of ", fixThreshold, ".\n", sep="")
			cat("From simulations, appropriate AICc threshold for tree of ", treeSize, " tips would be: ", threshold, ".\n", sep="");
			return(fixThreshold);
		}
	} else {
		if (stop == "threshold") {
			cat("Using AIC-threshold as analysis-terminating criterion.\n");
			cat("Appropriate threshold for tree of ", treeSize, " tips is: ", threshold, ".\n\n", sep="");
		}
		return(threshold);
	}
}


## The makeCacheMedusa function is like the first half of the original splitEdgeMatrix().
## It works through and reorders the edges, then works out start and end times of these
## based on the phylogeny's branching times.
##
## In addition, every node's descendants are also calculated.  The element 'desc' is a list.
## $desc[i] contains the indices within $edge, $t.start, etc., of all descendants of node 'i'
## (in ape node numbering format).
makeCacheMedusa <- function (phy, richness, all.nodes, shiftCut, mc, numCores, verbose=TRUE)
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
	
	if (verbose) cat("  Gathering descendant node information...");
	if (mc)
	{
		if (shiftCut == "both" || shiftCut == "stem")
		{
			desc.stem <- mclapply(seq_len(max(all.edges)), descendantsCutAtStem.idx, all.edges=all.edges, mc.cores=numCores);
		}
		if (shiftCut == "both" || shiftCut == "node")
		{
			if (!is.null(desc.stem))
			{
				root <- min(z[,"anc"]);
				desc.node <- mclapply(desc.stem, stripStem, mc.cores=numCores);
				desc.node[root] <- desc.stem[root];
			} else {
				desc.node <- mclapply(seq_len(max(all.edges)), descendantsCutAtNode.idx, all.edges=all.edges, mc.cores=numCores);
			}
		}
	} else {
		if (shiftCut == "both" || shiftCut == "stem")
		{
			desc.stem <- lapply(seq_len(max(all.edges)), descendantsCutAtStem.idx, all.edges=all.edges);
		}
		if (shiftCut == "both" || shiftCut == "node")
		{
			if (!is.null(desc.stem))
			{
				root <- min(z[,"anc"]);
				desc.node <- lapply(desc.stem, stripStem);
				desc.node[root] <- desc.stem[root];
			} else {
				desc.node <- lapply(seq_len(max(all.edges)), descendantsCutAtNode.idx, all.edges=all.edges);
			}
		}
	}
	if (verbose) cat(" done.\n");
	
## Needed downstream; don't recalculate
 ## Gives the number of tips associated with an internal node; determines whether a node is 'virgin' or not
	num.tips <- list()
	if (verbose) cat("  Gathering tip richness information...");
	if (mc)
	{
		num.tips <- mclapply(all.nodes, getNumTips, phy=phy, totalTips=n.tips, mc.cores=numCores);
	} else {
		num.tips <- lapply(all.nodes, getNumTips, phy=phy, totalTips=n.tips);
	}
	if (verbose) cat(" done.\n");
	
	res <- list(z=z, desc.stem=desc.stem, desc.node=desc.node, num.tips=num.tips);
	return(res);
}


## Get the number of tips descended from internal node 'node'.
## Needed for determining whether nodes are virgin nodes.
## Uses code from geiger functions 'node.sons' and 'node.leaves'.
getNumTips <- function (node, phy, totalTips=NULL)
{
	if (is.null(totalTips)) totalTips <- length(phy$tip.label);
	if (node <= totalTips) return(1);
	
	n <- 0;
	d <- phy$edge[which(phy$edge[, 1] == node), 2];
	for (j in d) {
		if (j <= totalTips)
			n <- n + 1
		else n <- n + getNumTips(j, phy, totalTips)
	}
	return(n);
}





# *** MAYBE RECODE descendantsCutAtStem IN C++ AS IT IS SLOW FOR LARGE TREES ***

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
	return(ans);
}


## The function 'descendants' returns the indices of all descendants within the edge matrix.
descendantsCutAtStem.idx <- function (node.list, all.edges)
{
	which(all.edges[,1] == node.list | all.edges[,2] %in% descendantsCutAtStem(node.list, all.edges));
}


# Remove stem node from previously calculated set of descendants; about a billion times faster than descendantsCutAtNode
stripStem <- function (x)
{
	y <- unlist(x);
	return(y[-1]);
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
	return(ans);
}


## The function 'descendants' returns the indices of all descendants within the edge matrix.
descendantsCutAtNode.idx <- function (node.list, all.edges)
{
	which(all.edges[,1] == node.list | all.edges[,2] %in% descendantsCutAtNode(node.list, all.edges));
}


## Check that provided arguments are valid
checkValidArguments <- function (phy, richness, model, modelLimit, stop, shiftCut, criterion, stepBack,
	preserveModelFlavour, epsilon, r, b, d, fixThreshold, initialR, initialE,
	verbose, mc, numCores)
{
	if (class(phy) != "phylo" && class(phy) != "multiPhylo") {cat("\n\nWARNING: tree is not of class \"phylo\". Stopping.\n"); stop;}
	
## String arguments
	model=match.arg(model, choices=c("mixed", "bd", "yule"));
	stop=match.arg(stop, choices=c("threshold","modelLimit"));
	shiftCut=match.arg(shiftCut, choices=c("both", "stem", "node"));
	criterion=match.arg(criterion, choices=c("aicc", "aic"));
	
## Boolean arguments
	if (class(stepBack) != "logical") {cat("\n\nWARNING: argument \"stepBack\"is not of class \"logical\". Stopping.\n"); stop;}
	if (class(preserveModelFlavour) != "logical") {cat("\n\nWARNING: argument \"preserveModelFlavour\"is not of class \"logical\". Stopping.\n"); stop;}
	if (class(verbose) != "logical") {cat("\n\nWARNING: argument \"verbose\"is not of class \"logical\". Stopping.\n"); stop;}
	if (class(mc) != "logical") {cat("\n\nWARNING: argument \"mc\"is not of class \"logical\". Stopping.\n"); stop;}
	
## Numeric arguments
	if (class(modelLimit) != "numeric" && class(modelLimit) != "NULL") {cat("\n\nWARNING: argument \"modelLimit\"is not valid. Expecting 'NULL' or numeric. Stopping.\n"); stop;}
	if (class(epsilon) != "numeric" && class(epsilon) != "NULL") {cat("\n\nWARNING: argument \"epsilon\"is not valid. Expecting 'NULL' or numeric. Stopping.\n"); stop;}
	if (class(r) != "numeric" && class(r) != "NULL") {cat("\n\nWARNING: argument \"r\"is not valid. Expecting 'NULL' or numeric. Stopping.\n"); stop;}
	if (class(b) != "numeric" && class(b) != "NULL") {cat("\n\nWARNING: argument \"b\"is not valid. Expecting 'NULL' or numeric. Stopping.\n"); stop;}
	if (class(d) != "numeric" && class(d) != "NULL") {cat("\n\nWARNING: argument \"d\"is not valid. Expecting 'NULL' or numeric. Stopping.\n"); stop;}
	if (class(fixThreshold) != "numeric" && class(fixThreshold) != "NULL") {cat("\n\nWARNING: argument \"fixThreshold\"is not valid. Expecting 'NULL' or numeric. Stopping.\n"); stop;}
	if (class(initialR) != "numeric" && class(initialR) != "NULL") {cat("\n\nWARNING: argument \"initialR\"is not valid. Expecting 'NULL' or numeric. Stopping.\n"); stop;}
	if (class(initialE) != "numeric" && class(initialE) != "NULL") {cat("\n\nWARNING: argument \"initialE\"is not valid. Expecting 'NULL' or numeric. Stopping.\n"); stop;}
	if (class(numCores) != "numeric" && class(numCores) != "NULL") {cat("\n\nWARNING: argument \"numCores\"is invalid. Expecting 'NULL' or numeric. Stopping.\n"); stop;}
}