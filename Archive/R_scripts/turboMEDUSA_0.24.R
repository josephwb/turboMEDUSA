require(geiger)

# You wil want to source this file:
# source("turboMEDUSA_0.23.R")

# Three commands will be of interest:

# 1. Run it
# results <- runTurboMEDUSA (phy, richness=NULL, model.limit=20, stop="model.limit", model="bd",
#	criterion="aicc", initial.r=0.05, initial.e=0.5, plotFig=FALSE, nexus=FALSE,
#	verbose=TRUE, ...)

# 2. Summarize results
# treeParameters <- summarizeTurboMEDUSA (results, model=NULL, cutoff="threshold", criterion="aicc", plotTree=TRUE,
#	time=TRUE, node.labels=TRUE, cex=0.5, plotSurface=FALSE, n.points=100, ...)

# - plotting parameters for optimal tree are stored in treeParameters

# 3. Plot pretty tree
# plotPrettyTree (treeParameters, time=TRUE, node.labels=FALSE, cex=0.5, ...)

# Please email me if you have any questions: josephwb@uidaho.edu

runTurboMEDUSA <-
function(phy, richness=NULL, model.limit=20, stop="model.limit", model="bd",
	criterion="aicc", initial.r=0.05, initial.e=0.5, plotFig=FALSE, nexus=FALSE,
	verbose=TRUE, ...)
{
	if (nexus) phy <- read.nexus(phy)
	if (is.null(richness))  # Assume tree represents single species tips and is completely sampled
	{
		richness <- data.frame(taxon=phy$tip.label, n.taxa=1)
	} else {
## Before determining model.limit, prune tree as necessary (from 'taxon' information in 'richness')
		phyData <- prune.tree.merge.data(phy, richness, verbose)
		phy <- phyData$phy
		richness <- phyData$richness
	}
	
## Limit on number of piecewise models fitted; based on tree size, aicc correction factor, 
## and flavour of model fitted (i.e. # parameters estimated; birth-death or pure-birth)
	model.limit <- get.max.model.limit(richness, model.limit, model, stop, verbose)
	
## Determine correct AICc threshold from tree size (based on simulations)
 ## Should be used for interpreting model-fit
	threshold <- get.threshold(length(phy$tip.label))
	cat("Appropriate AICc threshold for tree of ", length(phy$tip.label), " tips is: ", threshold, ".\n\n", sep="")
	
## Store pertinent information: branch times, richness, ancestors
	cat("Preparing data for analysis... ")
	obj <- make.cache.medusa(phy, richness)
	cat("done.\n")
	
## Keep track of all nodes, internal and pendant (for keeping track of breakpoints)
	pend.nodes <- seq_len(length(phy$tip.label))   # Calculate pendant splits just once, keep track through various models
	int.nodes <- (length(phy$tip.label)+2):max(phy$edge) # Omit root node
	root.node <- length(phy$tip.label)+1
	all.nodes <- c(pend.nodes, root.node, int.nodes)
	
	anc <- obj$anc
	z <- obj$z
	z.orig <- z # Save for summarizing models later
	
## Pre-fit pendant edges so these values need not be re(re(re))calculated; amounts to ~25% of all calculations
 ## Will show particular performance gain for edges with many fossil observations
	cat("Optimizing parameters for pendant edges... ")
	tips.bd <- list(); tips.yule <- list();
	
	if (model == "bd" | model == "mixed")
	{
		tips.bd <- lapply(pend.nodes, medusa.ml.prefit, z, anc, initial.r, initial.e, model="bd")
	}
	if (model == "yule" | model == "mixed")
	{
		tips.yule <- lapply(pend.nodes, medusa.ml.prefit, z, anc, initial.r, initial.e, model="yule")
	}
	tips <- list(bd=tips.bd, yule=tips.yule)
	cat("done.\n")
	
## Pre-fit virgin internal nodes; should deliver performance gain for early models, and especially for large trees
 ## Remain useful until a spilt is accepted within the clade
	cat("Pre-calculating parameters for virgin internal nodes... ")
	virgin.nodes.bd <- list(); virgin.nodes.yule <- list();
	
	if (model == "bd" | model == "mixed")
	{
		virgin.nodes.bd <- lapply(int.nodes, medusa.ml.prefit, z, anc, initial.r, initial.e, model="bd")
	}
	if (model == "yule" | model == "mixed")
	{
		virgin.nodes.yule <- lapply(int.nodes, medusa.ml.prefit, z, anc, initial.r, initial.e, model="yule")
	}
	virgin.nodes <- list(bd=virgin.nodes.bd, yule=virgin.nodes.yule)
	cat("done.\n\n")
	
	prefit <- list(tips=tips, virgin.nodes=virgin.nodes)
	
## Needed downstream; don't recalculate
 ## Gives the number of tips associated with an internal node; determines whether a node is 'virgin' or not
	num.tips <- list()
	num.tips <- lapply(all.nodes, get.num.tips, phy)
	
## 'fit' holds current results; useful for initializing subsequent models
	if (model == "mixed")
	{
		fit.bd <- medusa.ml.initial(z, initial.r, initial.e, model="bd")
		fit.yule <- medusa.ml.initial(z, initial.r, initial.e, model="yule")
		if (fit.bd[[criterion]] < fit.yule[[criterion]]) {
			fit <- fit.bd
		} else {
			fit <- fit.yule
		}
	} else {
		fit <- medusa.ml.initial(z, initial.r, initial.e, model)
	}
	
	models <- list(fit)
	
	if (stop == "model.limit")
	{
		cat("Step 1 (of ", model.limit, "): best likelihood = ", models[[1]]$lnLik, "; AICc = ", models[[1]]$aicc, "\n", sep="")
		for (i in seq_len(model.limit-1))
		{
			node.list <- all.nodes[-fit$split.at]
			res <- lapply(node.list, medusa.ml.update, z, anc, fit, prefit, num.tips, root.node, model, criterion)
			
# Select model with best score according to the specific criterion employed (default aicc)
			best <- which.min(unlist(lapply(res, "[[", criterion)))
			models <- c(models, res[best])
			z <- medusa.split(node.list[best], z, anc)$z
			fit <- res[[best]]   # keep track of '$split.at' i.e. nodes already considered
			
			cat("Step ", i+1, " (of ", model.limit, "): best likelihood = ", models[[i+1]]$lnLik, "; AICc = ", models[[i+1]]$aicc,
				"; break at node ", models[[i+1]]$split.at[i+1],"\n", sep="")
		}
	} else if (stop == "threshold") {
		i <- 1
		done <- FALSE
		cat("Step 1: best likelihood = ", models[[1]]$lnLik, "; AICc = ", models[[1]]$aicc, "\n", sep="")
		while (!done & i < model.limit)
		{
			node.list <- all.nodes[-fit$split.at]
			res <- lapply(node.list, medusa.ml.update, z, anc, fit, prefit, num.tips, root.node, model, criterion)
			
# Select model with best score according to the specific criterion employed (default aicc)
			best <- which.min(unlist(lapply(res, "[[", criterion)))
	# Compare last accepted model to current best model
			if (as.numeric(models[[length(models)]][criterion]) - as.numeric(res[[best]][criterion]) < threshold)
			{
				cat("\nNo significant increase in ", criterion, " score. Disregarding subsequent piecewise models.\n", sep="")
				done <- TRUE
				break;
			}
			models <- c(models, res[best])
			z <- medusa.split(node.list[best], z, anc)$z
			fit <- res[[best]]   # keep track of '$split.at' i.e. nodes already considered
			
			cat("Step ", i+1, ": best likelihood = ", models[[i+1]]$lnLik, "; AICc = ", models[[i+1]]$aicc,
				"; break at node ", models[[i+1]]$split.at[i+1],"\n", sep="")
			i <- i+1
		}
	}
	
	modelSummary <- calculate.model.fit.summary(models, phy, plotFig=ifelse(length(models) > 1 & plotFig, TRUE, FALSE))
	if (verbose)
	{
		cat("\n", "Model fit summary:", "\n\n", sep="")
		print(modelSummary)
	}
	results <- list(z=z.orig, anc=anc, models=models, phy=phy, threshold=threshold, modelSummary=modelSummary)
	
	return(results)
}




summarizeTurboMEDUSA <-
function(results, model=NULL, cutoff="threshold", criterion="aicc", plotTree=TRUE, time=TRUE, node.labels=TRUE, cex=0.5,
	plotSurface=FALSE, n.points=100, ...)
{
# Desirables:
#  1. table listing parameter values of selected model
#  2. list parameters of base model
#  3. tree printed with colour-coded edges, node labels to indicate split position(s)
#  4. plot likelihood surface
	
# Extract constituent components from results
	fit <- results$models
	phy <- results$phy
	z <- results$z
	anc <- results$anc
	modelSummary <- results$modelSummary
	threshold <- results$threshold
	
	fit; phy; z; anc; modelSummary; threshold;
	
# First, determine which model is desired
	model.id <- 0
	if (!is.null(model))
	{
		model.id <- model
	} else {   # Find best model using some criterion (threshold or user-defined)
		if (cutoff != "threshold") {threshold <- cutoff}
		else {cat("\nSelecting model based on corrected threshold (improvement in information theoretic score of ", threshold, " units).\n", sep="")}
		model.id <- 1
		while (1)
		{
			if ((model.id + 1) > length(fit)) break;
			if ((unlist(fit[[model.id]][criterion]) - unlist(fit[[model.id+1]][criterion])) < threshold) break;
			model.id <- model.id + 1
		}
	}
		
	break.pts <- fit[[model.id]]$split.at
	opt.model <- as.data.frame(cbind(Split.node=break.pts, fit[[model.id]]$par, LnLik.part=fit[[model.id]]$lnLik.part))
	base.model <- as.data.frame(fit[[1]]$par)
	
	cat("\nEstimated parameter values for model #", model.id, ":\n\n", sep="")
	print.data.frame(opt.model, digits=5)
	opt.weight <- 0
	opt.fit <- 0
	base.weight <- 0
	base.fit <- 0
	if (criterion == "aicc")
	{
		opt.weight <- modelSummary$w.aicc[model.id]
		opt.fit <- modelSummary$aicc[model.id]
		base.weight <- modelSummary$w.aicc[1]
		base.fit <- modelSummary$aicc[1]
	} else { # aic used
		opt.weight <- modelSummary$w.aic[model.id]
		opt.fit <- modelSummary$aic[model.id]
		base.weight <- modelSummary$w.aic[1]
		base.fit <- modelSummary$aic[1]
	}
	cat("\nModel fit summary for model #", model.id, ":\n\n", sep="")
	cat("\tLog-likelihood = ", as.numeric(results$models[[model.id]]["lnLik"]), "\n", sep="")
	cat("\t", criterion, " = ", opt.fit, "\n", sep="")
	cat("\t", criterion, " weight = ", opt.weight, "\n\n", sep="")
	
	if (model.id != 1)
	{
		if (modelSummary$N.Param[1] == 1) {model <- "Yule"} else {model <- "BD"}
		cat("\nFor comparison, estimated values for the base (single homogeneous-", model, ") model are:\n\n", sep="")
		print.data.frame(base.model, digits=5, row.names=FALSE)
		cat("\nModel fit summary for base model:\n\n", sep="")
		cat("\tLog-likelihood = ", as.numeric(results$models[[1]]["lnLik"]), "\n", sep="")
		cat("\t", criterion, " = ", base.fit, "\n", sep="")
		cat("\t", criterion, " weight = ", base.weight, "\n\n", sep="")
	}
	
# Get desired tree-model conformation
	if (length(break.pts) > 1)
	{
		for (i in 2:length(break.pts))
		{
			tmp <- medusa.split(break.pts[i], z, anc)
			z <- tmp$z
		}
	}
	mm <- integer();
# Plot tree with purdy colours and labelled nodes (to better map between tree and table)
	if (plotTree)
	{
		dev.new()
		margin <- FALSE
		
# This need to be changed to reflect new structure
		mm <- match(phy$edge[,2], z[,"dec"])
		if (time) {margin=TRUE}
		plot.phylo(phy, edge.color=z[mm,"partition"], no.margin=!margin, cex=cex, ...)
		if (time)
		{
			axisPhylo(cex.axis=0.75)
			mtext("Divergence Time (MYA)", at=(max(get("last_plot.phylo", envir = .PlotPhyloEnv)$xx)*0.5), side = 1, line = 2, cex=0.75)
		}
		if (node.labels)
		{
			for (i in  1:length(break.pts))
			{
				nodelabels(i, node= break.pts[i], frame = "c", font = 1, cex=0.5)
			}
		}
	}
	
## Due to my current ineptitude in plotting in R, plotting surface is currently [0,1] for
#  both r and epsilon, even if values of interest are clustered in some subregion.

## *** NEED TO MAKE CHECKS IF CURRENT MODEL IS YULE, AND THEN FIGURE HOW TO PLOT IT ***
	if (plotSurface)
	{
		n.pieces <- length(opt.model[,1])
		dev.new()
		
		if ((sqrt(n.pieces)) == round(sqrt(n.pieces)))  #make square plotting surface
		{
			op <- par(mfrow=c(sqrt(n.pieces),sqrt(n.pieces)))
		} else {
			layout(matrix(1:n.pieces))
		}
		
		if (n.pieces > 1) {cat("Computing surfaces...\n")} else {cat("Computing surface...\n")}
		for (k in 1: n.pieces)
		{
			partition <- k
	## *** Need to indiciate model flavour here
			lik <- make.lik.medusa.part(z[z[,"partition"] == partition,,drop=FALSE], model="bd")
			
			lik.vals <- matrix(nrow=n.points, ncol=n.points)
			r.vals <- seq(from=1e-10, to=1.0, length.out=n.points)
			eps.vals <- seq(from=1e-10, to=1.0, length.out=n.points)
			
			for (i in 1:length(r.vals))
			{
				for (j in 1:(length(eps.vals))) {lik.vals[i,j] <- lik(pars=c(r.vals[i], eps.vals[j]))}
			}
			if (n.pieces > 1) {cat("Completed computing surface for piecewise model #", k, "\n", sep="")}
		
# Contour plot
			max.lik <- as.numeric(opt.model[partition,"LnLik.part"])  # MLE
			lines <- c(max.lik-0.5, max.lik-1, max.lik-2, max.lik-3, max.lik-4, max.lik-5, max.lik-10, max.lik-50, max.lik-100)
			contour(lik.vals, levels=lines, labels=c(0.5, 1, 2, 3, 4, 5, 10, 50, 100), axes=FALSE, xlab="r (b-d)", ylab="epsilon (d/b)")
			tics<-floor(c(1, n.points/4, n.points/2, n.points*3/4, n.points))
			axis(1, at=c(0, 0.25, 0.5, 0.75, 1), labels=round(r.vals[tics], 3))
			axis(2, at=c(0, 0.25, 0.5, 0.75, 1), labels=round(eps.vals[tics], 3))
			points(x=opt.model[partition,"r"], y=opt.model[partition,"epsilon"], pch=16, col="red") # inferred ML values
			legend("topright", "ML", pch=16, col="red", inset = .05, cex=0.75, bty="n")
			if (n.pieces > 1) {title(main=paste("Piecewise model #", k, sep=""))}
		}
	}
	treeParameters <- list(mm=mm, break.pts=break.pts, phy=phy, z=z)
	return(treeParameters)
}




## Function to prune tree using 'richness' information, assumed to have minimally two columns, "taxon" and "n.taxa"
##   Perhaps relax on these column names, may cause too many problems
## May also include 'exemplar' column; in that case, rename relevant tip.label before pruning.
prune.tree.merge.data <- function (phy, richness, verbose)
{
# Rename exemplar taxa with taxon name in richness file
	if (!is.null(richness$exemplar))
	{
# Change relevant tip.labels in phy; individual 'exemplar' may be NA, use original tip.label.
# Ordering in richness file should NOT be assumed to match order of tip.labels
		i.na <- is.na(richness$exemplar)
		phy$tip.label[match(richness$exemplar[!i.na], phy$tip.label)] <- as.character(richness$taxon[!i.na])
	}
	
# make sure things are in the correct order and of correct format
	if (length(richness[1,]) == 2)
	{
		if (colnames(richness)[1] != "taxon" || colnames(richness)[2] != "n.taxa")
		{
			if (class(richness[,1]) == "factor" & class(richness[,2]) == "integer")
			{
				colnames(richness) = c("taxon", "n.taxa")
			} else if (class(richness[,1]) == "integer" & class(richness[,2]) == "factor")
			{
				colnames(richness) = c("n.taxa", "taxon")
			} else {
				cat("turboMEDUSA thinks your richness data is in an incorrect format. See ?runTurboMEDUSA.\n")
				stop
			}
		}
	}
	
# checking for typo; if same size, nothing should be dropped
	check <- FALSE
	if (length(phy$tip.label) == length(richness[,1])) check <- TRUE

# Prune tree down to lineages with assigned richnesses
	temp <- richness[, "n.taxa"]
	names(temp) <- richness[, "taxon"]
	pruned <- treedata(phy, temp, warnings=verbose)  # geiger function calling ape (namecheck)
	if (check) {
		if (length(phy$tip.label) != length(pruned$phy$tip.label)) {
			cat("turboMEDUSA thinks there is a typo in either the tree or richness files.\n")
			stop
		}
	}
	phy <- pruned$phy
# Check the tree
	#	plotNN(phy)					# Node numbers (ape-style) plotted
	
	return(list(phy=phy, richness=richness))
}



## Original default was to fit 20 models (or less if the tree was small).
## Changing to a stop-criterion (stop="model.limit") e.g. when k = n-1 (i.e. when denominator of aicc correction is undefined).
## k <- (3*i-1) # when both birth and death are estimated, where i is the number of piecewise models
  ## This occurs when i = n/3
  ## If Yule, max i = n/2
## n <- (2*num.taxa - 1) == (2*length(richness[,1]) - 1) # i.e. total number of nodes in tree (internal + pendant)
## Alternatively use aicc threshold itself as a stopping criterion (stop="threshold").
get.max.model.limit <- function (richness, model.limit, model, stop, verbose)
{
	samp.size <- (2*length(richness[,1]) - 1)
	if (model == "bd")
	{
		max.model.limit <- as.integer(samp.size/3) - ((!(samp.size %% 3)) * 1)
	} else {
		max.model.limit <- as.integer(samp.size/2) - ((!(samp.size %% 2)) * 1)
	}
	
	if (stop == "model.limit")
	{
		if (model.limit > max.model.limit) {model.limit <- max.model.limit}
	} else if (stop == "threshold") {
		model.limit <- max.model.limit
	}
	
	if (verbose)
	{
		cat("\nLimiting consideration to a maximum of ", model.limit, " piecewise", sep="")
		if (model == "bd") {cat(" BD models")} else {cat(" pure-birth (Yule) models")}
		if (stop == "threshold") {cat(" (or until threshold is not satisfied)")}
		cat(".\n\n")
	}
	
	return(model.limit)
}



## Fitted curve from random b-d simulations
## Value corresponds to 95th percentile of AICc(split) - AICc(no-split) for no-split simulations
## x-shifted power function
get.threshold <- function (x)
{
	a = -3.5941052380332650E+01
	b =  6.7372587299747000E+00
	c = -1.0061508340754866E-01
	Offset =  2.7516678664333408E+01
	y <- a * (x-b)^c + Offset
	if (y < 0) y <- 0
	return(y)
}



## The make.cache.medusa function is like the first half of the original splitEdgeMatrix().
## It works through and reorders the edges, then works out start and end times of these
## based on the phylogeny's branching times.
##
## In addition, every node's ancestors are also calculated.  The element 'anc' is a list.
## $anc[i] contains the indices within $edge, $t.start, etc., of all ancestors of node 'i'
## (in ape node numbering format).
make.cache.medusa <- function (phy, richness)
{
	n.tips <- length(phy$tip.label)
	n.int <- nrow(phy$edge) - n.tips
	
## Ape numbers the tips first
	i.int <- seq_len(n.int)
	interior <- phy$edge[,2] %in% phy$edge[,1]
	bt <- branching.times(phy)
	
# Consider only internal edges first
	edges.int <- phy$edge[interior,]
	colnames(edges.int) <- c("anc", "dec")
	
	t.0 <- bt[match(edges.int[,1], (n.tips+1):max(edges.int))]
	t.1 <- c(t.0[i.int] - phy$edge.length[interior])
	
	z.internal <- cbind(edges.int, t.0, t.1, t.len=t.0 - t.1,
		n.0=rep(1, n.int), n.t=rep(NA, n.int))
	
# Now, pendant edges; 
	edges.pendant <- phy$edge[match(seq_len(n.tips), phy$edge[,2]),]
	colnames(edges.pendant) <- c("anc", "dec")
	
	t.0 <- bt[match(edges.pendant[,1], (n.tips+1):max(edges.pendant))]
	t.1 <- rep(0, n.tips)
# cannot assume richness ordering necessarily matches that of tip labels
	ext.richness <- richness$n.taxa[match(phy$tip.label, richness$taxon)]
	
	z.pendant <- cbind(edges.pendant, t.0, t.1, t.len=t.0 - t.1,
		n.0=rep(1, n.tips), n.t=ext.richness)
	
	z <- rbind(z.internal, z.pendant)
	z <- cbind(z,partition=rep(1, length(z[,1]))) # Stores piecewise model structure
	rownames(z) <- NULL
	
# Used for identifying ancestral nodes below i.e. tracking breakpoints
	all.edges <- as.matrix(z[,c("anc","dec")])
	
	list(z=z, anc=lapply(seq_len(max(all.edges)), ancestors.idx, all.edges))
# And, we're good to go...
}



## This generates the indices of all ancestors of a node, using ape's edge matrix.
## Deals with row numbers of the edge matrix rather than node numbers of the tree.
ancestors <- function (node, all.edges)
{
	ans <- node
	repeat
	{
		node <- all.edges[all.edges[,1] %in% node,2]
		if (length(node) > 0) {ans <- c(ans, node)} else {break}
	}
	unlist(ans)
}



## The function 'ancestors' returns the indices of all ancestors within the edge matrix.
ancestors.idx <- function (node.list, all.edges)
{
	which(all.edges[,1] == node.list | all.edges[,2] %in% ancestors(node.list, all.edges))
}



## Needed for determining whther nodes are virgin nodes
get.num.tips <- function (node, phy)
{
	n <- length(node.leaves(phy,node))
	return(n)
}



## Only used for base model
medusa.ml.initial <- function (z, initial.r, initial.e, model)
{
	rootnode <- min(z[,"anc"])
	obj <- medusa.ml.fit.partition(1, z, sp=c(initial.r, initial.e), model)
	
	model.fit <- calculate.model.fit(fit=obj, z)
	
	list(par=matrix(obj$par, nrow=1, dimnames=list(NULL,c("r", "epsilon"))), lnLik.part=obj$lnLik, 
	   lnLik=obj$lnLik, split.at=rootnode, aic=model.fit[1], aicc=model.fit[2], num.par=model.fit[3])
}



## Pre-fit values for pendant edges; DON'T recalculate later; should account for ~25% of all calculations
medusa.ml.prefit <- function (node, z, anc, initial.r, initial.e, model)
{
	obj <- medusa.split(node, z, anc)
	z <- obj$z
# Partition '2' represents the clade/edge of interest
	fitted <- medusa.ml.fit.partition(2, z, sp=c(initial.r, initial.e), model)
	
	return(fitted)
}



## 'fit' contains parameter values from previous model, used to initialize subsequent model.
## Pass in pre-fitted values for pendant edges and virgin nodes (in 'prefit'); DON'T recalculate.
## Need to consider the possibility of birth-death, yule, or mixed models.
medusa.ml.update <- function (node, z, anc, fit, prefit, num.tips, root.node, model, criterion)
{
	obj <- medusa.split(node, z, anc)
	z <- obj$z
	aff <- obj$affected
	
	op <- fit$par
	sp <- op[aff[1],] # Use previously fit parameter values from clade that is currently being split
	
## In mixed models, want to conserve flavour of previously fit model (right?)
	if (model == "mixed")
	{
		if (sum(!is.na(sp)) < 2)
		{
			fit1 <- medusa.ml.fit.partition(aff[1], z, sp, model="yule")
		} else {
			fit1 <- medusa.ml.fit.partition(aff[1], z, sp, model="bd")
		}
	} else {
		fit1 <- medusa.ml.fit.partition(aff[1], z, sp, model)
	}
	op[aff[1],] <- fit1$par # Replace parameters with new values for diminished clade
	
	fit2 <- NULL
	
## Check if pendant; calculations already done
	if (node < root.node)
	{
		if (model == "bd" | model == "mixed")
		{
			fit2.bd <- prefit$tips$bd[[node]]
		} else {
			fit2.bd <- NULL
		}
		if (model == "yule" | model == "mixed")
		{
			fit2.yule <- prefit$tips$yule[[node]]
		} else {
			fit2.yule <- NULL
		}
## Check if virgin node; save more time!
	} else if (length(unique(z[(z[,"partition"] == aff[2] & z[,"dec"] < root.node),"dec"])) == num.tips[[node]])
	{
		if (model == "bd" | model == "mixed")
		{
			fit2.bd <- prefit$virgin.nodes$bd[[node - root.node]]
		} else {
			fit2.bd <- NULL
		}
		if (model == "yule" | model == "mixed")
		{
			fit2.yule <- prefit$virgin.nodes$yule[[node - root.node]]
		} else {
			fit2.yule <- NULL
		}
## Novel arrangement; need to calculate
 ## Figure out which flavour of model wins below (medusa.ml.fit.partition)
	} else {
		if (model == "bd" | model == "mixed")
		{
			if (is.na(sp[2])) {sp[2] <- 0.5}
			fit2.bd <- medusa.ml.fit.partition(aff[2], z, sp, model="bd")
		} else {
			fit2.bd <- NULL
		}
		if (model == "yule" | model == "mixed")
		{
			fit2.yule <- medusa.ml.fit.partition(aff[2], z, sp, model="yule")
		} else {
			fit2.yule <- NULL
		}
	}
	
## Check which flavour of model is needed; check fit using desired criterion
	if (is.null(fit2.bd) & !is.null(fit2.yule))
	{
		fit2 <- fit2.yule
	} else if (is.null(fit2.yule) & !is.null(fit2.bd)) {
		fit2 <- fit2.bd
	} else {
## Dealing with a 'mixed' model here; need to consider number of parameters.
		bd.fit <- fit
		yule.fit <- fit
		
		bd.fit$par <- rbind(op, fit2.bd$par)
		bd.fit$lnLik.part[aff] <- c(fit1$lnLik, fit2.bd$lnLik)
		bd.fit$lnLik <- sum(bd.fit$lnLik.part)
		bd.model.fit <- calculate.model.fit(bd.fit, z)
		
		yule.fit$par <- rbind(op, fit2.yule$par)
		yule.fit$lnLik.part[aff] <- c(fit1$lnLik, fit2.yule$lnLik)
		yule.fit$lnLik <- sum(yule.fit$lnLik.part)
		yule.model.fit <- calculate.model.fit(yule.fit, z)
		
		if (criterion == "aic") {element <- 1} else {element <- 2}
		if (bd.model.fit[[element]] < yule.model.fit[[element]])
		{
			fit2 <- fit2.bd
		} else {
			fit2 <- fit2.yule
		}
	}
	
	fit$par <- rbind(op, fit2$par)
	fit$lnLik.part[aff] <- c(fit1$lnLik, fit2$lnLik) # Replace parameters with new values for diminished clade
	fit$split.at <- c(fit$split.at, node)
	fit$lnLik <- sum(fit$lnLik.part)
	
	model.fit <- calculate.model.fit(fit, z)
	
	fit$aic <- model.fit[1]
	fit$aicc <- model.fit[2]
	fit$num.par <- model.fit[3]
	
	return(fit)
}



## Split the edge matrix 'z' by adding a partition rooted at node 'node'.
##   Note: in original MEDUSA parlance, this is cutAtStem=T.
## The list 'anc' is a list of ancestors (see make.cache.medusa, above).
## Returns a list with elements:
##   z: new medusa matrix, with the new partition added
##   affected: indices of the partitions affected by the split (n == 2).
medusa.split <- function (node, z, anc)
{
	part <- z[,"partition"]
	base <- min(part[z[,1] == node | z[,2] == node])
	tag <- max(part) + 1

	i <- anc[[node]]
	idx <- i[part[i] == base]
	z[idx,"partition"] <- tag
	
	z[which(z["dec"] == node),"partition"] <- tag # Possible to have several edges to consider

	list(z=z, affected=c(unique(part[idx]), tag))
}



## sp = initializing values for r & epsilon
## Default values should never be used (except for first model), as the values from the previous model are passed in
medusa.ml.fit.partition <- function (partition, z, sp=c(0.1, 0.05), model)
{
# Construct likelihood function:
	lik <- make.lik.medusa.part(z[z[,"partition"] == partition,,drop=FALSE], model)
	
	if (model == "bd")
	{
		fit <- optim(fn=lik, par=sp, method="N", control=list(fnscale=-1)) # last argument connotes maximization
		list(par=fit$par, lnLik=fit$value)
	} else {
		fit <- optimize(f=lik, interval=c(0, 1), maximum=TRUE)
		par <- c(fit$maximum, NA)
		list(par=par, lnLik=fit$objective)
	}
#	list(par=fit$par, lnLik=fit$value)
}



## make.lik.medusa.part: generate a likelihood function for a single partition.
make.lik.medusa.part <- function (partition, model)
{

# Handle internal and pendant edges separately
	is.int <- is.na(partition[,"n.t"])
	is.pend <- !is.int
	
	n.int <- sum(is.int)
	n.pend <- sum(is.pend)
	
	if (n.int + n.pend != length(partition[,1])) stop("You messed up, yo.")
	
## Internal and pendant calculations differ; split'em up
	int  <- partition[is.int,,drop=FALSE]
	pend <- partition[is.pend,,drop=FALSE]
	
	sum.int.t.len <- sum(int[,"t.len"])  # Simply sum all internal edges
	int.t.0 <- int[,"t.0"]
	
# 'n.0' = Foote's 'a', initial diversity; 'n.t' = Foote's 'n', final diversity
	pend.n.0 <- pend[,"n.0"] # Foote's 'a': initial diversity
	pend.n.t <- pend[,"n.t"] # Foote's 'n': final diversity
	pend.t.len <- pend[,"t.len"]
	
# User may pass in epsilon; don't change it, just estimate r
	f <- function(pars)
	{
		if (model == "bd")
		{
			r <- pars[1]
			epsilon <- pars[2]
			
			if (r < 0 | epsilon <= 0 | epsilon >= 1) {return(-Inf)}
		} else if (model == "yule") {
			r <- pars[1]
			epsilon <- 0
			
			if (r < 0) {return(-Inf)}
		}
			
#		if (r < 0 | epsilon <= 0 | epsilon >= 1) {return(-Inf)}
		
		if (n.int == 0) {l.int <- 0} else {
## Likelihood of internal edges from Rabosky et al. (2007) equation (2.3):
			l.int <- n.int * log(r) - r * sum.int.t.len - sum(log(1 - (epsilon * exp(-r * int.t.0))))
		}
		
		if (n.pend == 0) {l.pend <- 0} else {
## Calculations are from the following:
## Rabosky et al. 2007. Proc. Roy. Soc. 274: 2915-2923.
## Foote et al. 1999. Science. 283: 1310-1314
## Raup. 1985. Paleobiology 11: 42-52 [Foote et al. correct the equation [A18] where a > 1]
## Bailey. 1964. The Elements Of Stochastic Processes, With Applications To The Natural Sciences
## Kendall. 1948. Ann. Math. Stat. 19: 1–15.
##
## A = probability of extinction of one lineage over time 't'
## B = A * (lambda/mu)
##
## When there is a single lineage at time 0 (a = 1), the calculation is
##   log(1 - A) + log(1 - B) + (n - 1)*log(B)
## but this is conditioned on survival by dividing by (1-A)
## (subtracting log(1 - A) on a log scale) which cancels to give:
##   log(1 - B) + (n - 1)*log(B)
##      - for n.t == 1, reduces further to log(1-B)
##
## A = mu*(exp((lambda - mu)*t) - 1)) / (lambda*exp((lambda - mu)*t) - mu)
##  let r = (lambda - mu); ert = exp((lambda - mu)*t)
## A = mu*(ert - 1)/(lambda*ert - mu)
##
## B = A * (lambda/mu)
##   = [mu*(ert - 1)/(lambda*ert - mu)] * (lambda/mu)
##   = (lambda*(ert - 1))/(lambda*ert - mu)
##   = (lambda*(ert - 1))/(lambda(ert - mu/lambda))
##   = (ert - 1) / (ert - epsilon)

## All pendant nodes begin with richness '1'; calculations simple.
#			i.pend.n.t.1 <- which(pend.n.t == 1)   # calculations even simpler: log(1-B)
#			i.pend.n.t.n1 <- which(pend.n.t != 1)
			
			ert <- exp(r * pend.t.len)
			B <- (ert - 1) / (ert - epsilon) # Equivalently: B <- (bert - b) / (bert - d)
			
			l.pend <- sum(log(1 - B) + (pend.n.t - 1)*log(B))
		}
		l.int + l.pend
	}
}



## 'fit' contains '$par' and '$lnlik'
calculate.model.fit <- function (fit, z)
{
## Sample size taken (for now) as the total num.nodes in the tree (internal + pendant)
  # num.nodes = (2*length(phy$tip.label) - 1) == (2*length(richness[,1]) - 1) == length(z[,1]) + 1
#	n <- (length(z[,1]) + 1) + sum(!is.na(z[,"n.f"]))
	
# Since each edge defines a node (i.e. an 'observation'), need only add root node as final obervation
	n <- (length(z[,1]) + 1)
	
 # Includes both formal parameters AND number of breaks. Note: first model does not involve a break.
## Models where all parameters are estimated (i.e. BD model):
  # 2 parameters for base model (no breakpoint) + 3 parameters (r, eps, breakpoint) for each subsequent model
  
  
# Determine number of piecewise models currently involved
	if (length(fit$par) < 3) # i.e. base model
	{
		num.models <- 1
	} else {
		num.models <- length(fit$par[,1])
	}
	
# Updated for more general models: check how many parameter values != NA
#	k <- 2 + (3 * (num.models - 1))
	k <- sum(!is.na(fit$par)) + (num.models - 1) # number of estimated parameters + number of breaks
	
	lnLik <- fit$lnLik
	
	aic <- (-2 * lnLik) + (2*k)
	aicc <- aic + 2*k*(k+1)/(n-k-1)
	
	model.fit <- c(aic, aicc, k)
	return(model.fit)
}



## Prints out a table of likelihoods, parameters, aic scores, and aic weights (delta-aics are also available, if desired)
calculate.model.fit.summary <- function (models, phy, plotFig, fig.title=NULL, ...)
{
	tmp <- matrix(nrow=(length(models)), ncol=6)
	colnames(tmp) <- c("N.Models", "Break.Node", "Ln.Lik", "N.Param", "aic", "aicc")
	
	w.aic <- numeric(length(models))
	w.aicc <- numeric(length(models))
	
	for (i in 1:length(tmp[,1]))
	{
		tmp[i,] <- c(i, as.integer(models[[i]]$split.at[i]), models[[i]]$lnLik, models[[i]]$num.par, models[[i]]$aic, models[[i]]$aicc)
	}
	
	all.res <- as.data.frame(tmp)
	all.res[1,2] <- NA # root node for base model
	
	w.aic <- calculate.model.weights(all.res$aic)
	w.aicc <- calculate.model.weights(all.res$aicc)
	
	all.res <- cbind(all.res[,c(1:5)], w.aic=w.aic$w, aicc=all.res$aicc, w.aicc=w.aicc$w)
	
	if (plotFig)
	{
		dev.new()
		plotModelFit(all.res)
		if (!is.null(fig.title)) {title(main=fig.title, cex.main=0.75)}
	}
	return(all.res)
}



## Self explanatory
calculate.model.weights <- function (fit)
{
	best <- min(fit)
	delta <- fit-best
	sumDelta <- sum(exp(-0.5 * delta))
	w <- (exp(-0.5 * delta)/sumDelta)
	
	results <- data.frame(fit=fit,delta=delta,w=w)
	
	return(results)
}



## Create a plot of model-fit vs. model-size
plotModelFit <- function (all.res)
{
	ylim <- c(min(all.res[,"aic"],all.res[,"aicc"]), max(all.res[,"aic"],all.res[,"aicc"]))
	plot(all.res[,"N.Models"],all.res[,"aicc"], xlab="Number of Piecewise Models", ylab="Model Fit", ylim=ylim, type="l", col="blue")
	points(all.res[,"N.Models"],all.res[,"aicc"], col="blue", pch=21, bg="white")
	points(all.res[,"N.Models"],all.res[,"aic"], col="black", type="l")
	points(all.res[,"N.Models"],all.res[,"aic"], col="black", pch=21, bg="white")
	
	legend("topleft", c("aicc","aic"), pch=21, pt.bg="white", lty=1, col=c("blue", "black"), inset = .05, cex=0.75, bty="n") # 'bottomright' also works
}


# treeParameters <- list(mm=mm, break.pts=break.pts, phy=phy, z=z)
plotPrettyTree <- function (treeParameters, time=TRUE, node.labels=FALSE, cex=0.5, ...)
{
	mm <- treeParameters$mm
	break.pts <- treeParameters$break.pts
	phy <- treeParameters$phy
	z <- treeParameters$z
	
	dev.new()
	margin <- FALSE
	
# This need to be changed to reflect new structure
	mm <- match(phy$edge[,2], z[,"dec"])
	if (time) {margin=TRUE}
	plot.phylo(phy, edge.color=z[mm,"partition"], no.margin=!margin, cex=cex, ...)
	if (time)
	{
		axisPhylo(cex.axis=0.75)
		mtext("Divergence Time (MYA)", at=(max(get("last_plot.phylo", envir = .PlotPhyloEnv)$xx)*0.5), side = 1, line = 2, cex=0.75)
	}
	if (node.labels)
	{
		for (i in  1:length(break.pts))
		{
			nodelabels(i, node= break.pts[i], frame = "c", font = 1, cex=0.5)
		}
	}
}




## Get b and d values from r (b-d) and epsilson (d/b)
## Used in previous version of program; now in terms of r and epsilon
## Possibly of use to users wishing to translate results
get.b.d <- function (r, epsilon)
{
	b <- r/(1-epsilon)
	d <- b-r   # Alternatively: d <- eps*r/(1-eps)
	return(list(b=b, d=d))
}

## Print out tree with ape-style node-numbering
## Possibly of interest for users to identify numbers of node(s) off interest
 ## If this is the case, make sure to pass in pruned tree
plotNN <- function (phy, time=TRUE, margin=TRUE, label.offset=0.5, cex=0.5, ...) 
{
	phy$node.label <- (length(phy$tip.label) + 1):max(phy$edge)
	plot.phylo(phy, show.node.label=TRUE, no.margin=!margin, label.offset=label.offset, cex=cex, ...)
	if (time && !margin) cat("Cannot plot time axis without a margin.\n")
	else if (time && margin) axisPhylo(cex.axis=0.75)
}