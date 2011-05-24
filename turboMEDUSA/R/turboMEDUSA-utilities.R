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
	
	if (length(phy$tip.label) != length(richness[,1]))
	{
# Prune tree down to lineages with assigned richnesses
		temp <- richness[, "n.taxa"]
		names(temp) <- richness[, "taxon"]
		pruned <- treedata(phy, temp, warnings=verbose)  # geiger function calling ape (namecheck)
		phy <- pruned$phy
# Check the tree
	#	plotNN(phy)					# Node numbers (ape-style) plotted
	}
	return(list(phy=phy, richness=richness))
}



## Original default was to fit 20 models (or less if the tree was small).
## Changing to a stop-criterion (stop="model.limit") e.g. when k = n-1 (i.e. when denominator of aicc correction is undefined).
## k <- (3*i-1) # when both birth and death are estimated, where i is the number of piecewise models
  ## This occurs when i = n/3
## n <- (2*num.taxa - 1) == (2*length(richness[,1]) - 1) # i.e. total number of nodes in tree (internal + pendant)
## Alternatively use aicc itself as a stopping criterion (stop="threshold").
get.max.model.limit <- function (richness, model.limit, stop, verbose)
{
	samp.size <- (2*length(richness[,1]) - 1)
	
	max.model.limit <- as.integer(samp.size/3) - ((!(samp.size %% 3)) * 1)
	if (stop == "model.limit")
	{
		if (model.limit > max.model.limit) {model.limit <- max.model.limit}
		if (verbose) {cat("\nLimiting consideration to ", model.limit, " piecewise BD models\n\n", sep="")}
	} else if (stop == "threshold") {
		model.limit <- max.model.limit
		if (verbose) {cat("\nConsidering a maximum of ", model.limit, " piecewise BD models (or until threshold is not satisfied)\n\n", sep="")}
	}
	
	return(model.limit)
}



## Fitted curve from random b-d simulations
## Value corresponds to 95th percentile of AICc(split) - AICc(no-split) for no-split simulations
get.threshold <- function (x)
{
	a = -3.5941052380332650E+01
	b =  6.7372587299747000E+00
	c = -1.0061508340754866E-01
	Offset =  2.7516678664333408E+01
	y <- a * (x-b)^c + Offset
	return(y)
}



## The make.cache.medusa function is like the first half of the original splitEdgeMatrix().
## It works through and reorders the edges, then works out start and end times of these
## based on the phylogeny's branching times.
##
## In addition, every node's ancestors are also calculated.  The element 'anc' is a list.
## $anc[i] contains the indices within $edge, $t.start, etc., of all ancestors of node 'i'
## (in ape node numbering format).
make.cache.medusa <- function (phy, richness, mc, num.cores)
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
	
	if (mc)
	{
		list(z=z, anc=mclapply(seq_len(max(all.edges)), ancestors.idx, all.edges), mc.cores=num.cores)
	} else {
		list(z=z, anc=lapply(seq_len(max(all.edges)), ancestors.idx, all.edges))
	} # And, we're good to go...
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
medusa.ml.initial <- function (z, initial.r, initial.e)
{
	rootnode <- min(z[,"anc"])
	obj <- medusa.ml.fit.partition(1, z, sp=c(initial.r, initial.e))
	
	model.fit <- calculate.model.fit(fit=obj, z)
	
	list(par=matrix(obj$par, nrow=1, dimnames=list(NULL,c("r", "epsilon"))), lnLik.part=obj$lnLik, 
	   lnLik=obj$lnLik, split.at=rootnode, aic=model.fit[1], aicc=model.fit[2], num.par=model.fit[3])
}



## Pre-fit values for pendant edges; DON'T recalculate later; should account for ~25% of all calculations
medusa.ml.prefit <- function (node, z, anc, initial.r, initial.e)
{
	obj <- medusa.split(node, z, anc)
	z <- obj$z
# Partition '2' represents the clade/edge of interest
	fitted <- medusa.ml.fit.partition(2, z, sp=c(initial.r, initial.e))
	
	return(fitted)
}



## 'fit' contains parameter values from previous model, used to initialize subsequent model.
## Pass in pre-fitted values for pendant edges and virgin nodes; DON'T recalculate.
medusa.ml.update <- function (node, z, anc, fit, tips, virgin.nodes, num.tips, root.node)
{
	obj <- medusa.split(node, z, anc)
	z <- obj$z
	aff <- obj$affected
	
	op <- fit$par
	sp <- op[aff[1],] # Use previously fit parameter values from clade that is currently being split
	
	fit1 <- medusa.ml.fit.partition(aff[1], z, sp)
	fit2 <- 0
	
## Check if pendant; calculations already done
	if (node < root.node)
	{
		fit2 <- tips[[node]]
## Check if virgin node; save more time!
	} else if (length(unique(z[(z[,"partition"] == aff[2] & z[,"dec"] < root.node),"dec"])) == num.tips[[node]])
	{
		fit2 <- virgin.nodes[[node - root.node]]
## Novel arrangement; need to calculate
	} else {
		fit2 <- medusa.ml.fit.partition(aff[2], z, sp)
	}
	
	op[aff[1],] <- fit1$par # Replace parameters with new values for diminished clade
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
medusa.ml.fit.partition <- function (partition, z, sp=c(0.1, 0.05))
{
# Construct likelihood function:
	lik <- make.lik.medusa.part(z[z[,"partition"] == partition,,drop=FALSE])
	
	fit <- optim(fn=lik, par=sp, method="N", control=list(fnscale=-1))
	list(par=fit$par, lnLik=fit$value)
}



## make.lik.medusa.part: generate a likelihood function for a single partition.
make.lik.medusa.part <- function (partition)
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
		r <- pars[1]
		epsilon <- pars[2]
			
		if (r < 0 | epsilon <= 0 | epsilon >= 1) {return(-Inf)}
		
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
## Models where all parameters are estimated:
  # 2 parameters for base model (no breakpoint) + 3 parameters (r, eps, breakpoint) for each subsequent model
	
	if (length(fit$par) < 3) # i.e. base model
	{
		num.models <- 1
	} else {
		num.models <- length(fit$par[,1])
	}
	
	k <- 2 + (3 * (num.models - 1))
	
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



## Get b and d values from r (b-d) and epsilson (d/b)
## Used in previous version of program; now in terms of r and epsilon
## Possibly of use to users wishing to translate results
get.b.d <- function (r, epsilon)
{
	b <- r/(1-epsilon)
	d <- b-r   # Alternatively: d <- eps*r/(1-eps)
	return(list(b=b, d=d))
}