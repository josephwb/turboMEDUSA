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
	
# checking for typo; if same size, nothing should be dropped
	check <- FALSE;
	if (length(phy$tip.label) == length(richness[,1])) check <- TRUE;

# Prune tree down to lineages with assigned richnesses
	temp <- richness[, "n.taxa"];
	names(temp) <- richness[, "taxon"];
	pruned <- treedata(phy, temp, warnings=verbose)  # geiger function calling ape (namecheck)
	if (check) {
		if (length(phy$tip.label) != length(pruned$phy$tip.label)) {
			cat("turboMEDUSA thinks there is a typo in either the tree or richness files.\n");
			stop
		}
	}
	phy <- pruned$phy;
# Check the tree
	#	plotNN(phy)					# Node numbers (ape-style) plotted
	
	return(list(phy=phy, richness=richness));
}


## Original default was to fit 20 models (or less if the tree was small).
## Changing to a stop-criterion (stop="model.limit") e.g. when k = n-1 (i.e. when denominator of aicc correction is undefined).
## k <- (3*i-1) # when both birth and death are estimated, where i is the number of piecewise models
  ## This occurs when i = n/3
  ## If Yule, max i = n/2
## n <- (2*num.taxa - 1) == (2*length(richness[,1]) - 1) # i.e. total number of nodes in tree (internal + pendant)
## Alternatively use aicc threshold itself as a stopping criterion (stop="threshold").
# AICc = AIC + 2*k*(k+1)/(n-k-1);
get.max.model.limit <- function (richness, model.limit, model, stop, verbose)
{
	samp.size <- (2*length(richness[,1]) - 1)
	if (model == "bd" || model == "mixed")
	{
		max.model.limit <- as.integer(samp.size/3) - ((!(samp.size %% 3)) * 1);
	} else {
		max.model.limit <- as.integer(samp.size/2) - ((!(samp.size %% 2)) * 1);
	}
	
	if (stop == "model.limit")
	{
		if (model.limit > max.model.limit) {model.limit <- max.model.limit;}
	} else {
		model.limit <- max.model.limit;
	}
	
	if (verbose)
	{
		cat("\nLimiting consideration to a maximum of ", model.limit, " piecewise", sep="");
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
	
	return(model.limit);
}


## Fitted curve from random b-d simulations
## Value corresponds to 95th percentile of AICc(split) - AICc(no-split) for no-split simulations
## x-shifted power function
get.threshold <- function (x)
{
	a = -3.5941052380332650E+01;
	b =  6.7372587299747000E+00;
	c = -1.0061508340754866E-01;
	Offset =  2.7516678664333408E+01;
	y <- a * (x-b)^c + Offset;
	if (y < 0) y <- 0;
	return(y);
}


## The make.cache.medusa function is like the first half of the original splitEdgeMatrix().
## It works through and reorders the edges, then works out start and end times of these
## based on the phylogeny's branching times.
##
## In addition, every node's descendants are also calculated.  The element 'desc' is a list.
## $desc[i] contains the indices within $edge, $t.start, etc., of all descendants of node 'i'
## (in ape node numbering format).
make.cache.medusa <- function (phy, richness, all.nodes, mc, num.cores)
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
		desc.stem <- mclapply(seq_len(max(all.edges)), descendants.cutAtStem.idx, all.edges=all.edges, mc.cores=num.cores);
		desc.node <- mclapply(seq_len(max(all.edges)), descendants.cutAtNode.idx, all.edges=all.edges, mc.cores=num.cores);
	} else {
		desc.stem <- lapply(seq_len(max(all.edges)), descendants.cutAtStem.idx, all.edges=all.edges);
		desc.node <- lapply(seq_len(max(all.edges)), descendants.cutAtNode.idx, all.edges=all.edges);
	}
	
	for (i in 1:length(desc.node)) # vectorize this shit!
	{
		if (length(desc.node[[i]]) == 0) # tips
		{
			desc.node[[i]] = descendants.cutAtStem.idx(node.list=i, all.edges=all.edges)
		}
	}
	
## Needed downstream; don't recalculate
 ## Gives the number of tips associated with an internal node; determines whether a node is 'virgin' or not
	num.tips <- list()
	if (mc)
	{
		num.tips <- mclapply(all.nodes, get.num.tips, phy=phy, mc.cores=num.cores);
	} else {
		num.tips <- lapply(all.nodes, get.num.tips, phy=phy);
	}
	
	res <- list(z=z, desc.stem=desc.stem, desc.node=desc.node, num.tips=num.tips);
	return(res);
}


## Needed for determining whether nodes are virgin nodes
get.num.tips <- function (node, phy)
{
	n <- length(node.leaves(phy,node));
	return(n);
}


## This generates the indices of all descendants of a node, using ape's edge matrix.
## Deals with row numbers of the edge matrix rather than node numbers of the tree.
descendants.cutAtStem <- function (node, all.edges)
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
descendants.cutAtStem.idx <- function (node.list, all.edges)
{
	which(all.edges[,1] == node.list | all.edges[,2] %in% descendants.cutAtStem(node.list, all.edges));
}

## This generates the indices of all descendants of a node, using ape's edge matrix.
## Deals with row numbers of the edge matrix rather than node numbers of the tree.
descendants.cutAtNode <- function (node, all.edges)
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
descendants.cutAtNode.idx <- function (node.list, all.edges)
{
	which(all.edges[,1] == node.list | all.edges[,2] %in% descendants.cutAtNode(node.list, all.edges));
}


## Only used for base model
medusa.ml.initial <- function (z, initialR, initialE, model)
{
	rootnode <- min(z[,"anc"]);
	obj <- medusa.ml.fit.partition(partition=1, z=z, sp=c(initialR, initialE), model=model);
	
	model.fit <- calculate.model.fit(fit=obj, z=z);
	
	return(list(par=matrix(obj$par, nrow=1, dimnames=list(NULL,c("r", "epsilon"))), lnLik.part=obj$lnLik, 
	   lnLik=obj$lnLik, split.at=rootnode, aic=round(model.fit[1], digits=7), aicc=round(model.fit[2], digits=7), num.par=model.fit[3], cut.at="node"));
}


## Pre-fit values for pendant edges; DON'T recalculate later; should account for ~25% of all calculations
## Also cache values for virgin nodes; useful until subsetted.

## sp = initializing values for r & epsilon
## Default values should never be used (except for first model), as the values from the previous model are passed in
medusa.ml.fit.partition <- function (partition, z, sp=c(0.1, 0.05), model)
{
# Construct likelihood function:
	lik <- make.lik.medusa.part(partition=(z[z[,"partition"] == partition,,drop=FALSE]), model=model);	
	foo <- function (x) {-lik(pars=exp(x));} # work with parameters in log-space to preserve precision
	
	if (model == "bd")
	{
		fit <- optim(fn=foo, par=log(sp), method="N");
		return(list(par=exp(fit$par), lnLik=-fit$value));
	} else {
		fit <- optimize(f=foo, interval=c(-25, 1));
		par <- c(exp(fit$minimum), NA);
		return(list(par=par, lnLik=-fit$objective));
	}
}



## sp = initializing values for r & epsilon
## Default values should never be used (except for first model), as the values from the previous model are passed in
# medusa.ml.fit.partition <- function (partition, z, sp=c(0.1, 0.05), model)
# {
# # Construct likelihood function:
	# lik <- make.lik.medusa.part(partition=(z[z[,"partition"] == partition,,drop=FALSE]), model=model);	
	# foo <- function (sp) {-lik(pars=exp(sp));} # work with parameters in log-space to preserve precision
	
	# if (model == "bd")
	# {
		# fit <- optim(fn=foo, par=log(sp), method="N");
		
		# return(list(par=exp(fit$par), lnLik=-fit$value));
	# } else {
		
# # use age and tip richness to inform interval
		# interval <- NULL;
		# lower <- -30;
		# upper <- 1;
		# node.richness <- sum(z[z[,"partition"] == partition,"n.t"], na.rm=TRUE);
		# # if (node.richness > 0)
		# # {
			# # sum.time <- sum(z[z[,"partition"] == partition,"t.len"], na.rm=TRUE);
	# # #		best.guess <- (node.richness - 2 + (1 + 2^-50)) / (sum.time + (1 + 2^-50));
			
			
			# # best.guess <- (log(node.richness) + sqrt(.Machine$double.eps))/(sum.time + sqrt(.Machine$double.eps));
# # #			cat("\nbest guess: ", best.guess, "\n")
# # #			cat("richness:", node.richness, "; time: ", sum.time, "\n", sep="");
# # #			cat("richness: ", z[z[,"partition"] == partition, "n.t"], "\n")
# # #			cat("anc: ", z[z[,"partition"] == partition, "anc"], "\n")
			
			# # lower <- log(best.guess/1000);
			# # upper <- log(best.guess*10);
	# # #		upper <- log(min((best.guess*10),1));
			# # interval <- c(lower, upper);
		# # } else {interval <- c(-30, 1);}
		# interval <- c(-30, 1);
# #		fit <- optimize(f=vf, interval=interval);
		# fit <- optimize(f=foo, interval=interval);
# #		fit <- optimize(f=foo, interval=c(-30, 1));
		# par <- c(exp(fit$minimum), NA);
		
		# if (exp(lower) > par[1] || exp(upper) < par[1])
		# {
			# cat("\nrichness:", node.richness, "; time: ", sum.time, "\n", sep="");
			# cat("lower: ", exp(lower), "; upper:", exp(upper), "; guess: ", best.guess, "\n", sep="");
			# cat("interval: ", exp(interval), "\n");
			# cat("anc: ", z[z[,"partition"] == partition, "anc"], "\n")
			# cat("dec: ", z[z[,"partition"] == partition, "dec"], "\n")
			# cat("LnLik: ", -fit$objective, "\n")
			# cat(par, "\n");
		# }
		
		# # if (lower < exp(-25) || upper > exp(1))
		# # {
			# # cat("\n", node.richness, sum.time, best.guess, "\n");
			# # cat(interval, "\n");
			# # cat(par, "\n");
		# # }
		
		# return(list(par=par, lnLik=-fit$objective));
	# }
# }



## Split the edge matrix 'z' by adding a partition rooted at node 'node'.
##   Note: in original MEDUSA parlance, this is cutAtStem=T.
## The list 'desc' is a list of descendants (see make.cache.medusa, above).
## Returns a list with elements:
##   z: new medusa matrix, with the new partition added
##   affected: indices of the partitions affected by the split (n == 2).
## This is where 'shiftCut' matters
medusa.split <- function (node, z, desc, shiftCut)
{
	descendants <- NULL;
	if (shiftCut == "stem")
	{
		descendants <- desc$desc.stem;
	} else {
		descendants <- desc$desc.node;
	}
	part <- z[,"partition"];
	base <- min(part[z[,"anc"] == node | z[,"dec"] == node]);
	tag <- max(part) + 1;
	i <- descendants[[node]];
	idx <- i[part[i] == base];
	z[idx,"partition"] <- tag;
	z[which(z["dec"] == node),"partition"] <- tag; # Possible to have several edges to consider
	return(list(z=z, affected=c(unique(part[idx]), tag)));
}



## 'fit' contains parameter values from previous model, used to initialize subsequent model.
## Pass in pre-fitted values for pendant edges and virgin nodes (in 'prefit'); DON'T recalculate.
## Need to consider the possibility of birth-death, yule, or mixed models.
## Need to consider where shft is placed (shiftCut). Placement affects not only new clade, but
## also the size of the split clade. Only relevant if shiftCut = "both".
## fit1 model is already logged; only need to record fit2 model, and only non-prefitted nodes
medusa.ml.update <- function (node, z, desc, fit, prefit, num.tips, root.node, model, criterion, shiftCut)
{
## various combinations possible
	fit1.stem <- NULL;
	fit1.node <- NULL;
	fit2.stem <- NULL;
	fit2.node <- NULL;
	cut.at <- NULL;
	
	sp <- NULL;
	aff <- NULL;
	op <- fit$par;
	cut.at <- NULL;
	
	fit1 <- NULL;
	fit2 <- NULL;
	
	if (shiftCut == "stem" || shiftCut == "both")
	{
## First, diminshed clade
		obj.stem <- medusa.split(node=node, z=z, desc=desc, shiftCut="stem");
		z.stem <- obj.stem$z;
		aff <- obj.stem$affected;
		sp <- op[aff[1],]; # Use previously fit parameter values from clade that is currently being split
		
		if (model == "mixed") ## In mixed models, want to conserve flavour of previously fit model (right?)
		{
			if (sum(!is.na(sp)) < 2) # yule
			{
				fit1.stem <- medusa.ml.fit.partition(partition=aff[1], z=z.stem, sp=sp, model="yule");
			} else {
				fit1.stem <- medusa.ml.fit.partition(partition=aff[1], z=z.stem, sp=sp, model="bd");
			}
		} else {
			fit1.stem <- medusa.ml.fit.partition(partition=aff[1], z=z.stem, sp=sp, model=model);
			fit1.stem$model <- model;
		}
## Second, new clade
		if (node < root.node) # tip, already calculated
		{
			fit2.stem <- prefit$tips[[node]];
		} else if (length(unique(z.stem[(z.stem[,"partition"] == aff[2] & z.stem[,"dec"] < root.node),"dec"])) == num.tips[[node]]) {
			fit2.stem <- prefit$virgin.nodes$stem[[node - root.node]];
		} else { # novel shift
			fit2.stem.bd <- NULL;
			fit2.stem.yule <- NULL;
			
			if (model == "yule" || model == "mixed")
			{
				fit2.stem.yule <- medusa.ml.fit.partition(aff[2], z.stem, sp=sp, model="yule");
				fit2.stem.yule$model <- "yule";
			}
			if (model == "bd" || model == "mixed")
			{
				if (is.na(sp[2])) {sp[2] <- 0.5;}
				fit2.stem.bd <- medusa.ml.fit.partition(aff[2], z.stem, sp=sp, model="bd");
				fit2.stem.bd$model <- "bd";
			}
## Figure out which model fits best
			if (is.null(fit2.stem.bd))
			{
				fit2.stem <- fit2.stem.yule;
			} else if (is.null(fit2.stem.yule)) {
				fit2.stem <- fit2.stem.bd;
			} else {
## Considering both places for a shift
				fit2.stem.bd.val <- calculate.model.fit(fit=fit2.stem.bd, z=z);
				fit2.stem.yule.val <- calculate.model.fit(fit=fit2.stem.yule, z=z);
				
				if (criterion == "aic") {element <- 1;} else {element <- 2;}
				
				if (fit2.stem.bd.val[[element]] < fit2.stem.yule.val[[element]])
				{
					fit2.stem <- fit2.stem.bd;
				} else {
					fit2.stem <- fit2.stem.yule;
				}
			}
		}
	}
	if ((shiftCut == "node" || shiftCut == "both") && (node >= root.node))
	{
## First, diminshed clade
		obj.node <- medusa.split(node=node, z=z, desc=desc, shiftCut="node");
		z.node <- obj.node$z;
		aff <- obj.node$affected;
		sp <- op[aff[1],]; # Use previously fit parameter values from clade that is currently being split
		
		if (model == "mixed") ## In mixed models, want to conserve flavour of previously fit model (right?)
		{
			if (sum(!is.na(sp)) < 2) # yule
			{
				fit1.node <- medusa.ml.fit.partition(partition=aff[1], z=z.node, sp=sp, model="yule");
				fit1.node$model <- "yule";
			} else {
				fit1.node <- medusa.ml.fit.partition(partition=aff[1], z=z.node, sp=sp, model="bd");
				fit1.node$model <- "bd";
			}
		} else {
			fit1.node <- medusa.ml.fit.partition(partition=aff[1], z=z.node, sp=sp, model=model);
			fit1.node$model <- model;
		}
## Second, new clade; check if virgin node (tips cannot be cut at node)
		if (length(unique(z.node[(z.node[,"partition"] == aff[2] & z.node[,"dec"] < root.node),"dec"])) == num.tips[[node]]) # virgin node
		{
			fit2.node <- prefit$virgin.nodes$node[[node - root.node]];
		}
		else { # novel shift
			fit2.node.bd <- NULL;
			fit2.node.yule <- NULL;
			
			if (model == "yule" || model == "mixed")
			{
				fit2.node.yule <- medusa.ml.fit.partition(aff[2], z.node, sp, model="yule");
				fit2.node.yule$model <- "yule";
			}
			if (model == "bd" || model == "mixed")
			{
				if (is.na(sp[2])) {sp[2] <- 0.5;}
				fit2.node.bd <- medusa.ml.fit.partition(aff[2], z.node, sp, model="bd");
				fit2.node.bd$model <- "bd";
			}
## Figure out which model fits best
			if (is.null(fit2.node.bd))
			{
				fit2.node <- fit2.node.yule;
			} else if (is.null(fit2.node.yule)) {
				fit2.node <- fit2.node.bd;
			} else {
## Considering both places for a shift
				fit2.node.bd.val <- calculate.model.fit(fit=fit2.node.bd, z=z);
				fit2.node.yule.val <- calculate.model.fit(fit=fit2.node.yule, z=z);
				
				if (criterion == "aic") {element <- 1;} else {element <- 2;}
				
				if (fit2.node.bd.val[[element]] < fit2.node.yule.val[[element]])
				{
					fit2.node <- fit2.node.bd;
				} else {
					fit2.node <- fit2.node.yule;
				}
			}	
		}
	}
	
## Now, figure out which shift position is optimal	
	if (is.null(fit2.node))
	{
		fit1 <- fit1.stem;
		fit2 <- fit2.stem;
		cut.at <- "stem";
	} else if (is.null(fit1.stem)) {
		fit1 <- fit1.node;
		fit2 <- fit2.node;
		cut.at <- "node";
	} else {
## Considering both places for a shift
		stem.lik <- (fit1.stem$lnLik + fit2.stem$lnLik);
		stem.par <- rbind(fit1.stem$par, fit2.stem$par)
		stem.val <- list(lnLik=stem.lik, par=stem.par);
		stem.fit <- calculate.model.fit(fit=stem.val, z=z);
		
		node.lik <- (fit1.node$lnLik + fit2.node$lnLik);
		node.par <- rbind(fit1.node$par, fit2.node$par)
		node.val <- list(lnLik=node.lik, par=node.par);
		node.fit <- calculate.model.fit(fit=node.val, z=z);
		
		element <- NULL;
		if (criterion == "aic") {element <- 1;} else {element <- 2;}
		
		if (stem.fit[[element]] < node.fit[[element]])
		{
			fit1 <- fit1.stem;
			fit2 <- fit2.stem;
			cut.at <- "stem";
#			cat("Stem wins!\n");
		} else {
			fit1 <- fit1.node;
			fit2 <- fit2.node;
			cut.at <- "node";
#			cat("Node wins!\n");
		}
	}
	op[aff[1],] <- fit1$par; # Replace parameters with new values for diminished clade
	
	fit$par <- rbind(op, fit2$par);
	fit$lnLik.part[aff] <- c(fit1$lnLik, fit2$lnLik); # Replace parameters with new values for diminished clade
	fit$split.at <- c(fit$split.at, node);
	fit$lnLik <- sum(fit$lnLik.part);
	
	model.fit <- calculate.model.fit(fit=fit, z=z);
	
	fit$aic <- model.fit[1];
	fit$aicc <- model.fit[2];
	fit$num.par <- model.fit[3];
	
	fit$cut.at <- c(fit$cut.at, cut.at);
	fit$model <- c(fit$model, fit2$model);
	
	return(fit);
}

