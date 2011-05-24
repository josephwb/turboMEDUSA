runTurboMEDUSA <-
function(phy, richness=NULL, model.limit=20, stop="model.limit",
	criterion="aicc", initial.r=0.05, initial.e=0.5, plotFig=FALSE, nexus=FALSE, verbose=TRUE, mc=FALSE, num.cores=NULL, ...)
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
## and flavour of model fitted (i.e. # parameters estimated; at the moment only birth-death considered)
	model.limit <- get.max.model.limit(richness, model.limit, stop, verbose)
	
## Determine correct AICc threshold from tree size (based on simulations)
 ## Should be used for interpreting model-fit
	threshold <- get.threshold(length(phy$tip.label))
	cat("Appropriate AICc threshold for tree of ", length(phy$tip.label), " tips is: ", threshold, ".\n\n", sep="")
	
## Store pertinent information: branch times, richness, ancestors
	cat("Preparing data for analysis... ")
	obj <- make.cache.medusa(phy, richness, mc, num.cores)
	cat("done.\n")
	
## Keep track of all nodes, internal and pendant (for keeping track of breakpoints)
	pend.nodes <- seq_len(length(phy$tip.label))   # Calculate pendant splits just once, keep track through various models
	int.nodes <- (length(phy$tip.label)+2):max(phy$edge) # Omit root node
	root.node <- length(phy$tip.label)+1
	all.nodes <- c(pend.nodes, root.node, int.nodes)
	
	anc <- obj$anc
	z <- obj$z
	z.orig <- z # Save for summarizing models
	
## Pre-fit pendant edges so these values need not be re(re(re))calculated; amounts to ~25% of all calculations
 ## Will show particular performance gain for edges with many fossil observations
	tips <- list()
	cat("Optimizing parameters for pendant edges... ")
	if (mc)
	{
		tips <- mclapply(pend.nodes, medusa.ml.prefit, z, anc, initial.r, initial.e, mc.cores=num.cores)
	} else {
		tips <- lapply(pend.nodes, medusa.ml.prefit, z, anc, initial.r, initial.e)
	}
	cat("done.\n")
	
## Pre-fit virgin internal nodes; should deliver performance gain for early models, and especially for large trees
 ## Remain useful until a spilt is accepted within the clade
	virgin.nodes <- list()
	cat("Pre-calculating parameters for virgin internal nodes... ")
	if (mc)
	{
		virgin.nodes <- mclapply(int.nodes, medusa.ml.prefit, z, anc, initial.r, initial.e, mc.cores=num.cores)
	} else {
		virgin.nodes <- lapply(int.nodes, medusa.ml.prefit, z, anc, initial.r, initial.e)
	}
	cat("done.\n\n")
	
## Needed downstream; don't recalculate
	num.tips <- list()
	if (mc)
	{
		num.tips <- mclapply(all.nodes, get.num.tips, phy, mc.cores=num.cores)
	} else {
		num.tips <- lapply(all.nodes, get.num.tips, phy)
	}
	
## 'fit' holds current results; useful for initializing subsequent models
	fit <- medusa.ml.initial(z, initial.r, initial.e)
	models <- list(fit)
	
	if (stop == "model.limit")
	{
		cat("Step 1 (of ", model.limit, "): best likelihood = ", models[[1]]$lnLik, "; AICc = ", models[[1]]$aicc, "\n", sep="")
		for (i in seq_len(model.limit-1))
		{
			node.list <- all.nodes[-fit$split.at]
			if (mc)  # multicore (i.e. multithreaded) processing. No GUI, and not at all on Windows
			{
				res <- mclapply(node.list, medusa.ml.update, z, anc, fit, tips, virgin.nodes, num.tips, root.node, mc.cores=num.cores)
			} else {
				res <- lapply(node.list, medusa.ml.update, z, anc, fit, tips, virgin.nodes, num.tips, root.node)
			}
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
			if (mc)  # multicore (i.e. multithreaded) processing. No GUI, and not at all on Windows
			{
				res <- mclapply(node.list, medusa.ml.update, z, anc, fit, tips, virgin.nodes, num.tips, root.node, mc.cores=num.cores)
			} else {
				res <- lapply(node.list, medusa.ml.update, z, anc, fit, tips, virgin.nodes, num.tips, root.node)
			}
# Select model with best score according to the specific criterion employed (default aicc)
			best <- which.min(unlist(lapply(res, "[[", criterion)))
			
			if (as.numeric(res[[best]][criterion]) - as.numeric(models[[length(models)]][criterion]) > threshold)
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
	
	model.summary <- calculate.model.fit.summary(models, phy, plotFig=ifelse(length(models) > 1 & plotFig & !mc, TRUE, FALSE))
	if (verbose)
	{
		cat("\n", "Model fit summary:", "\n\n", sep="")
		print(model.summary)
	}
	results <- list(z=z.orig, anc=anc, models=models, phy=phy, threshold=threshold, model.summary=model.summary)
	
	return(results)
}

