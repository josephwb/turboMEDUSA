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


