runTurboMEDUSA <-
function(phy, richness=NULL, model.limit=20, stop="model.limit", model="mixed",
	criterion="aicc", shiftCut="both", initialR=0.05, initialE=0.5, plotFig=FALSE, nexus=FALSE,
	verbose=TRUE, mc=FALSE, num.cores=NULL, ...)
{
	if (nexus) phy <- read.nexus(phy);
	if (is.null(richness))  # Assume tree represents single species tips and is completely sampled
	{
		richness <- data.frame(taxon=phy$tip.label, n.taxa=1);
	} else {
## Before determining model.limit, prune tree as necessary (from 'taxon' information in 'richness')
		phyData <- prune.tree.merge.data(phy=phy, richness=richness, verbose=verbose);
		phy <- phyData$phy;
		richness <- phyData$richness;
	}
	
## Limit on number of piecewise models fitted; based on tree size, aicc correction factor, 
## and flavour of model fitted (i.e. # parameters estimated; birth-death or pure-birth)
	model.limit <- get.max.model.limit(richness=richness, model.limit=model.limit, model=model, stop=stop, verbose=verbose);
	
## Determine correct AICc threshold from tree size (based on simulations)
 ## Should be used for interpreting model-fit
	threshold <- get.threshold(length(phy$tip.label));
	cat("Appropriate AICc threshold for tree of ", length(phy$tip.label), " tips is: ", threshold, ".\n\n", sep="");
	
## Store pertinent information: branch times, richness, descendants
	cat("Preparing data for analysis... ");
	
## Keep track of all nodes, internal and pendant (for keeping track of breakpoints)
	pend.nodes <- seq_len(length(phy$tip.label));   # Calculate pendant splits just once, keep track through various models
	int.nodes <- (length(phy$tip.label)+2):max(phy$edge); # Omit root node
	root.node <- length(phy$tip.label) + 1;
	all.nodes <- c(pend.nodes, root.node, int.nodes);
	
## The important bits. Set up z, get descendants and number of tips per node
	obj <- make.cache.medusa(phy=phy, richness=richness, all.nodes=all.nodes, mc=mc, num.cores=num.cores);
	
	desc <- list(desc.stem=obj$desc.stem, desc.node=obj$desc.node);
	z <- obj$z;
	z.orig <- z; # Save for summarizing models later
	num.tips <- obj$num.tips;
	
	cat("done.\n");
	
## Pre-fit pendant edges so these values need not be re(re(re))calculated; amounts to ~25% of all calculations
 ## Will show particular performance gain for edges with many fossil observations
	cat("Optimizing parameters for pendant edges... ");
	tips <- NULL;
# Will always be shiftCut="stem"; if mixed model, keep only best fit and throw out other in medusa.ml.prefit
	if (mc)
	{
		tips <- mclapply(pend.nodes, medusa.ml.prefit, z=z, desc=desc, initialR=initialR, initialE=initialE, model=model, shiftCut="stem", criterion=criterion, mc.cores=num.cores);
	} else {
		tips <- lapply(pend.nodes, medusa.ml.prefit, z=z, desc=desc, initialR=initialR, initialE=initialE, model=model, shiftCut="stem", criterion=criterion);
	}
	cat("done.\n");
	
## Pre-fit virgin internal nodes; should deliver performance gain for early models, and especially for large trees
 ## Remain useful until a spilt is accepted within the clade
 ## Need to incorporate cutAtStem here
	cat("Pre-calculating parameters for virgin internal nodes... ");
	virgin.stem <- list(); virgin.node <- list();
	if (mc)
	{
		if (shiftCut == "stem" || shiftCut == "both")
		{
			virgin.stem <- mclapply(int.nodes, medusa.ml.prefit, z=z, desc=desc, initialR=initialR, initialE=initialE, model=model, shiftCut="stem", criterion=criterion, mc.cores=num.cores);
		}
		if (shiftCut == "node" || shiftCut == "both")
		{
			virgin.node <- mclapply(int.nodes, medusa.ml.prefit, z=z, desc=desc, initialR=initialR, initialE=initialE, model=model, shiftCut="node", criterion=criterion, mc.cores=num.cores);
		}
	} else {
		if (shiftCut == "stem" || shiftCut == "both")
		{
			virgin.stem <- lapply(int.nodes, medusa.ml.prefit, z=z, desc=desc, initialR=initialR, initialE=initialE, model=model, shiftCut="stem", criterion=criterion);
		}
		if (shiftCut == "node" || shiftCut == "both")
		{
			virgin.node <- lapply(int.nodes, medusa.ml.prefit, z=z, desc=desc, initialR=initialR, initialE=initialE, model=model, shiftCut="node", criterion=criterion);
		}
	}
	virgin.nodes <- list(stem=virgin.stem, node=virgin.node);
	cat("done.\n\n");
	
	prefit <- list(tips=tips, virgin.nodes=virgin.nodes);
	
## Fit the base model
## 'fit' holds current results; useful for initializing subsequent models
	fit <- list();
	if (model == "mixed")
	{
		fit.bd <- medusa.ml.initial(z=z, initialR=initialR, initialE=initialE, model="bd");
		fit.yule <- medusa.ml.initial(z=z, initialR=initialR, initialE=initialE, model="yule");
		if (fit.bd[[criterion]] < fit.yule[[criterion]]) {
			fit <- fit.bd;
			fit$model <- "bd";
		} else {
			fit <- fit.yule;
			fit$model <- "yule";
		}
	} else if (model == "bd") {
		fit <- medusa.ml.initial(z=z, initialR=initialR, initialE=initialE, model="bd");
		fit$model <- "bd";
	} else {
		fit <- medusa.ml.initial(z=z, initialR=initialR, initialE=initialE, model="yule");
		fit$model <- "yule";
	}
	models <- list(fit);
	
	if (stop == "model.limit")
	{
		cat("Step 1 (of ", model.limit, "): best likelihood = ", models[[1]]$lnLik, "; AICc = ", models[[1]]$aicc, "; model = ", models[[1]]$model, "\n", sep="");
		for (i in seq_len(model.limit-1))
		{
			node.list <- all.nodes[-fit$split.at]; # this seems inefficient; create list and delete from there; will need -fit$split.at[i+1]
			
			
			
			if (mc)  # multicore (i.e. multithreaded) processing. No GUI, and not at all on Windows
			{
				res <- mclapply(node.list, medusa.ml.update, z=z, desc=desc, fit=fit, prefit=prefit, num.tips=num.tips, root.node=root.node, model=model, criterion=criterion, shiftCut=shiftCut, mc.cores=num.cores);
			} else {
				res <- lapply(node.list, medusa.ml.update, z=z, desc=desc, fit=fit, prefit=prefit, num.tips=num.tips, root.node=root.node, model=model, criterion=criterion, shiftCut=shiftCut);
			}
# Select model with best score according to the specific criterion employed (default aicc)
			best <- which.min(unlist(lapply(res, "[[", criterion)));
			models <- c(models, res[best]);
			fit <- res[[best]];   # keep track of '$split.at' i.e. nodes already considered
			
			z <- medusa.split(node=node.list[best], z=z, desc=desc, shiftCut=fit$cut.at[i+1])$z;
			
			cat("Step ", i+1, " (of ", model.limit, "): best likelihood = ", round(models[[i+1]]$lnLik, digits=7), "; AICc = ", models[[i+1]]$aicc,
				"; shift at node ", models[[i+1]]$split.at[i+1], "; model=", models[[i+1]]$model[i+1], "; cut=", models[[i+1]]$cut.at[i+1], "\n", sep="");
		}
	} else if (stop == "threshold") {
		i <- 1;
		done <- FALSE;
		cat("Step 1: best likelihood = ", round(models[[i+1]]$lnLik, digits=7), "; AICc = ", round(models[[i+1]]$aicc, digits=7), "; model = ", models[[1]]$model[i+1], "\n", sep="");
		while (!done & i < model.limit)
		{
			node.list <- all.nodes[-fit$split.at]; # this seems inefficient; create list and delete from there; will need -fit$split.at[i+1]
			
			
			
			if (mc)  # multicore (i.e. multithreaded) processing. No GUI, and not at all on Windows
			{
				res <- mclapply(node.list, medusa.ml.update, z=z, desc=desc, fit=fit, prefit=prefit, num.tips=num.tips, root.node=root.node, model=model, criterion=criterion, shiftCut=shiftCut, mc.cores=num.cores);
			} else {
				res <- lapply(node.list, medusa.ml.update, z=z, desc=desc, fit=fit, prefit=prefit, num.tips=num.tips, root.node=root.node, model=model, criterion=criterion, shiftCut=shiftCut);
			}
# Select model with best score according to the specific criterion employed (default aicc)
			best <- which.min(unlist(lapply(res, "[[", criterion)));
	# Compare last accepted model to current best model
			if (as.numeric(models[[length(models)]][criterion]) - as.numeric(res[[best]][criterion]) < threshold)
			{
				cat("\nNo significant increase in ", criterion, " score. Disregarding subsequent piecewise models.\n", sep="");
				done <- TRUE;
				break;
			}
			models <- c(models, res[best]);
			fit <- res[[best]];   # keep track of '$split.at' i.e. nodes already considered
			
			z <- medusa.split(node=node.list[best], z=z, desc=desc, shiftCut=fit$cut.at[i+1])$z;
			
			cat("Step ", i+1, ": best likelihood = ", round(models[[i+1]]$lnLik, digits=7), "; AICc = ", round(models[[i+1]]$aicc, digits=7),
				"; shift at node ", models[[i+1]]$split.at[i+1], "; model=", models[[i+1]]$model[i+1], "; cut=", models[[i+1]]$cut.at[i+1], "\n", sep="");
			i <- i+1;
		}
	}
	
	modelSummary <- calculate.model.fit.summary(models=models, phy=phy, plotFig=ifelse(length(models) > 1 & plotFig & !mc, TRUE, FALSE), threshold=threshold);
	if (verbose)
	{
		cat("\n", "Model fit summary:", "\n\n", sep="");
		print(modelSummary);
		if (threshold > 0)
		{
			cat("\nAIC weights are not reported, as they are meaningless when using a threshold criterion.\n")
		}
	}
	results <- list(z=z.orig, desc=desc, models=models, phy=phy, threshold=threshold, modelSummary=modelSummary);
	
#	return(results);
}

