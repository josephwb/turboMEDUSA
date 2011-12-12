runTurboMEDUSA <-
function(phy, richness=NULL, model.limit=20, stop="model.limit", model="bd",
	criterion="aicc", cutAtStem=TRUE, initialR=0.05, initialE=0.5, plotFig=FALSE, nexus=FALSE,
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
	obj <- make.cache.medusa(phy=phy, richness=richness, mc=mc, num.cores=num.cores);
	cat("done.\n");
	
## Keep track of all nodes, internal and pendant (for keeping track of breakpoints)
	pend.nodes <- seq_len(length(phy$tip.label));   # Calculate pendant splits just once, keep track through various models
	int.nodes <- (length(phy$tip.label)+2):max(phy$edge); # Omit root node
	root.node <- length(phy$tip.label) + 1;
	all.nodes <- c(pend.nodes, root.node, int.nodes);
	
	desc <- NULL;
	desc <- list(desc.stem=obj$desc.stem, desc.node=obj$desc.node)
	
	z <- obj$z;
	z.orig <- z; # Save for summarizing models later
	
## Pre-fit pendant edges so these values need not be re(re(re))calculated; amounts to ~25% of all calculations
 ## Will show particular performance gain for edges with many fossil observations
	cat("Optimizing parameters for pendant edges... ");
	tips <- NULL;
# Will always be "cutAtStem=TRUE"; if mixed model, keep only best fit and throw out other in medusa.ml.prefit
	if (mc)
	{
		tips <- mclapply(pend.nodes, medusa.ml.prefit, z=z, desc=desc, initialR=initialR, initialE=initialE, model=model, cutAtStem=TRUE, criterion=criterion, mc.cores=num.cores);
	} else {
		tips <- lapply(pend.nodes, medusa.ml.prefit, z=z, desc=desc, initialR=initialR, initialE=initialE, model=model, cutAtStem=TRUE, criterion=criterion);
	}
	cat("done.\n");
	
## Pre-fit virgin internal nodes; should deliver performance gain for early models, and especially for large trees
 ## Remain useful until a spilt is accepted within the clade
 ## Need to incorporate cutAtStem here
	cat("Pre-calculating parameters for virgin internal nodes... ");
	virgin.stem <- list(); virgin.node <- list();
	if (mc)
	{
		if (model == "bd" | model == "mixed")
		{
			virgin.nodes.bd <- mclapply(int.nodes, medusa.ml.prefit, z=z, desc=desc, initialR=initialR, initialE=initialE, model="bd", cutAtStem=cutAtStem, criterion=criterion, mc.cores=num.cores);
		}
		if (model == "yule" | model == "mixed")
		{
			virgin.nodes.yule <- mclapply(int.nodes, medusa.ml.prefit, z=z, desc=desc, initialR=initialR, initialE=initialE, model="yule", cutAtStem=cutAtStem, criterion=criterion, mc.cores=num.cores);
		}
	} else {
		if (cutAtStem == TRUE || cutAtStem == "both")
		{
			virgin.stem <- lapply(int.nodes, medusa.ml.prefit, z=z, desc=desc, initialR=initialR, initialE=initialE, model=model, cutAtStem=TRUE, criterion=criterion);
		}
		if (cutAtStem == FALSE || cutAtStem == "both")
		{
			virgin.node <- lapply(int.nodes, medusa.ml.prefit, z=z, desc=desc, initialR=initialR, initialE=initialE, model=model, cutAtStem=FALSE, criterion=criterion);
		}
	}
	virgin.nodes <- list(stem=virgin.stem, node=virgin.node);
	cat("done.\n\n");
	
	prefit <- list(tips=tips, virgin.nodes=virgin.nodes);
	
## Needed downstream; don't recalculate
 ## Gives the number of tips associated with an internal node; determines whether a node is 'virgin' or not
	num.tips <- list()
	if (mc)
	{
		num.tips <- mclapply(all.nodes, get.num.tips, phy=phy, mc.cores=num.cores);
	} else {
		num.tips <- lapply(all.nodes, get.num.tips, phy=phy);
	}
	
## Fit the base model
## 'fit' holds current results; useful for initializing subsequent models
	fit <- list();
	if (model == "mixed")
	{
		fit.bd <- medusa.ml.initial(z=z, initialR=initialR, initialE=initialE, model="bd");
		fit.yule <- medusa.ml.initial(z=z, initialR=initialR, initialE=initialE, model="yule");
		if (fit.bd[[criterion]] < fit.yule[[criterion]]) {
			fit <- fit.bd;
		} else {
			fit <- fit.yule;
		}
	} else if (model == "bd") {
		fit <- medusa.ml.initial(z=z, initialR=initialR, initialE=initialE, model="bd");
	} else if (model == "yule") {
		fit <- medusa.ml.initial(z=z, initialR=initialR, initialE=initialE, model="yule");
	}
	models <- list(fit);
	
	if (stop == "model.limit")
	{
		cat("Step 1 (of ", model.limit, "): best likelihood = ", models[[1]]$lnLik, "; AICc = ", models[[1]]$aicc, "\n", sep="");
		for (i in seq_len(model.limit-1))
		{
			node.list <- all.nodes[-fit$split.at];
			if (mc)  # multicore (i.e. multithreaded) processing. No GUI, and not at all on Windows
			{
				res <- mclapply(node.list, medusa.ml.update, z=z, desc=desc, fit=fit, prefit=prefit, num.tips=num.tips, root.node=root.node, model=model, criterion=criterion, cutAtStem=cutAtStem, mc.cores=num.cores);
			} else {
				res <- lapply(node.list, medusa.ml.update, z=z, desc=desc, fit=fit, prefit=prefit, num.tips=num.tips, root.node=root.node, model=model, criterion=criterion, cutAtStem=cutAtStem);
			}
# Select model with best score according to the specific criterion employed (default aicc)
			best <- which.min(unlist(lapply(res, "[[", criterion)));
			models <- c(models, res[best]);
			fit <- res[[best]];   # keep track of '$split.at' i.e. nodes already considered
			
			z <- medusa.split(node=node.list[best], z=z, desc=desc, stemCut=ifelse(fit$cut.at == "stem", TRUE, FALSE))$z;
			
			cat("Step ", i+1, " (of ", model.limit, "): best likelihood = ", models[[i+1]]$lnLik, "; AICc = ", models[[i+1]]$aicc,
				"; break at node ", models[[i+1]]$split.at[i+1], "; cut=", models[[i+1]]$cut.at, "\n", sep="");
		}
	} else if (stop == "threshold") {
		i <- 1;
		done <- FALSE;
		cat("Step 1: best likelihood = ", models[[1]]$lnLik, "; AICc = ", models[[1]]$aicc, "\n", sep="");
		while (!done & i < model.limit)
		{
			node.list <- all.nodes[-fit$split.at];
			if (mc)  # multicore (i.e. multithreaded) processing. No GUI, and not at all on Windows
			{
				res <- mclapply(node.list, medusa.ml.update, z=z, desc=desc, fit=fit, prefit=prefit, num.tips=num.tips, root.node=root.node, model=model, criterion=criterion, cutAtStem=cutAtStem, mc.cores=num.cores);
			} else {
				res <- lapply(node.list, medusa.ml.update, z=z, desc=desc, fit=fit, prefit=prefit, num.tips=num.tips, root.node=root.node, model=model, criterion=criterion, cutAtStem=cutAtStem);
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
			
			
			z <- medusa.split(node=node.list[best], z=z, desc=desc, stemCut=ifelse(fit$cut.at == "stem", TRUE, FALSE))$z;
			
			
			cat("Step ", i+1, ": best likelihood = ", models[[i+1]]$lnLik, "; AICc = ", models[[i+1]]$aicc,
				"; break at node ", models[[i+1]]$split.at[i+1], "; cut=", models[[i+1]]$cut.at, "\n", sep="");
			i <- i+1;
		}
	}
	
	modelSummary <- calculate.model.fit.summary(models=models, phy=phy, plotFig=ifelse(length(models) > 1 & plotFig & !mc, TRUE, FALSE));
	if (verbose)
	{
		cat("\n", "Model fit summary:", "\n\n", sep="");
		print(modelSummary);
	}
	results <- list(z=z.orig, desc=desc, models=models, phy=phy, threshold=threshold, modelSummary=modelSummary);
	
	return(results);
}

