runTurboMEDUSA <- function(phy=phy, richness=NULL, modelLimit=20, stop="modelLimit", model="mixed",
	fixedEpsilon=NULL, criterion="aicc", shiftCut="both", initialR=0.05, initialE=0.5,
	preserveModelFlavour=TRUE, plotFig=FALSE, verbose=TRUE, mc=FALSE, numCores=NULL, ...)
{
	phyData <- prepareData(phy=phy, richness=richness, verbose=verbose);
	phy <- phyData$phy;
	richness <- phyData$richness;
	sp <- c(initialR, initialE);
		
## Determine correct AICc threshold from tree size (based on simulations)
 ## Should be used for interpreting model-fit
	threshold <- getThreshold(length(phy$tip.label));
	cat("Appropriate AICc threshold for tree of ", length(phy$tip.label), " tips is: ", threshold, ".\n", sep="");
	
## Limit on number of piecewise models fitted; based on tree size, aicc correction factor, 
## and flavour of model fitted (i.e. # parameters estimated; birth-death or pure-birth)
	modelLimit <- getMaxModelLimit(richness=richness, modelLimit=modelLimit, model=model, stop=stop, verbose=verbose);
	
## Store pertinent information: branch times, richness, descendants
	cat("Preparing data for analysis... ");
## Keep track of all nodes, internal and pendant (for keeping track of breakpoints)
	pend.nodes <- seq_len(length(phy$tip.label));   # Calculate pendant splits just once, keep track through various models
	int.nodes <- (length(phy$tip.label)+2):max(phy$edge); # Omit root node
	root.node <- length(phy$tip.label) + 1;
	all.nodes <- c(pend.nodes, root.node, int.nodes);
	
## The important bits. Set up z, get descendants and number of tips per node
	obj <- makeCacheMedusa(phy=phy, richness=richness, all.nodes=all.nodes, mc=mc, numCores=numCores);
	
	desc <- list(desc.stem=obj$desc.stem, desc.node=obj$desc.node);
	z <- obj$z;
	z.orig <- z; # Save for summarizing models later
	num.tips <- obj$num.tips;
	
	cat("done.\n");
		
## Pre-fit pendant edges so these values need not be re(re(re))calculated; amounts to ~25% of all calculations
## Will show particular performance gain for edges with many fossil observations
	cat("Optimizing parameters for pendant edges... ");
	tips <- NULL;
# Will always be shiftCut="stem"; if mixed model, keep only best fit and throw out other in medusaMLPrefit
	if (mc) {
		tips <- mclapply(pend.nodes, medusaMLPrefitStem, z=z, desc=desc$desc.stem, sp=sp,
			model=model, criterion=criterion, mc.cores=numCores);
	} else {
		tips <- lapply(pend.nodes, medusaMLPrefitStem, z=z, desc=desc$desc.stem, sp=sp,
			model=model, criterion=criterion);
	}
	cat("done.\n");
	
## Pre-fit virgin internal nodes; should deliver performance gain for early models, and especially for large trees
 ## Remain useful until a spilt is accepted within the clade
 ## Need to incorporate cutAtStem here
	cat("Pre-calculating parameters for virgin internal nodes... ");
	virgin.stem <- list(); virgin.node <- list();
	if (mc) {
		if (shiftCut == "stem" || shiftCut == "both") {
			virgin.stem <- mclapply(int.nodes, medusaMLPrefitStem, z=z, desc=desc$desc.stem, sp=sp, model=model,
				criterion=criterion, mc.cores=numCores);
		}
		if (shiftCut == "node" || shiftCut == "both") {
			virgin.node <- mclapply(int.nodes, medusaMLPrefitNode, z=z, desc=desc$desc.node, sp=sp, model=model,
				criterion=criterion, mc.cores=numCores);
		}
	} else {
		if (shiftCut == "stem" || shiftCut == "both") {
			virgin.stem <- lapply(int.nodes, medusaMLPrefitStem, z=z, desc=desc$desc.stem, sp=sp, model=model,
				criterion=criterion);
		}
		if (shiftCut == "node" || shiftCut == "both") {
			virgin.node <- lapply(int.nodes, medusaMLPrefitNode, z=z, desc=desc$desc.node, sp=sp, model=model,
				criterion=criterion);
		}
	}
	virgin.nodes <- list(stem=virgin.stem, node=virgin.node);
	cat("done.\n\n");
	
	prefit <- list(tips=tips, virgin.nodes=virgin.nodes, num.tips=num.tips);
	
## Fit the base model
## 'fit' holds current results; useful for initializing subsequent models
	fit <- list();
	fit <- medusaMLFitBase(z=z, sp=sp, model=model, criterion=criterion);
	models <- list(fit);

	if (stop == "modelLimit") {
		
		cat("Step 1 (of ", modelLimit, "): best likelihood = ", models[[1]]$lnLik, "; AICc = ", models[[1]]$aicc,
			"; model = ", models[[1]]$model, "\n", sep="");
		
		for (i in seq_len(modelLimit-1)) {
			node.list <- all.nodes[-fit$split.at]; # this seems inefficient; create list and delete from there; will need -fit$split.at[i+1]
			
			if (mc) { # multicore (i.e. multithreaded) processing. No GUI, and not at all on Windows
				
				res <- mclapply(node.list, medusaMLUpdate, z=z, desc=desc, fit=fit, prefit=prefit, root.node=root.node,
					model=model, criterion=criterion, shiftCut=shiftCut, preserveModelFlavour=preserveModelFlavour, mc.cores=numCores);
					
			} else {
				res <- lapply(node.list, medusaMLUpdate, z=z, desc=desc, fit=fit, prefit=prefit, root.node=root.node,
					model=model, criterion=criterion, shiftCut=shiftCut, preserveModelFlavour=preserveModelFlavour);
			}
# Select model with best score according to the specific criterion employed (default aicc)
			best <- which.min(unlist(lapply(res, "[[", criterion)));
			models <- c(models, res[best]);
			fit <- res[[best]];   # keep track of '$split.at' i.e. nodes already considered
			
			z <- medusaSplit(node=node.list[best], z=z, desc=desc, shiftCut=fit$cut.at[i+1])$z;
			
			cat("Step ", i+1, " (of ", modelLimit, "): best likelihood = ", round(models[[i+1]]$lnLik, digits=7),
				"; AICc = ", models[[i+1]]$aicc, "; shift at node ", models[[i+1]]$split.at[i+1], "; model=",
				models[[i+1]]$model[i+1], "; cut=", models[[i+1]]$cut.at[i+1], "\n", sep="");
		}
	} else if (stop == "threshold") {
		
		i <- 1;
		done <- FALSE;
		
		cat("Step 1: best likelihood = ", round(models[[1]]$lnLik, digits=7), "; AICc = ", round(models[[1]]$aicc, digits=7),
			"; model = ", models[[1]]$model[1], "\n", sep="");
		
		while (!done & i < modelLimit) {
			node.list <- all.nodes[-fit$split.at]; # this seems inefficient; create list and delete from there; will need -fit$split.at[i+1]
			if (mc)  # multicore (i.e. multithreaded) processing. No GUI, and not at all on Windows
			{
				res <- mclapply(node.list, medusaMLUpdate, z=z, desc=desc, fit=fit, prefit=prefit, root.node=root.node,
					model=model, criterion=criterion, shiftCut=shiftCut, preserveModelFlavour=preserveModelFlavour, mc.cores=numCores);
			} else {
				res <- lapply(node.list, medusaMLUpdate, z=z, desc=desc, fit=fit, prefit=prefit, root.node=root.node,
					model=model, criterion=criterion, shiftCut=shiftCut, preserveModelFlavour=preserveModelFlavour);
			}
# Select model with best score according to the specific criterion employed (default aicc)
			best <- which.min(unlist(lapply(res, "[[", criterion)));
	# Compare last accepted model to current best model
			if (as.numeric(models[[length(models)]][criterion]) - as.numeric(res[[best]][criterion]) < threshold) {
				cat("\nNo significant increase in ", criterion, " score. Disregarding subsequent piecewise models.\n", sep="");
				done <- TRUE;
				break;
			}
			models <- c(models, res[best]);
			fit <- res[[best]];   # keep track of '$split.at' i.e. nodes already considered
			
			z <- medusaSplit(node=node.list[best], z=z, desc=desc, shiftCut=fit$cut.at[i+1])$z;
			
			cat("Step ", i+1, ": best likelihood = ", round(models[[i+1]]$lnLik, digits=7), "; AICc = ",
				round(models[[i+1]]$aicc, digits=7), "; shift at node ", models[[i+1]]$split.at[i+1],
				"; model=", models[[i+1]]$model[i+1], "; cut=", models[[i+1]]$cut.at[i+1], "\n", sep="");
			i <- i+1;
		}
	}
	
	modelSummary <- calculateModelFitSummary(models=models, phy=phy, plotFig=ifelse(length(models) > 1 & plotFig & !mc, TRUE, FALSE),
		threshold=threshold);
	if (verbose) {
		cat("\n", "Model fit summary:", "\n\n", sep="");
		print(modelSummary);
		# if (threshold > 0)
		# {
			# cat("\nAIC weights are not reported, as they are meaningless when using a threshold criterion.\n")
		# }
	}
	results <- list(z=z.orig, desc=desc, models=models, phy=phy, threshold=threshold, modelSummary=modelSummary);
}