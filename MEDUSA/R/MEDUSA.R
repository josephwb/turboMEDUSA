MEDUSA <- function(phy, richness=NULL, model="mixed", modelLimit=20, stop="threshold",
	shiftCut="both", criterion="aicc", stepBack=TRUE, preserveModelFlavour=FALSE, epsilon=NULL, r=NULL,
	b=NULL, d=NULL, fixThreshold=NULL, initialR=0.05, initialE=0.5, verbose=TRUE, mc=FALSE, numCores=NULL,
	resolveTree=FALSE, ...)
{
	checkValidArguments(phy, richness, model, modelLimit, stop, shiftCut, criterion, stepBack,
		preserveModelFlavour, epsilon, r, b, d, fixThreshold, initialR, initialE,
		verbose, mc, numCores, resolveTree);
	
	conf <- configureModel(model=model, epsilon=epsilon, r=r, b=b, d=d, initialR=initialR, initialE=initialE);
	sp <- conf$sp;
	model <- conf$model;
	fixPar <- conf$fixPar;
	richness <- formatRichness(richness, phy); # do this once, instead of once per tree
	
	# Determine correct AICc threshold from tree size (based on simulations)
	 # Should be used for interpreting model-fit
	threshold <- getThreshold(phy, fixThreshold=fixThreshold, stop=stop);
		
	# Limit on number of piecewise models fitted; based on tree size, aicc correction factor, 
	# and flavour of model fitted (i.e. # parameters estimated; birth-death or pure-birth)
	modelLimit <- getMaxModelLimit(richness=richness, modelLimit=modelLimit, model=model, stop=stop);
	
	runMEDUSA <- function (phy, richness, multiTree=FALSE, verbose, ...) { # wtf is multiTree?!?
		
		phyData <- prepareData(phy=phy, richness=richness, verbose=verbose, resolveTree);
		phy <- phyData$phy;
		richness <- phyData$richness;
		
	# Store pertinent information: branch times, richness, descendants
		cat("Preparing data for analysis:\n");
	# Keep track of all nodes, internal and pendant (for keeping track of breakpoints)
		pend.nodes <- seq_len(length(phy$tip.label));   # Calculate pendant splits just once, keep track through various models
		int.nodes <- unique(phy$edge[,1])[-1]; # Omit root node
		root.node <- length(phy$tip.label) + 1;
		all.nodes <- c(pend.nodes, root.node, int.nodes);
		
	# The important bits. Set up z, get descendants and number of tips per node
		obj <- makeCacheMedusa(phy=phy, richness=richness, all.nodes=all.nodes, shiftCut=shiftCut, mc=mc, numCores=numCores);
		desc <- list(stem=obj$desc.stem, node=obj$desc.node);
		z <- obj$z;
		num.tips <- obj$num.tips;
		
		cat("done.\n\n");
		
	# Fit the base model. This is done first as it will catch potential invalid fixed parameter values.
	# 'fit' holds current results; useful for initializing subsequent models
		baseFit <- list();
		baseFit <- medusaMLFitBase(z=z, sp=sp, model=model, fixPar=fixPar, criterion=criterion);
		if (baseFit$lnLik == -Inf && !is.null(fixPar))
		{stop("\n\nConstrained model cannot be fit to data with current fixed parameter value. Stopping analysis.\n\n");}
		
	# If only one model is desired (i.e. base model), don't bother with all of the precalculations.
		if (modelLimit != 1) {
	# Pre-fit pendant edges so these values need not be re(re(re))calculated; amounts to ~25% of all calculations
	# Will show particular performance gain for edges with many fossil observations
			cat("Optimizing parameters for pendant edges... ");
			tips <- NULL;
	# Will always be shiftCut="stem"; if mixed model, keep only best fit and throw out other in medusaMLPrefit
			tips <- prefitTips(pend.nodes=pend.nodes, z=z, sp=sp, model=model, fixPar=fixPar, criterion=criterion,
				mc=mc, numCores=numCores);
			cat("done.\n");
			
	# Pre-fit virgin internal nodes; should deliver performance gain for early models, and especially for large trees
	 # Remain useful until a spilt is accepted within the clade
			if (length(int.nodes) > 0) { # hmm, when would this be false?
				cat("Pre-calculating parameters for internal nodes... ");
				virgin.stem <- list(); virgin.node <- list();
				if (mc) {
					if (shiftCut == "stem" || shiftCut == "both") {
						virgin.stem <- mclapply(int.nodes, medusaMLPrefitStem, z=z, desc=desc$stem, sp=sp, model=model,
							fixPar=fixPar, criterion=criterion, mc.cores=numCores);
					}
					if (shiftCut == "node" || shiftCut == "both") {
						virgin.node <- mclapply(int.nodes, medusaMLPrefitNode, z=z, desc=desc$node, sp=sp, model=model,
							fixPar=fixPar, criterion=criterion, mc.cores=numCores);
					}
				} else {
					if (shiftCut == "stem" || shiftCut == "both") {
						virgin.stem <- lapply(int.nodes, medusaMLPrefitStem, z=z, desc=desc$stem, sp=sp, model=model,
							fixPar=fixPar, criterion=criterion);
					}
					if (shiftCut == "node" || shiftCut == "both") {
						virgin.node <- lapply(int.nodes, medusaMLPrefitNode, z=z, desc=desc$node, sp=sp, model=model,
							fixPar=fixPar, criterion=criterion);
					}
				}
				virgin.nodes <- list(stem=virgin.stem, node=virgin.node);
			} else {
				virgin.nodes <- NULL;
			}
			
			cat("done.\n\n");
			
			prefit <- list(tips=tips, virgin.nodes=virgin.nodes, num.tips=num.tips);
		}
		
		stopOnLimit <- function(fit) { # not likely to be used. get rid of it?
			optModel <- fit;
			
			cat("Step 1 (of ", modelLimit, "): lnLik=", optModel$lnLik, "; AICc=", optModel$aicc,
				"; model=", optModel$model, "\n", sep="");
			
			for (i in seq_len(modelLimit-1)) {
				node.list <- all.nodes[-fit$split.at];
				
				if (mc) { # multicore (i.e. multithreaded) processing. No GUI, and not at all on Windows
					
					res <- mclapply(node.list, medusaMLUpdate, z=z, desc=desc, fit=fit, prefit=prefit, root.node=root.node,
						model=model, fixPar=fixPar, criterion=criterion, shiftCut=shiftCut, preserveModelFlavour=preserveModelFlavour, mc.cores=numCores);
						
				} else {
					res <- lapply(node.list, medusaMLUpdate, z=z, desc=desc, fit=fit, prefit=prefit, root.node=root.node,
						model=model, fixPar=fixPar, criterion=criterion, shiftCut=shiftCut, preserveModelFlavour=preserveModelFlavour);
				}
	# Select model with best score according to the specific criterion employed (default aicc)
				best <- res[[which.min(unlist(lapply(res, "[[", criterion)))]];
				
				node <- best$split.at;
				cut <- best$cut.at;
				
			# re-optimize on best break node/position
				fit <- medusaFitOptimal(node=node, z=z, desc=desc, fit=fit, prefit=prefit, root.node=root.node, model=model,
					fixPar=fixPar, criterion=criterion, shiftCut=cut, preserveModelFlavour=preserveModelFlavour);
				
				z <- medusaSplit(node=node, z=z, desc=desc, shiftCut=tail(fit$cut.at,1))$z;
				step <- rbind(optModel$step, c("add", tail(fit$split.at,1)));
				
	# Consider parameter removal
				if (stepBack) {
					backFit <- backStep(currentModel=fit, z=z, step=step, model=model, fixPar=fixPar, criterion=criterion);
					fit <- backFit$fit;
					z <- backFit$z;
					step <- backFit$step;
				}
				
				fit$z <- z;
				fit$step <- step;
				optModel <- fit;
				
				cat("Step ", i+1, " (of ", modelLimit, "): lnLik=", round(fit$lnLik, digits=7),
					"; AICc=", fit$aicc, "; shift at node ", tail(fit$split.at,1), "; model=",
					tail(fit$model,1), "; cut=", tail(fit$cut.at,1), "; # shifts=", length(fit$split.at) - 1, "\n", sep="");
			}
			return(optModel);
		}
		
		stopOnThreshold <- function (fit) {
			optModel <- fit;
			i <- 1;
			done <- FALSE;
			
			cat("Step 1: lnLik=", round(optModel$lnLik, digits=7), "; AICc=", round(optModel$aicc, digits=7),
				"; model=", optModel$model[1], "\n", sep="");
			
			while (!done && i) {
				node.list <- all.nodes[-fit$split.at];
				if (mc)  # multicore (i.e. multithreaded) processing. No GUI, and not at all on Windows
				{
					res <- mclapply(node.list, medusaMLUpdate, z=z, desc=desc, fit=fit, prefit=prefit, root.node=root.node,
						model=model, fixPar=fixPar, criterion=criterion, shiftCut=shiftCut, preserveModelFlavour=preserveModelFlavour,
						mc.cores=numCores);
				} else {
					res <- lapply(node.list, medusaMLUpdate, z=z, desc=desc, fit=fit, prefit=prefit, root.node=root.node,
						model=model, fixPar=fixPar, criterion=criterion, shiftCut=shiftCut, preserveModelFlavour=preserveModelFlavour);
				}
	# Select model with best score according to the specific criterion employed (default aicc)
				best <- res[[which.min(unlist(lapply(res, "[[", criterion)))]];
				
				node <- best$split.at;
				cut <- best$cut.at;
				
			# re-optimize on best break node/position
				fit <- medusaFitOptimal(node=node, z=z, desc=desc, fit=fit, prefit=prefit, root.node=root.node, model=model,
					fixPar=fixPar, criterion=criterion, shiftCut=cut, preserveModelFlavour=preserveModelFlavour);
				
				z <- medusaSplit(node=node, z=z, desc=desc, shiftCut=tail(fit$cut.at,1))$z;
				step <- rbind(optModel$step, c("add", tail(fit$split.at,1)));
				
	# Consider parameter removal
				if (stepBack) {
					backFit <- backStep(currentModel=fit, z=z, step=step, model=model, fixPar=fixPar, criterion=criterion);
					fit <- backFit$fit;
					z <- backFit$z;
					step <- backFit$step;
				}
				
				fit$z <- z;
				fit$step <- step;
				
				if (as.numeric(optModel[criterion]) - as.numeric(fit[criterion]) < threshold) {
					if (verbose) cat("\nNo significant increase in ", criterion, " score. Disregarding subsequent piecewise models.\n", sep="");
					done <- TRUE;
					break;
				}
				
				if (!is.null(backFit$remove)) {printRemovedShifts(remove=backFit$remove);}
				optModel <- fit;
				
				cat("Step ", i+1, ": lnLik=", round(fit$lnLik, digits=7),
					"; AICc=", fit$aicc, "; shift at node ", tail(fit$split.at,1), "; model=",
					tail(fit$model,1), "; cut=", tail(fit$cut.at,1), "; # shifts=", length(fit$split.at) - 1, "\n", sep="");
				i <- i+1;
			}
			return(optModel);
		}
		
		if (stop == "modelLimit") optModel <- stopOnLimit(fit=baseFit) else optModel <- stopOnThreshold(fit=baseFit);
		
		modelSummary <- optModelSummary(optModel);
		
		results <- list(desc=desc, optModel=optModel, phy=phy, fixPar=fixPar, criterion=criterion,
			stop=stop, threshold=threshold, modelSummary=modelSummary);
		class(results) <- "medusa";
		
		if (verbose) {cat("\n"); cat("Summary of optimal MEDUSA model:\n"); print(modelSummary)};
		cat("\n");
		return(results);
	}
	
	if (class(phy) == "multiPhylo") {
		phy <- manageTipLabels(phy=phy, mc=mc, numCores=numCores);
		results <- lapply(phy, runMEDUSA, richness=richness, multiTree=TRUE, verbose=FALSE, ...); # prevent extraneous bits from being printed to screen
		results <- list(results=results, richness=richness, medusaVersion=packageVersion("MEDUSA"));
		class(results) <- "multiMedusa";
	} else {
		results <- runMEDUSA(phy=phy, richness=richness, verbose=verbose, ...);
		results$richness <- richness;
		results$medusaVersion <- packageVersion("MEDUSA");
	}
	
	invisible(results);
}