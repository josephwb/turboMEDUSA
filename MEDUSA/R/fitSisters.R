fitSisters <- function (phy, richness=NULL, node=NULL, tips=NULL, model="mixed",
	criterion="aicc", epsilon=NULL, r=NULL, b=NULL, d=NULL, plotSurface=FALSE,
	initialR=0.05, initialE=0.5, verbose=TRUE, mc=FALSE, numCores=NULL, ...) {
	if (is.null(tips) && is.null(node) && node != "root") stop("\n\nWarning: need to provide either a node number or set of 2 taxa whose MRCA is the desired node. Stopping.\n");
	
## Set up model	
	conf <- configureModel(model=model, epsilon=epsilon, r=r, b=b, d=d,
		initialR=initialR, initialE=initialE);
	sp <- conf$sp;
	model <- conf$model;
	fixPar <- conf$fixPar;
	
	runFitSisters <- function(phy, richness, ...) {
## Use richness information to prune tree. Could be several trees.
		phyData <- extractSubTree(phy=phy, tips=tips, node=node, richness=richness, verbose=verbose);
		phy <- phyData$phy;
		richness <- phyData$richness;
		
## Identify nodes
		pend.nodes <- seq_len(length(phy$tip.label));
		int.nodes <- (length(phy$tip.label)+2):max(phy$edge); # Omit root node
		root.node <- length(phy$tip.label) + 1;
		all.nodes <- c(pend.nodes, root.node, int.nodes);
		
## The important bits. Set up z, get descendants and number of tips per node
		obj <- makeCacheMedusa(phy=phy, richness=richness, all.nodes=all.nodes, shiftCut ="stem", mc=mc, numCores=numCores, verbose=verbose);
		desc <- obj$desc.stem;
		z <- obj$z;
		
## Get immediate descendant nodes of root. In case of multiple trees, make sure ordering is consistent.
		sisNodes <- z[z[,"anc"] == root.node,"dec"];
		
		fit <- medusaMLFitBase(z=z, sp=sp, model=model, fixPar=fixPar, criterion=criterion);
		models <- list(fit);
		sisterModel <- sisterFit(node=sisNodes[1], z=z, desc=desc, fit=fit, model=model, fixPar=fixPar, criterion=criterion);
		models <- c(models, sisterModel);
		
		return(list(models=models, phy=phy));
	}
	
	if (class(phy) == "phylo") {
		res <- runFitSisters(phy, richness);
	} else if (class(phy) == "multiPhylo") {
		if (mc)
		{
			res <- parallel::mclapply(phy, runFitSisters, richness);
		} else {
			res <- lapply(phy, runFitSisters, richness);
		}
	}
	
## Summary
	summ <- sisterFitSummary(results=res, criterion=criterion, fixPar=fixPar, plotSurface=plotSurface);
	
	results <- list(parameterSummary=summ$parSumm, modelSummary=summ$modSumm, models=res$models, fixPar=fixPar, nTrees=summ$nTrees);
	class(results) <- "sisterFit";
	
	invisible(results);
}

## 'fit' contains parameter values from previous model, used to initialize subsequent model.
## Pass in pre-fitted values for pendant edges and virgin nodes (in 'prefit'); DON'T recalculate.
## Need to consider the possibility of birth-death, yule, mixed, or constrained models.
## Need to consider where shft is placed (shiftCut). Placement affects not only new clade, but
## also the size of the split clade. Only relevant if shiftCut = "both".
sisterFit <- function (node, z, desc, fit, model, fixPar, criterion) {
	sp <- NULL;
	aff <- NULL;
	op <- fit$par; # store previously fit parameter values
	cut.at <- NULL;
	cool <- TRUE;
	
	fit1 <- NULL;
	fit2 <- NULL;
	
## First, diminshed clade
	obj <- medusaSplitStem(node=node, z=z, desc=desc);
	z.stem <- obj$z;
	aff <- obj$affected;

## Ensure that neither partition is empty; can occur with "node" or "both" cutting. If so, kill it.
	if (sum(z.stem[,"partition"] == aff[1]) == 0 || sum(z.stem[,"partition"] == aff[2]) == 0) {
		fit$lnLik <- -Inf;
		fit$aic <- Inf;
		fit$aicc <- Inf;
		return(fit);
	}
	
## Everything is cool; proceed.
	sp <- op[aff[1],]; # Use previously fit parameter values from clade that is currently being split
	
## first, consider diminshed clade. may result in a clade that has been cached previously
	dimClade <- z.stem[z.stem[,"partition"] == aff[1],,drop=FALSE];
	
	fit1 <- getOptimalModelFlavour(z=dimClade, sp=sp, model=model, fixPar=fixPar, criterion=criterion);
	
## Sister clade
	newClade <- z.stem[z.stem[,"partition"] == aff[2],,drop=FALSE];
	fit2 <- getOptimalModelFlavour(z=newClade, sp=sp, model=model, fixPar=fixPar, criterion=criterion);

	cut.at <- "stem";
	
	op[aff[1],] <- fit1$par; # Replace parameters with new values for diminished clade
	
	fit$model[aff[1]] <- fit1$model;
	fit$par <- rbind(op, fit2$par);
	fit$lnLik.part[aff] <- c(fit1$lnLik, fit2$lnLik); # Replace parameters with new values for diminished clade
	fit$split.at <- c(fit$split.at, node);
	fit$lnLik <- sum(fit$lnLik.part);
	
	model.fit <- calculateModelFit(fit=fit, z=z);
	
	fit$aic <- model.fit[1];
	fit$aicc <- model.fit[2];
	fit$num.par <- model.fit[3];
	fit$cut.at <- c(fit$cut.at, cut.at);
	fit$model <- c(fit$model, fit2$model);
	fit$z <- z.stem;
	fit$step <- rbind(fit$step, c("add", node));
	
	return(list(fit));
}

findMrca <- function(phy, tips) {
	tt <- match(tips, phy$tip.label);
	if (sum(!is.na(tt)) != length(tips)) {
		for (i in 1:length(tt)) {
			if (is.na(tt[i])) {
				cat("\n\nWarning: taxon '", tips[i], "' is not found in the tree. Stopping.\n", sep="");
				stop;
			}
		}
	}
	return(getMRCA(phy, tt)); # ape function
}

# 'tips' contains 2 (or more) tip labels, used for defining a MRCA node
extractSubTree <- function (phy, tips=NULL, node=NULL, richness=NULL, verbose=T) {
	if (is.null(tips) && is.null(node) && node != "root") stop("\n\nWarning: need to provide either a node number or set of 2 taxa whose MRCA is the desired node. Stopping.\n");
	
	if (!is.null(tips)) {
		node <- findMrca(phy, tips);
	} else if (node == "root") {
		node <- length(phy$tip.label) + 1;
	}
#	cat(node, "\n");
	phy <- extract.clade(phy, node); # ape function
	
## Prune tree with richness data as necessary.
	phyData <- prepareData(phy=phy, richness=richness, verbose=verbose);
	
	return(phyData);
}


sisterFitSummary <- function (results, criterion, fixPar, plotSurface=FALSE) {
	summarizeSingleTree <- function ()
	{
	## Get profile likelihoods for base and sister model separately
		models <- results$models;
		z <- results$z;
		
		base <- models[[1]];
		sis <- models[[2]];
		
		baseProfLikes <- getProfileLikelihoods(z=base$z, parm=base$par, models=base$model, fixPar=fixPar);
		sisterProfLikes <- getProfileLikelihoods(z=sis$z, parm=sis$par, models=sis$model, fixPar=fixPar);
		profLikes <- rbind(baseProfLikes, sisterProfLikes);
		par <- rbind(base$par, sis$par);
		modelFlavour <- c(base$model, sis$model[1], sis$model[2]);
		richness <- c(sum(z[,"n.t"], na.rm=T), sum(sis$z[sis$z[,"partition"] == 1,"n.t"], na.rm=T), sum(sis$z[sis$z[,"partition"] == 2,"n.t"], na.rm=T));
		
		clades <- c("Base", "Clade1", "Clade2");
		indLikes <- c(base$lnLik, sis$lnLik.part[1], sis$lnLik.part[2]);
		modelLikes <- c(base$lnLik, sis$lnLik);
		numPar <- c(base$num.par, sis$num.par);
		
	## Calculate model weights; 'fit' contains just AIC scores. Returns: data.frame(fit=fit, delta=delta, w=w);
		scores <- c(as.numeric(models[[1]][criterion]), as.numeric(models[[2]][criterion]));
		weights <- calculateModelWeights(fit=scores);
		
		if (all(modelFlavour == "yule")) {
			parameterSummary <- cbind(Clade=clades, Richness=richness, lnLik=indLikes, model=modelFlavour, r=par[,1], profLikes[,1:2]);
		} else {
			parameterSummary <- cbind(Clade=clades, Richness=richness,lnLik=indLikes, model=modelFlavour, par, profLikes);
		}
		
		modelSummary <- cbind(Model=c("Base", "Sisters"), lnLik=modelLikes, numPar=numPar, weights);
		
		cat("\nParameter Summary:\n\n");
		print(parameterSummary);
		
		cat("\nModel Fit Summary:\n\n");
		print(modelSummary);
		
		summary <- list(parSumm=parameterSummary, modSumm=modelSummary);
		return(summary);
	}
	
	summarizeMultipleTrees <- function () {
		cladeNames <-c("Base", "Clade1", "Clade2");
		basePar <- NULL;
		clade1Par <- NULL;
		clade2Par <- NULL;
		AICscores <- list();
		
		for (i in 1:length(results)) {
			models <- results[[i]]$models;
			base <- models[[1]];
			sis <- models[[2]];
			basePar <- rbind(basePar, base$par);
			clade1Par <- rbind(clade1Par, sis$par[1,]);
			clade2Par <- rbind(clade2Par, sis$par[2,]);
			AICscores <- c(AICscores, list(c(as.numeric(base[criterion]), as.numeric(sis[criterion]))));
		}
		summBasePar <- summarizeParameters(basePar);
		summClade1Par <- summarizeParameters(clade1Par);
		summClade2Par <- summarizeParameters(clade2Par);
		summPar <- rbind(summBasePar, summClade1Par, summClade2Par);
		summPar <- as.data.frame(summPar); rownames(summPar) <- NULL;
		summPar <- cbind(Clade=cladeNames, summPar);
		
		if (all(is.na(summPar[,"Mean.eps"]))) { # yule; forget epsilon
			summPar <- summPar[,c(1:5)];
			parameterSummary <- summPar[,c(1:5)];
		} else {
			parameterSummary <- summPar;
		}
		
		modelSummary <- summarizeMultipleSisterModelFit(AICscores);
		
		cat("\nResults from the analysis of ", length(results), " phylogenetic trees.\n", sep="");
		
		cat("\nParameter Summary:\n\n");
		print(parameterSummary);
		
		cat("\nModel Fit Summary:\n\n");
		print(modelSummary);
		
		summary <- list(parSumm=parameterSummary, modSumm=modelSummary);
		return(summary);
	}
	
	if (length(results) == 2) { # single tree
		summary <- summarizeSingleTree();
		summary$nTrees <- 1;
	} else { # multiple trees
		summary <- summarizeMultipleTrees();
		summary$nTrees <- length(results);
	}
	
	return(summary);
}


## Passed-in 'p' will be a two-column matrix
summarizeParameters <- function(p) {
## Want: mean, min, max, st.dev
	p1 <- c(Mean.r=mean(p[,1]), Min=min(p[,1]), Max=max(p[,1]), StDev=sd(p[,1]));
	p2 <- c(Mean.eps=mean(p[,2]), Min=min(p[,2]), Max=max(p[,2]), StDev=sd(p[,2]));
	
	parSummary <- c(p1, p2);
	
	return(parSummary);
}


## 'AICscores' is a list of AIC scores for base and sister models
summarizeMultipleSisterModelFit <- function(AICscores) {
	delta <- NULL;
	w <- NULL;
	
	for (i in 1:length(AICscores)) {
		x <- calculateModelWeights(AICscores[[i]]);
		delta <- rbind(delta, x$delta);
		w <- rbind(w, x$w);
	}
	
	res <- matrix(nrow=2, ncol=6);
	colnames(res) <- c("Mean.delta", "SD.delta", "Mean.w", "SD.w", "Min.w", "Max.w");
	
	res[1,] <- c(mean(delta[,1]), sd(delta[,1]), mean(w[,1]), sd(w[,1]), min(w[,1]), max(w[,1]));
	res[2,] <- c(mean(delta[,2]), sd(delta[,2]), mean(w[,2]), sd(w[,2]), min(w[,2]), max(w[,2]));
	res <- as.data.frame(res);
	
	best <- sum(delta[,1] == 0) / length(delta[,1]);
	best <- c(best, 1-best);
	
	modelSummary <- cbind(Model=c("Base", "Sisters"), Prop.best=best, res);
	
	return(modelSummary);
}

print.sisterFit <- function(x, ...) {
	if (x$nTrees > 1) {cat("\nResults from the analysis of ", x$nTrees, " phylogenetic trees.\n", sep="");}
	
	cat("\nParameter Summary:\n\n");
	print(x$parameterSummary);
	if (x$nTrees == 1) {cat("\n95% confidence intervals calculated from profile likelihoods.\n");}
	
	cat("\nModel Fit Summary:\n\n");
	print(x$modelSummary);
	cat("\n");
}
