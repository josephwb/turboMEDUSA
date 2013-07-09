## Only used for base model
medusaMLFitBase <- function (z, sp, model, fixPar, criterion) {
	fit <- getOptimalModelFlavour(z=z, sp=sp, model=model, fixPar=fixPar, criterion=criterion);
	model.fit <- calculateModelFit(fit=fit, z=z);
	
	steps <- c("add", as.numeric(min(z[,"anc"])));
	
	return(list(par=matrix(fit$par, nrow=1, dimnames=list(NULL,c("r", "epsilon"))), lnLik.part=fit$lnLik, lnLik=fit$lnLik,
		split.at=min(z[,"anc"]), aic=round(model.fit[1], digits=7), aicc=round(model.fit[2], digits=7), num.par=model.fit[3],
		cut.at="node", model=fit$model, z=z, step=matrix(steps, nrow=1, dimnames=list(NULL,c("step", "node")))));
}


## Split the edge matrix 'z' by adding a partition at node 'node'.
## The list 'desc' is a list of descendants (see makeCacheMedusa).
## Returns a list with elements:
##   z: new medusa matrix, with the new partition added
##   affected: indices of the partitions affected by the split (n == 2).
medusaSplitStem <- function (node, z, desc, extract=FALSE) {
	part <- z[,"partition"];
	base <- min(part[z[,"anc"] == node | z[,"dec"] == node]);
	tag <- max(part) + 1;
	i <- desc[[node]];
	idx <- i[part[i] == base];
	z[idx,"partition"] <- tag;
	if (extract) {z <- z[idx,,drop=FALSE];}
	return(list(z=z, affected=c(unique(part[idx]), tag)));
}

medusaSplitNode <- function (node, z, desc, extract=FALSE) {
	part <- z[,"partition"];
	base <- min(part[z[,"anc"] == node]);
	tag <- max(part) + 1;
	i <- desc[[node]];
	idx <- i[part[i] == base];
	z[idx,"partition"] <- tag;
	if (extract) {z <- z[idx,,drop=FALSE];}
	return(list(z=z, affected=c(unique(part[idx]), tag)));
}

## The general function for when the flavour of shiftCut is unknown
medusaSplit <- function (node, z, desc, shiftCut) {
	descendants <- desc$stem; # will be correct half the time
	if (shiftCut == "node") {descendants <- desc$node;}
	
	part <- z[,"partition"];
	base <- min(part[z[,1] == node | z[,2] == node]);
	tag <- max(part) + 1;
	i <- descendants[[node]];
	idx <- i[part[i] == base];
	z[idx,"partition"] <- tag;
	return(list(z=z, affected=c(unique(part[idx]), tag)));
}


## Pre-fit values for pendant edges; DON'T recalculate later; should account for ~25% of all calculations
## Also cache values for virgin nodes; useful until subsetted, especially for large trees.
## Tips with richness of 1 have not undergone anything along edge. r = 0, lik = 1.
prefitTips <- function (pend.nodes, z, sp, model, fixPar, criterion, mc, numCores) {
	# for case where all tips have richness of 1 i.e. a fully sampled tree.
	if (all(z[z[,"dec"] %in% pend.nodes,"n.t"] == 1) && (model == "yule" || model == "mixed")) {
		fit <- list(list(par=c(0, NA), lnLik=0, model="yule"));
		fit <- rep(fit, length(pend.nodes));
		return(fit);
	}
	
	# yule will always have a better AIC for tip, regardless of richness
	if (model == "mixed" || model == "yule") {
		if (mc) {
			tips <- mclapply(pend.nodes, medusaMLPrefitTip, z=z, sp=sp,
				model="yule", fixPar=fixPar, criterion=criterion, mc.cores=numCores);
		} else {
			tips <- lapply(pend.nodes, medusaMLPrefitTip, z=z, sp=sp,
				model="yule", fixPar=fixPar, criterion=criterion);
		}
	} else {
		if (mc) {
			tips <- mclapply(pend.nodes, medusaMLPrefitTip, z=z, sp=sp,
				model=model, fixPar=fixPar, criterion=criterion, mc.cores=numCores);
		} else {
			tips <- lapply(pend.nodes, medusaMLPrefitTip, z=z, sp=sp,
				model=model, fixPar=fixPar, criterion=criterion);
		}
	}
	return(tips);
}

medusaMLPrefitTip <- function (node, z, sp, model, fixPar, criterion) {
	z.tip <- z[z[,"dec"] == node,,drop=FALSE];
	
	if ((z.tip[,"n.t"] == 1) && (model == "yule" || model == "mixed")) {
		return(list(par=c(0, NA), lnLik=0, model="yule"));
	}
	
	# tips are always better fit by yule. don't bother with BD.
	if (model == "yule" || model == "mixed") {
		if (z.tip[,"n.t"] == 1) { # single tip. nothing has happened along edge, so r = 0
			return(list(par=c(0, NA), lnLik=0, model="yule"));
		} else { # unresolved clade
			# no t.len information available, only depth
			# MLE birth rate is: log(n.t) / depth
			r <- as.numeric(log(z.tip[,"n.t"]) / z.tip[,"t.len"]);
			lik <- as.numeric(-z.tip[,"t.len"] * r + (z.tip[,"n.t"] - 1) * log(1 - exp(-z.tip[,"t.len"] * r)));
			return(list(par=c(r, NA), lnLik=lik, model="yule"));
		}
	} else {
		fit <- getOptimalModelFlavour(z=z.tip, sp=sp, model=model, fixPar=fixPar, criterion=criterion);
		return(fit);
	}
}

medusaMLPrefitNode <- function (node, z, desc, sp, model, fixPar, criterion) {
	z.node <- medusaSplitNode(node=node, z=z, desc=desc, extract=TRUE)$z;
	if (all(z.node[,"n.t"] == 1, na.rm=T)) {
		fit <- getOptimalModelFlavour(z=z.node, sp=sp, model=model, fixPar=fixPar, criterion=criterion, simple=TRUE);
	} else {
		fit <- getOptimalModelFlavour(z=z.node, sp=sp, model=model, fixPar=fixPar, criterion=criterion);
	}
	return(fit);
}

medusaMLPrefitStem <- function (node, z, desc, sp, model, fixPar, criterion) {
	z.stem <- medusaSplitStem(node=node, z=z, desc=desc, extract=TRUE)$z;
	if (all(z.stem[,"n.t"] == 1, na.rm=T)) {
		fit <- getOptimalModelFlavour(z=z.stem, sp=sp, model=model, fixPar=fixPar, criterion=criterion, simple=TRUE);
	} else {
		fit <- getOptimalModelFlavour(z=z.stem, sp=sp, model=model, fixPar=fixPar, criterion=criterion);
	}
	return(fit);
}

# solution known when subtree is completely sampled
simpleYule <- function (z) {
	n.int <- sum(is.na(z[,"n.t"])); # the number of speciation events
	if (n.int == 0) {
		return(list(par=c(0, NA), lnLik=0));
	}
	sum.t <- sum(z[,"t.len"]);
	r <- (n.int) / sum.t;
	lik <- n.int * log(r) - r * sum.t;
	par <- c(r, NA);
	return(list(par=c(r, NA), lnLik=lik));
}

## When model == mixed, fit both and find optimal flavour
getOptimalModelFlavour <- function (z, sp, model, fixPar, criterion, simple=FALSE) {
	fit.bd <- NULL;
	fit.yule <- NULL;
	fit <- NULL;
	
	if (model == "yule" | model == "mixed") {
		if (simple) {
			fit.yule <- simpleYule(z = z);
		} else {
			fit.yule <- medusaMLFitPartition(z=z, sp=sp, model="yule");
		}
		fit.yule$model <- "yule";
	}
	if (model == "bd" | model == "mixed") {
		if (is.na(sp[2])) {sp[2] <- 0.5;}
		fit.bd <- medusaMLFitPartition(z=z, sp=sp, model="bd");
		fit.bd$model <- "bd";
	}
	if (model != "mixed" && model != "bd" && model != "yule") { # i.e. the constrained models
		fit <- medusaMLFitPartition(z=z, sp=sp, model=model, fixPar=fixPar);
		fit$model <- model;
		return(fit);
	}
	
## Figure out which model fits best
	if (is.null(fit.bd)) {
		fit <- fit.yule;
	} else if (is.null(fit.yule)) {
		fit <- fit.bd;
	} else {
## Considering both models
		fit <- getBestPartialModel(fit1=fit.yule, fit2=fit.bd, z=z, criterion=criterion);
	}
	return(fit);
}


## Memory issues required a restructuring. medusaMLUpdate finds intermediate results *within* a particular model size.
##  - only pertinent values are recorded: split.at, aic, aicc, cut.at.
##  - the optimal model is later reoptimized to get all required values: above + likelihood, z, etc.
## 'fit' contains parameter values from previous model, used to initialize subsequent model.
## Pass in pre-fitted values for pendant edges and virgin nodes (in 'prefit'); DON'T recalculate.
## Need to consider the possibility of birth-death, yule, mixed, or constrained models.
## Need to consider where shft is placed (shiftCut). Placement affects not only new clade, but
## also the size of the split clade. Only relevant if shiftCut = "both".
#medusaMLUpdate <- function (node, z, desc, fit, prefit, num.tips, root.node, model, fixPar, criterion, shiftCut, preserveModelFlavour)
medusaMLUpdate <- function (node, z, desc, fit, prefit, root.node, model, fixPar, criterion, shiftCut, preserveModelFlavour) {
## various combinations possible
	fit1.stem <- NULL;
	fit1.node <- NULL;
	fit2.stem <- NULL;
	fit2.node <- NULL;
	cut.at <- NULL;
	
	sp <- NULL;
	aff <- NULL;
	op <- fit$par; # store previously fit parameter values
	cool <- TRUE;
	
	fit1 <- NULL;
	fit2 <- NULL;
	
	if (shiftCut == "stem" | shiftCut == "both" | node < root.node) { # can enter on "node" if a tip
## First, diminshed clade
		obj.stem <- medusaSplitStem(node=node, z=z, desc=desc$stem);
		z.stem <- obj.stem$z;
		aff <- obj.stem$affected;

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
		
## check if diminished clade conforms to a cached clade
		x <- which(dimClade[,"anc"] == min(dimClade[,"anc"]));
		if (length(x) == 1) { # a cut-at-stem scenario
			y <- as.numeric(dimClade[x, "dec"]);
			if (length(unique(dimClade[(dimClade[,"dec"] < root.node),"dec"])) == prefit$num.tips[[y]]) {
				if (y < root.node)
				{
					fit1.stem <- prefit$tips[[y]];
				} else {
					fit1.stem <- prefit$virgin.nodes$stem[[y - root.node]];
				}
			}
		}
				
		if (is.null(fit1.stem)) {
			if (model == "mixed") {
				if (preserveModelFlavour) { ## In mixed models, may want to conserve flavour of previously fit model
					if (sum(!is.na(sp)) < 2) # yule
					{
						fit1.stem <- medusaMLFitPartition(z=dimClade, sp=sp, model="yule");
					} else {
						fit1.stem <- medusaMLFitPartition(z=dimClade, sp=sp, model="bd");
					}
				} else {
## consider both model flavours
					fit1.stem <- getOptimalModelFlavour(z=dimClade, sp=sp, model=model, fixPar=fixPar, criterion=criterion);
				}
			} else {
				fit1.stem <- medusaMLFitPartition(z=dimClade, sp=sp, model=model, fixPar=fixPar);
				fit1.stem$model <- model;
			}
		}
		
## Second, new clade
		if (node < root.node) { # tip, already calculated
			fit2.stem <- prefit$tips[[node]];
		} else if (length(unique(z.stem[(z.stem[,"partition"] == aff[2] & z.stem[,"dec"] < root.node),"dec"])) == prefit$num.tips[[node]]) {
## vigin node, already calculated
			fit2.stem <- prefit$virgin.nodes$stem[[node - root.node]];
		} else {
## novel shift
			newClade <- z.stem[z.stem[,"partition"] == aff[2],,drop=FALSE];
			fit2.stem <- getOptimalModelFlavour(z=newClade, sp=sp, model=model, fixPar=fixPar, criterion=criterion);
		}
	}
	
	if ((shiftCut == "node" || shiftCut == "both") && (node > root.node)) { # never enter if tip
## First, diminshed clade
		obj.node <- medusaSplitNode(node=node, z=z, desc=desc$node);
		z.node <- obj.node$z;
		aff <- obj.node$affected;
		
## Need to check if cut is valid. May be inadmissable because of pattern of previous breaks (especially with shiftCut=both)
		if (is.na(aff[1]) || is.na(aff[2])) {
			cool <- FALSE;
		}		
## Check that partition is not empty; can occur with "node" or "both" cutting.
		if (sum(z.node[,"partition"] == aff[1]) == 0 || sum(z.node[,"partition"] == aff[2]) == 0 || !cool) {
			fit$lnLik <- -Inf;
			fit$aic <- Inf;
			fit$aicc <- Inf;
			return(fit);
		}
		
## Everything is cool; proceed.
		sp <- op[aff[1],]; # Use previously fit parameter values from clade that is currently being split
		
## first, consider diminshed clade. may result in a clade that has been cached
		dimClade <- z.node[z.node[,"partition"] == aff[1],,drop=FALSE];
		
		if (model == "mixed") {
			if (preserveModelFlavour) { ## In mixed models, may want to conserve flavour of previously fit model
				if (sum(!is.na(sp)) < 2) # yule
				{
					fit1.node <- medusaMLFitPartition(z=dimClade, sp=sp, model="yule"); # should this change? probably.
				} else {
					fit1.node <- medusaMLFitPartition(z=dimClade, sp=sp, model="bd");
				}
			} else {
## consider both model flavours
				fit1.node <- getOptimalModelFlavour(z=dimClade, sp=sp, model=model, fixPar=fixPar, criterion=criterion);
			}
		} else {
			fit1.node <- medusaMLFitPartition(z=dimClade, sp=sp, model=model, fixPar=fixPar);
			fit1.node$model <- model;
		}
## Second, new clade
		if (length(unique(z.node[(z.node[,"partition"] == aff[2] & z.node[,"dec"] < root.node),"dec"])) == prefit$num.tips[[node]]) {
## vigin node, already calculated
			fit2.node <- prefit$virgin.nodes$node[[node - root.node]];
		} else {
## novel shift
			newClade <- z.node[z.node[,"partition"] == aff[2],,drop=FALSE];
			fit2.node <- getOptimalModelFlavour (z=newClade, sp=sp, model=model, fixPar=fixPar, criterion=criterion);
		}
	}
	
## Now, figure out which shift position is optimal	
	if (is.null(fit2.node)) {
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
		stem.val <- list(lnLik=stem.lik, par=stem.par, model=model);
		stem.fit <- calculateModelFit(fit=stem.val, z=z);
		
		node.lik <- (fit1.node$lnLik + fit2.node$lnLik);
		node.par <- rbind(fit1.node$par, fit2.node$par)
		node.val <- list(lnLik=node.lik, par=node.par, model=model);
		node.fit <- calculateModelFit(fit=node.val, z=z);
		
		if (criterion == "aic") {element <- 1;} else {element <- 2;}
		
		if (stem.fit[[element]] < node.fit[[element]]) {
			fit1 <- fit1.stem;
			fit2 <- fit2.stem;
			cut.at <- "stem";
		} else {
			fit1 <- fit1.node;
			fit2 <- fit2.node;
			cut.at <- "node";
		}
	}
	op[aff[1],] <- fit1$par; # Replace parameters with new values for diminished clade
	
	if (!preserveModelFlavour) {fit$model[aff[1]] <- fit1$model;} # update altered model
	
	fit$par <- rbind(op, fit2$par);
	fit$lnLik.part[aff] <- c(fit1$lnLik, fit2$lnLik); # Replace parameters with new values for diminished clade
	fit$lnLik <- sum(fit$lnLik.part);
	
	model.fit <- calculateModelFit(fit=fit, z=z);
	
	intFit <- list(aic=model.fit[1], aicc=model.fit[2], split.at=node, cut.at=cut.at);
	
	return(intFit);
}


## this was the original general function, but presented a significant memory issue for large trees.
## re-fit optimal node shift to get all required values; significant memory issue for large trees if done otherwise.
## shiftCut and node are fixed. model is as well, but it seems tricky...
medusaFitOptimal <- function (node, z, desc, fit, prefit, root.node, model, fixPar, criterion, shiftCut, preserveModelFlavour) {
## various combinations possible
	fit1 <- NULL;
	fit1 <- NULL;
	fit2 <- NULL;
	fit2 <- NULL;
	cut.at <- shiftCut;
	
	sp <- NULL;
	aff <- NULL;
	op <- fit$par; # store previously fit parameter values
	cool <- TRUE;
	
	fit1 <- NULL;
	fit2 <- NULL;
	
	if (shiftCut == "stem") {
## First, diminshed clade
		obj <- medusaSplitStem(node=node, z=z, desc=desc$stem);
		z <- obj$z;
		aff <- obj$affected;

## Ensure that neither partition is empty; can occur with "node" or "both" cutting. If so, kill it.
		if (sum(z[,"partition"] == aff[1]) == 0 || sum(z[,"partition"] == aff[2]) == 0) {
			fit$lnLik <- -Inf;
			fit$aic <- Inf;
			fit$aicc <- Inf;
			return(fit);
		}
		
## Everything is cool; proceed.
		sp <- op[aff[1],]; # Use previously fit parameter values from clade that is currently being split
		
## first, consider diminshed clade. may result in a clade that has been cached previously
		dimClade <- z[z[,"partition"] == aff[1],,drop=FALSE];
		
## check if diminished clade conforms to a cached clade
		x <- which(dimClade[,"anc"] == min(dimClade[,"anc"]));
		if (length(x) == 1) { # a cut-at-stem scenario
			y <- as.numeric(dimClade[x, "dec"]);
			if (length(unique(dimClade[(dimClade[,"dec"] < root.node),"dec"])) == prefit$num.tips[[y]]) {
				if (y < root.node) {
					fit1 <- prefit$tips[[y]];
				} else {
					fit1 <- prefit$virgin.nodes$stem[[y - root.node]];
				}
			}
		}
				
		if (is.null(fit1)) {
			if (model == "mixed") {
				if (preserveModelFlavour) ## In mixed models, may want to conserve flavour of previously fit model
				{
					if (sum(!is.na(sp)) < 2) { # yule
						fit1 <- medusaMLFitPartition(z=dimClade, sp=sp, model="yule");
					} else {
						fit1 <- medusaMLFitPartition(z=dimClade, sp=sp, model="bd");
					}
				} else {
## consider both model flavours
					fit1 <- getOptimalModelFlavour(z=dimClade, sp=sp, model=model, fixPar=fixPar, criterion=criterion);
				}
			} else {
				fit1 <- medusaMLFitPartition(z=dimClade, sp=sp, model=model, fixPar=fixPar);
				fit1$model <- model;
			}
		}
		
## Second, new clade
		if (node < root.node) { # tip, already calculated
			fit2 <- prefit$tips[[node]];
		} else if (length(unique(z[(z[,"partition"] == aff[2] & z[,"dec"] < root.node),"dec"])) == prefit$num.tips[[node]]) {
## vigin node, already calculated
			fit2 <- prefit$virgin.nodes$stem[[node - root.node]];
		} else {
## novel shift
			newClade <- z[z[,"partition"] == aff[2],,drop=FALSE];
			fit2 <- getOptimalModelFlavour(z=newClade, sp=sp, model=model, fixPar=fixPar, criterion=criterion);
		}
	} else if (shiftCut == "node") {
## First, diminshed clade
		obj <- medusaSplitNode(node=node, z=z, desc=desc$node);
		z <- obj$z;
		aff <- obj$affected;
		
## Need to check if cut is valid. May be inadmissable because of pattern of previous breaks
		if (is.na(aff[1]) || is.na(aff[2])) {
			cool <- FALSE;
		}		
## Check that partition is not empty; can occur with "node" cutting.
		if (sum(z[,"partition"] == aff[1]) == 0 || sum(z[,"partition"] == aff[2]) == 0 || !cool) {
			fit$lnLik <- -Inf;
			fit$aic <- Inf;
			fit$aicc <- Inf;
			return(fit);
		}
		
## Everything is cool; proceed.
		sp <- op[aff[1],]; # Use previously fit parameter values from clade that is currently being split
		
## first, consider diminshed clade. may result in a clade that has been cached
		dimClade <- z[z[,"partition"] == aff[1],,drop=FALSE];
		
		if (model == "mixed") {
			if (preserveModelFlavour) { ## In mixed models, may want to conserve flavour of previously fit model
				if (sum(!is.na(sp)) < 2) { # yule
					fit1 <- medusaMLFitPartition(z=dimClade, sp=sp, model="yule"); # should this change? probably.
				} else {
					fit1 <- medusaMLFitPartition(z=dimClade, sp=sp, model="bd");
				}
			} else {
## consider both model flavours
				fit1 <- getOptimalModelFlavour(z=dimClade, sp=sp, model=model, fixPar=fixPar, criterion=criterion);
			}
		} else {
			fit1 <- medusaMLFitPartition(z=dimClade, sp=sp, model=model, fixPar=fixPar);
			fit1$model <- model;
		}
		
## Second, new clade
		if (length(unique(z[(z[,"partition"] == aff[2] & z[,"dec"] < root.node),"dec"])) == prefit$num.tips[[node]]) {
## vigin node, already calculated
			fit2 <- prefit$virgin.nodes$node[[node - root.node]];
		} else {
## novel shift
			newClade <- z[z[,"partition"] == aff[2],,drop=FALSE];
			fit2 <- getOptimalModelFlavour (z=newClade, sp=sp, model=model, fixPar=fixPar, criterion=criterion);
		}
	}
	
	op[aff[1],] <- fit1$par; # Replace parameters with new values for diminished clade
	
	if (!preserveModelFlavour) {fit$model[aff[1]] <- fit1$model;} # update altered model
	
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
	
	return(fit);
}


## sp = initializing values for r & epsilon. fixPar = fixed value for constrained model.
## Default values should never be used (except for first model), as the values from the previous model are passed in.
## The 'suppressWarnings' catches harmless warnings for when extreme parameter values are tried, 
## generating -Inf log-likelihoods. Depends on sclae of the tree.
medusaMLFitPartition <- function (z, sp=c(0.05, 0.5), model, fixPar=NULL) {
#	new.part <- z[z[,"partition"] == partition.id,,drop=FALSE]; # now subsetted earlier
	new.part <- z;
	
# Construct likelihood function:
	lik <- makePartitionLikelihood(partition=new.part, model=model, fixPar=fixPar);
	foo <- function (x) {-lik(pars=exp(x));} # work with parameters in log-space to preserve precision
	
	if (model == "bd") {
		fit <- optim(fn=foo, par=log(sp), method="N", control=list(maxit=5000));
		
		if (fit$convergence != 0) {cat("\nDidn't converge. Shit. Convergence = ", fit$convergence, "\n");}
		
		return(list(par=exp(fit$par), lnLik=-fit$value));
	}
		
## grab values to get an intelligent upper bound on b or r
		node.richness  <- sum(new.part[,"n.t"], na.rm=TRUE);
		depth <- max(new.part[,"t.0"]);
		
# use different intervals based on model flavour
	if (model == "yule") {
		
		n.int <- sum(is.na(new.part[,"n.t"]));
		
		if ((n.int == 0) && all(new.part[,"n.t"] == 1, na.rm = TRUE)) {
			return(list(par=c(0, NA), lnLik=-Inf));
		}
		
		maxVal <- (log(node.richness) / depth) * 5;
		if (node.richness <= 1) {maxVal <- 1e-5;}
		
		suppressWarnings(fit <- optimize(f=foo, interval=c(-25, log(maxVal))));
		par <- c(exp(fit$minimum), NA);
		
		while (par[1]/maxVal > 0.95) { # crash against boundary; doesn't seem to get used...
			maxVal <- par[1] * 3;
#			cat("Hit boundary\n")
			suppressWarnings(fit <- optimize(f=foo, interval=c(log(par[1]/2), log(maxVal))));
			par <- c(exp(fit$minimum), NA);
		}
		return(list(par=par, lnLik=-fit$objective));

		
	} else if (model == "fixedD") {
## this is the messiest model flavour; reasonable range depends on magnitude of fixD
		
		x <- (log(node.richness) / depth) * 5;
		if (node.richness <= 1) {x <- 0.01;} # MLE for single tip will be r ~ 0
		maxVal <- fixPar + x;
		
		suppressWarnings(fit <- optimize(f=foo, interval=c(log(fixPar), log(maxVal))));
		b <- exp(fit$minimum);
		par <- c((b - fixPar), (fixPar / b));
		
		while (fit$objective == Inf) {
			cat("maxVal =", maxVal, "; b =", b, "; r =", b - fixPar, "; epsilon =", fixPar/b, "\n")
			
			x <- x * 0.75;
			maxVal <- fixPar + x;
			suppressWarnings(fit <- optimize(f=foo, interval=c(log(fixPar), log(maxVal))));
			b <- exp(fit$minimum);
			par <- c((b - fixPar), (fixPar / b));
		}
		
		return(list(par=par, lnLik=-fit$objective));
		
	} else if (model == "fixedEpsilon") {
		
		maxVal <- (log(node.richness) / depth) * 5;
		if (node.richness <= 1) {maxVal <- 1;}
		
		suppressWarnings(fit <- optimize(f=foo, interval=c(-25, log(maxVal))));
		par <- c(exp(fit$minimum), sp[2]);
		
		while (par[1]/maxVal > 0.9) { # crash against boundary
			maxVal <- par[1] * 3;
			suppressWarnings(fit <- optimize(f=foo, interval=c(log(par[1]/2), log(maxVal))));
			par <- c(exp(fit$minimum), sp[2]);
		}
		
		return(list(par=par, lnLik=-fit$objective));
		
	} else if (model == "fixedB") {
		
		suppressWarnings(fit <- optimize(f=foo, interval=c(-50, log(fixPar))));
		d <- exp(fit$minimum);
		par <- c((fixPar - d), (d / fixPar));
		return(list(par=par, lnLik=-fit$objective));
		
	} else if (model == "fixedR") {
		
		suppressWarnings(fit <- optimize(f=foo, interval=c(-50, 0)));
		par <- c(sp[1], exp(fit$minimum));
		return(list(par=par, lnLik=-fit$objective));
	}
}


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

## makePartitionLikelihood: generate a likelihood function for a single partition.
makePartitionLikelihood <- function (partition, model, fixPar) {
	i.int <- is.na(partition[,"n.t"]);
	i.pend <- !(i.int);
	n.int <- sum(i.int);
	
	simple <- FALSE;
	if (all(partition[,"n.t"] == 1, na.rm = TRUE)) {simple <- TRUE;}
	
	if (simple) {
		if (model == "yule") {
			sum.t <- sum(partition[,"t.len"]);
			f <- function(pars) {
				if (pars <= 0) return(-Inf);
				r <- pars[1];
				return(n.int * log(r) - r * sum.t);
			}
			return(f)
		} 
		
		sum.int.t.len <- sum(partition[i.int,"t.len"]);
		int.t.0 <- partition[i.int,"t.0"];
		pend.t.len <- partition[i.pend,"t.len"];
		
		if (model == "bd") {
			f <- function(pars) {
				r <- pars[1];
				epsilon <- pars[2];
				if (r <= 0 | epsilon < 0 | epsilon >= 1) {return(-Inf);}
				l.int <- n.int * log(r) - r * sum.int.t.len - sum(log(1 - (epsilon * exp(-r * int.t.0))));
				l.pend <- sum(log(1 - epsilon) - log(exp(r * pend.t.len) - epsilon));
				return(l.int + l.pend);
			}
			return(f)
		} else if (model == "fixedD") {
			fixD <- fixPar;
			f <- function(pars) {
				b <- pars[1];
				r <- b - fixD;
				epsilon <- fixD/b;
				if (r <= 0 | epsilon < 0 | epsilon >= 1) {return(-Inf);}
				l.int <- n.int * log(r) - r * sum.int.t.len - sum(log(1 - (epsilon * exp(-r * int.t.0))));
				l.pend <- sum(log(1 - epsilon) - log(exp(r * pend.t.len) - epsilon));
				return(l.int + l.pend);
			}
			return(f)
		} else if (model == "fixedEpsilon") {
			epsilon <- fixPar;
			f <- function(pars) {
				r <- pars[1];
				if (r <= 0 | epsilon < 0 | epsilon >= 1) {return(-Inf);}
				l.int <- n.int * log(r) - r * sum.int.t.len - sum(log(1 - (epsilon * exp(-r * int.t.0))));
				l.pend <- sum(log(1 - epsilon) - log(exp(r * pend.t.len) - epsilon));
				return(l.int + l.pend);
			}
			return(f)
		} else if (model == "fixedR") {
			r <- fixPar;
			f <- function(pars) {
				epsilon <- pars[1];
				if (r <= 0 | epsilon < 0 | epsilon >= 1) {return(-Inf);}
				l.int <- n.int * log(r) - r * sum.int.t.len - sum(log(1 - (epsilon * exp(-r * int.t.0))));
				l.pend <- sum(log(1 - epsilon) - log(exp(r * pend.t.len) - epsilon));
				return(l.int + l.pend);
			}
			return(f)
		} else if (model == "fixedB") {
			fixB <- fixPar;
			f <- function(pars) {
				d <- pars[1];
				r <- fixB - d;
				epsilon <- d/fixB;
				if (r <= 0 | epsilon < 0 | epsilon >= 1) {return(-Inf);}
				l.int <- n.int * log(r) - r * sum.int.t.len - sum(log(1 - (epsilon * exp(-r * int.t.0))));
				l.pend <- sum(log(1 - epsilon) - log(exp(r * pend.t.len) - epsilon));
				return(l.int + l.pend);
			}
			return(f)
		}
	} else { # not simple: one or more tips represent > 1 extant species
		sum.int.t.len <- sum(partition[i.int,"t.len"]);
		pend.t.len <- partition[i.pend,"t.len"];
		pend.n.t.minus.1 <- partition[i.pend,"n.t"] - 1;
		
		if (model == "yule") {
			sum.pend.t.len <- sum(pend.t.len);
			f <- function(pars) {
				r <- pars[1];
				if (r <= 0) return(-Inf);
				return(n.int * log(r) - r * sum.int.t.len + sum(-pend.t.len * r + pend.n.t.minus.1*log(1 - exp(-pend.t.len * r))));
			}
			return(f)
		}
		
		int.t.0 <- partition[i.int,"t.0"];
		
		if (model == "bd") {
			f <- function(pars) {
				r <- pars[1];
				epsilon <- pars[2];
				if (r <= 0 | epsilon < 0 | epsilon >= 1) {return(-Inf);}
				l.int <- n.int * log(r) - r * sum.int.t.len - sum(log(1 - (epsilon * exp(-r * int.t.0))));
				ert <- exp(r * pend.t.len);
				B <- (ert - 1) / (ert - epsilon);
				l.pend <- sum(log(1 - B) + pend.n.t.minus.1*log(B));
				return(l.int + l.pend);
			}
			return(f)
		} else if (model == "fixedD") {
			fixD <- fixPar;
			f <- function(pars) {
				b <- pars[1];
				r <- b - fixD;
				epsilon <- fixD/b;
				if (r <= 0 | epsilon < 0 | epsilon >= 1) {return(-Inf);}
				l.int <- n.int * log(r) - r * sum.int.t.len - sum(log(1 - (epsilon * exp(-r * int.t.0))));
				ert <- exp(r * pend.t.len);
				B <- (ert - 1) / (ert - epsilon);
				l.pend <- sum(log(1 - B) + pend.n.t.minus.1*log(B));
				return(l.int + l.pend);
			}
			return(f)
		} else if (model == "fixedEpsilon") {
			epsilon <- fixPar;
			f <- function(pars) {
				r <- pars[1];
				if (r <= 0 | epsilon < 0 | epsilon >= 1) {return(-Inf);}
				l.int <- n.int * log(r) - r * sum.int.t.len - sum(log(1 - (epsilon * exp(-r * int.t.0))));
				ert <- exp(r * pend.t.len);
				B <- (ert - 1) / (ert - epsilon);
				l.pend <- sum(log(1 - B) + pend.n.t.minus.1*log(B));
				return(l.int + l.pend);
			}
			return(f)
		} else if (model == "fixedR") {
			r <- fixPar;
			f <- function(pars) {
				epsilon <- pars[1];
				if (r <= 0 | epsilon < 0 | epsilon >= 1) {return(-Inf);}
				l.int <- n.int * log(r) - r * sum.int.t.len - sum(log(1 - (epsilon * exp(-r * int.t.0))));
				ert <- exp(r * pend.t.len);
				B <- (ert - 1) / (ert - epsilon);
				l.pend <- sum(log(1 - B) + pend.n.t.minus.1*log(B));
				return(l.int + l.pend);
			}
			return(f)
		} else if (model == "fixedB") {
			fixB <- fixPar;
			f <- function(pars) {
				d <- pars[1];
				r <- fixB - d;
				epsilon <- d/fixB;
				if (r <= 0 | epsilon < 0 | epsilon >= 1) {return(-Inf);}
				l.int <- n.int * log(r) - r * sum.int.t.len - sum(log(1 - (epsilon * exp(-r * int.t.0))));
				ert <- exp(r * pend.t.len);
				B <- (ert - 1) / (ert - epsilon);
				l.pend <- sum(log(1 - B) + pend.n.t.minus.1*log(B));
				return(l.int + l.pend);
			}
			return(f)
		}
	}
}


## 'fit' contains '$par' and '$lnlik' and '$model'; check last value for presence of fixed parameters
calculateModelFit <- function (fit, z) {
## Sample size taken (for now) as the total num.nodes in the tree (internal + pendant)
  # num.nodes = (2*length(phy$tip.label) - 1) == (2*length(richness[,1]) - 1) == length(z[,1]) + 1
#	n <- (length(z[,1]) + 1) + sum(!is.na(z[,"n.f"]));
	
# Since each edge defines a node (i.e. an 'observation'), need only add root node as final obervation
	n <- (length(z[,1]) + 1);
	
 # Includes both formal parameters AND number of breaks. Note: first model does not involve a break.
## Models where all parameters are estimated (i.e. BD model):
  # 2 parameters for base model (no breakpoint) + 3 parameters (r, eps, breakpoint) for each subsequent model
  
# Determine number of piecewise models currently involved
	if (length(fit$par) < 3) { # i.e. base model
		num.models <- 1;
	} else {
		num.models <- length(fit$par[,1]);
	}
	model = fit$model[1];
	
	if (model == "fixedEpsilon" || model == "fixedR" || model == "fixedB" || model == "fixedD") {
		k <- 2 * num.models - 1;
	} else {
		k <- sum(!is.na(fit$par)) + num.models - 1; # number of estimated parameters + number of breaks
	}
	
	lnLik <- fit$lnLik;
	
	aic <- getAIC(lnLik, k);
	aicc <- getAICc(lnLik, k, n);
	
	model.fit <- c(aic, aicc, k);
	return(model.fit);
}


getAIC <- function (lnLik, k) {
	return(-2 * lnLik + 2*k);
}

getAICc <- function (lnLik, k, n) {
	return(getAIC(lnLik,k) + 2*k*(k+1)/(n-k-1));
}


## Used for comparing models fit to the same partition
getBestPartialModel <- function (fit1, fit2, z, criterion) {
	n <- (length(z[,1]) + 1); # the number of nodes involved
## Add '1' to parameters to account for break
	k1 <- 1 + sum(!is.na(fit1$par));
	k2 <- 1 + sum(!is.na(fit2$par));
	
	if (n - k1 <= 1 || n - k2 <= 1) { # deals with single edges, where AICc correction becomes undefined. use AIC.
		if (getAIC(fit1$lnLik,k1) < getAIC(fit2$lnLik,k2)) {
			return(fit1);
		} else {
			return(fit2);
		}
	} else {
		if (criterion == "aicc") {
			if (getAICc(fit1$lnLik, k1, n) < getAICc(fit2$lnLik, k2 ,n)) {
				return(fit1);
			} else {
				return(fit2);
			}
		} else {
			if (getAIC(fit1$lnLik, k1) < getAIC(fit2$lnLik, k2)) {
				return(fit1);
			} else {
				return(fit2);
			}
		}
	}
}


## Consider removing previously-fit rate shifts
backStep <- function (currentModel, z, step, model, fixPar, criterion) {
## As a first step, only consider removing entire shifts. Later deal with individual parameters.
	z.opt <- z;
	bestModel <- currentModel;
	bestScore <- as.numeric(bestModel[criterion]);
	allDeletedShifts <- NULL;
	bestRemoved <- NULL;
	improve <- T;
	
	while (improve) { # may be possible to remove > 1 previously fit shift
		allDeletedShifts <- c(allDeletedShifts, bestRemoved);
		currentModel <- bestModel;
		z <- z.opt;
		cuts <- bestModel$cut.at;
		nodes <- bestModel$split.at;
		pars <- bestModel$par;
		numModels <- length(bestModel$par)/2;
		improve <- F;
		
		if (numModels > 2) {
			for (i in 2:(numModels - 1)) { # don't waste time removing last shift
				fitModel <- currentModel;
				obj <- dissolveSplit(z, cut=cuts[i], node=nodes[i], aff=i);
				aff <- obj$affected;
				z.temp <- obj$z[obj$z[,"partition"] == aff,,drop=FALSE];
				
		## set par to weighted mean of 2 affected partitions
			## updated to reflect total path length rather than number of tips
			## really only influences weighted parameter starting values
			
				weights <- c(sum(z[which(z[,"partition"] == aff),"t.len"]), sum(z[which(z[,"partition"] == i),"t.len"]));
				sp <- c(weighted.mean(pars[c(aff,i),1], weights), weighted.mean(pars[c(aff,i),2], weights, na.rm=T));
				fit <- getOptimalModelFlavour(z=z.temp, sp=sp, model=model, fixPar=fixPar, criterion=criterion);
				
#				weights.old <- c(sum(z[,"partition"] == aff), sum(z[,"partition"] == i)); # weight from number of edges involved
#				sp.old <- c(weighted.mean(pars[c(aff,i),1], weights.old), weighted.mean(pars[c(aff,i),2], weights.old, na.rm=T));
#				fit.old <- getOptimalModelFlavour(z=z.temp, sp=sp.old, model=model, fixPar=fixPar, criterion=criterion);

		## Update fit values
				fitModel$par[aff,] <- fit$par;
				fitModel$par <- fitModel$par[-i,];
				fitModel$lnLik.part[aff] <- fit$lnLik;
				fitModel$lnLik.part <- fitModel$lnLik.part[-i];
				fitModel$lnLik <- sum(fitModel$lnLik.part);
				model.fit <- calculateModelFit(fit=fitModel, z=z);
				fitModel$aic <- model.fit[1];
				fitModel$aicc <- model.fit[2];
				fitModel$num.par <- model.fit[3];
				
				if (fitModel[criterion] < bestScore) {
					fitModel$split.at <- fitModel$split.at[-i];
					fitModel$model[aff] <- fit$model;
					fitModel$model <- fitModel$model[-i];
					fitModel$cut.at <- fitModel$cut.at[-i];
					bestModel <- fitModel;
					bestScore <- as.numeric(fitModel[criterion]);
					z.opt <- updateZ(z=obj$z, deletedPart=i);
					bestRemoved <- nodes[i];
					improve <- T;
				}
			}
			if (improve) {step <- rbind(step, c("remove", bestRemoved));}
		}
	}
	return(list(fit=bestModel, z=z.opt, step=step, remove=bestRemoved));
}


## Remove previously-fit rate shift
dissolveSplit <- function (z, cut, node, aff) {
## Grab ancestral branch partition membership
	anc <- z[which(z[,"dec"] == node)];
	root <- min(z[,"anc"]);
	tag <- NULL;
	
	if (cut == "node") {
		tag <- as.numeric(z[which(z[,"dec"] == node),"partition"]);
	} else if (cut == "stem" && anc > root) {
		tag <- as.numeric(z[which(z[,"dec"] == anc),"partition"]);
	} else if (cut == "stem" && anc == root) { # need to take other side of root
		dec <- z[which(z[,"anc"] == root),"dec"];
		tag <- as.numeric(z[which(z[,"dec"] == dec[which(dec != node)]),"partition"]); # ug. li.
	}
	
	idx <- which(z[,"partition"] == aff);
	z[idx,"partition"] <- tag;
	
	return(list(z=z, affected=tag));
}


## reduce partition IDs to reflect dissolved split
updateZ <- function (z, deletedPart) {
	idx <- z[,"partition"] > deletedPart;
	z[idx,"partition"] <- z[idx,"partition"] - 1;
	return(z);
}


## Only print if model improves AIC score
printRemovedShifts <- function (remove) {
	for (i in 1:length(remove)) {
		cat("  Removing shift at node #", remove[i], "\n", sep="");
	}
}