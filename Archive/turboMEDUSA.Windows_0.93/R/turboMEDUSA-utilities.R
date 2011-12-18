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
## shiftCut can only be "stem" or "node" (not "both"), as both are evaluated separately
medusa.ml.prefit <- function (node, z, desc, initialR, initialE, model, shiftCut, criterion)
{
	fitted.bd <- NULL;
	fitted.yule <- NULL;
#	fit <- NULL;
	
## if model is mixed, grab the optimally fitted one, drop the other; it is cutAtStem that matters, as it is the sum of 2 break likelihoods
  ## optimal alone may be different than optimal in tandem
	if (shiftCut == "stem")
	{
		obj <- medusa.split(node=node, z=z, desc=desc, shiftCut=shiftCut);
		if (model == "bd" || model == "mixed")
		{
			z.bd.stem <- obj$z;
			fitted.bd <- medusa.ml.fit.partition(partition=2, z=z.bd.stem, sp=c(initialR, initialE), model="bd");
			fitted.bd$model <- "bd";
			fitted.bd$cut.at <- "stem";
		}
		if (model == "yule" || model == "mixed")
		{
			z.yule.stem <- obj$z;
			fitted.yule <- medusa.ml.fit.partition(partition=2, z=z.yule.stem, sp=c(initialR, initialE), model="yule");
			fitted.yule$model <- "yule";
			fitted.yule$cut.at <- "stem";
		}
	} else if (shiftCut == "node")
	{
		obj <- medusa.split(node=node, z=z, desc=desc, shiftCut=shiftCut);
		if (model == "bd" || model == "mixed")
		{
			z.bd.node <- obj$z;
			fitted.bd <- medusa.ml.fit.partition(partition=2, z=z.bd.node, sp=c(initialR, initialE), model="bd");
			fitted.bd$model <- "bd";
			fitted.bd$cut.at <- "node";
		}
		if (model == "yule" || model == "mixed")
		{
			z.yule.node <- obj$z;
			fitted.yule <- medusa.ml.fit.partition(partition=2, z=z.yule.node, sp=c(initialR, initialE), model="yule");
			fitted.yule$model <- "yule";
			fitted.yule$cut.at <- "node";
		}
	}
## Check which flavour of model fits best
	if (is.null(fitted.bd) & !is.null(fitted.yule))
	{
		return(fitted.yule);
	} else if (is.null(fitted.yule) & !is.null(fitted.bd)) {
		return(fitted.bd);
	} else {
## Dealing with a 'mixed' model here; need to consider number of parameters.
		bd.model.fit <- calculate.model.fit(fit=fitted.bd, z=z);
		yule.model.fit <- calculate.model.fit(fit=fitted.yule, z=z);
		element <- NULL;
		if (criterion == "aic") {element <- 1;} else {element <- 2;}
		if (is.nan(yule.model.fit[[element]]))
		{
			return(fitted.bd);
		} else if (bd.model.fit[[element]] < yule.model.fit[[element]])
		{
			return(fitted.bd);
		} else {
			return(fitted.yule);
		}
	}
}






## sp = initializing values for r & epsilon
## Default values should never be used (except for first model), as the values from the previous model are passed in
medusa.ml.fit.partition <- function (partition, z, sp=c(0.05, 0.5), model)
{
	new.part <- z[z[,"partition"] == partition,,drop=FALSE];
# Check that partition is not empty; can occur with "node" or "both" cutting.
	if (length(new.part) == 0)
	{
		cat("\nSaved time, bitches!\n")
		par <- c(NA, NA);
		return(list(par=par, lnLik=-Inf));
	} else {
# Construct likelihood function:
		lik <- make.lik.medusa.part(partition=(z[z[,"partition"] == partition,,drop=FALSE]), model=model);
		foo <- function (x) {-lik(pars=exp(x));} # work with parameters in log-space to preserve precision
		
		if (model == "bd")
		{
			fit <- optim(fn=foo, par=log(sp), method="N", control=list(maxit=1000)); # last argument connotes maximization
			
			if(fit$convergence != 0) {cat("\nDidn't converge. Shit. Convergence = ", fit$convergence, "\n");}
			
			return(list(par=exp(fit$par), lnLik=-fit$value));
		} else {
			# node.richness <- sum(z[z[,"partition"] == partition,"n.t"], na.rm=TRUE);
			# sum.time <- sum(z[z[,"partition"] == partition,"t.len"], na.rm=TRUE);
			# best.guess <- (node.richness - 2 + (1 + 2^-50)) / (sum.time + (1 + 2^-50));
			# best.guess2 <- log(node.richness) / 
			
			# fit <- nlm(log(sp[1]), f=foo)
			# par <- c(exp(fit$estimate), NA);
			# return(list(par=par, lnLik=-fit$minimum));
			
			fit <- nlminb(log(sp[1]), objective=foo, upper=log(50), control=list(maxiter=1000))
			par <- c(exp(fit$par), NA);
			return(list(par=par, lnLik=-fit$objective));
	
			# fit <- optimize(f=foo, interval=c(-25, 1));
			# par <- c(exp(fit$minimum), NA);
			# return(list(par=par, lnLik=-fit$objective));
		}
	}
}






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
	base <- NULL;
	part <- z[,"partition"];
	
	if (shiftCut == "stem")
	{
		descendants <- desc$desc.stem;
	} else {
		descendants <- desc$desc.node;
	}
	
	base <- min(part[z[,1] == node | z[,2] == node]);
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
	op <- fit$par; # store previously fit parameter values
	cut.at <- NULL;
	cool <- TRUE;
	
	new.part.1 <- NULL;
	new.part.2 <- NULL;
	
	fit1 <- NULL;
	fit2 <- NULL;
	
	if (shiftCut == "stem" || shiftCut == "both" || node < root.node) # can enter on "node" if a tip
	{
## First, diminshed clade
		obj.stem <- medusa.split(node=node, z=z, desc=desc, shiftCut="stem");
		z.stem <- obj.stem$z;
		aff <- obj.stem$affected;

## Check that neither partition is not empty; can occur with "node" or "both" cutting. If so, kill it.
		new.part.1 <- z.stem[z.stem[,"partition"] == aff[1],,drop=FALSE];
		new.part.2 <- z.stem[z.stem[,"partition"] == aff[2],,drop=FALSE];
		
		if (length(new.part.1) == 0 || length(new.part.2) == 0)
		{
			fit$lnLik <- -Inf;
			fit$aic <- Inf;
			fit$aicc <- Inf;
			return(fit);
		}
		
## Everything is cool; proceed.
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
## vigin node, already calculated
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
	if ((shiftCut == "node" || shiftCut == "both") && !(node < root.node)) # never enter if tip
	{
## First, diminshed clade
		obj.node <- medusa.split(node=node, z=z, desc=desc, shiftCut="node");
		z.node <- obj.node$z;
		aff <- obj.node$affected;
		
## Need to check if cut is valid. May be inadmissable because of pattern of previous breaks (especially with shiftCut=both)
		if (is.na(aff[1]) || is.na(aff[2]))
		{
			cool <- FALSE;
		} else {
			new.part.1 <- z.node[z.node[,"partition"] == aff[1],,drop=FALSE];
			new.part.2 <- z.node[z.node[,"partition"] == aff[2],,drop=FALSE];
		}
		
## Check that partition is not empty; can occur with "node" or "both" cutting.
		if (length(new.part.1) == 0 || length(new.part.2) == 0 || !cool)
		{
			fit$lnLik <- -Inf;
			fit$aic <- Inf;
			fit$aicc <- Inf;
			return(fit);
		}
		
## Everything is cool; proceed.
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
## Second, new clade
		if (node < root.node) # tip, already calculated
		{
			fit2.node <- prefit$tips[[node]];
		} else if (length(unique(z.node[(z.node[,"partition"] == aff[2] & z.node[,"dec"] < root.node),"dec"])) == num.tips[[node]]) {
## vigin node, already calculated
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








## make.lik.medusa.part: generate a likelihood function for a single partition.
make.lik.medusa.part <- function (partition, model)
{
# Handle internal and pendant edges separately
	is.int <- is.na(partition[,"n.t"]);
	is.pend <- !is.int;
	
	n.int <- sum(is.int);
	n.pend <- sum(is.pend);
	
## Internal and pendant calculations differ; split'em up
	int  <- partition[is.int,,drop=FALSE];
	pend <- partition[is.pend,,drop=FALSE];
	
	sum.int.t.len <- sum(int[,"t.len"]);  # Simply sum all internal edges
	int.t.0 <- int[,"t.0"];
	
# 'n.0' = Foote's 'a', initial diversity; 'n.t' = Foote's 'n', final diversity
	pend.n.0 <- pend[,"n.0"]; # Foote's 'a': initial diversity
	pend.n.t <- pend[,"n.t"]; # Foote's 'n': final diversity
	pend.t.len <- pend[,"t.len"];
	
	if (n.int + n.pend != length(partition[,1])) stop("You messed up, yo.");
	
# User may pass in epsilon; don't change it, just estimate r
	f <- function(pars)
	{
		if (model == "bd")
		{
			r <- pars[1];
			epsilon <- pars[2];
			
			if (r <= 0 || epsilon <= 0 || epsilon >= 1) {return(-Inf);}
		} else {
			
			r <- pars[1];
			epsilon <- 0;
			
			if (r <= 0 || is.nan(r)) {return(-Inf);}
		}
		
		l.int <- numeric();
		l.pend <- numeric();
		
		if (n.int == 0) {l.int <- 0;} else {
## Likelihood of internal edges from Rabosky et al. (2007) equation (2.3):
			l.int <- n.int * log(r) - r * sum.int.t.len - sum(log(1 - (epsilon * exp(-r * int.t.0))));
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
			
			ert <- exp(r * pend.t.len);
			B <- (ert - 1) / (ert - epsilon); # Equivalently: B <- (bert - b) / (bert - d)
			
			l.pend <- sum(log(1 - B) + (pend.n.t - 1)*log(B));
		}
		return(l.int + l.pend);
	}
}


## 'fit' contains '$par' and '$lnlik'
calculate.model.fit <- function (fit, z)
{
## Sample size taken (for now) as the total num.nodes in the tree (internal + pendant)
  # num.nodes = (2*length(phy$tip.label) - 1) == (2*length(richness[,1]) - 1) == length(z[,1]) + 1
#	n <- (length(z[,1]) + 1) + sum(!is.na(z[,"n.f"]));
	
# Since each edge defines a node (i.e. an 'observation'), need only add root node as final obervation
	n <- (length(z[,1]) + 1);
	
 # Includes both formal parameters AND number of breaks. Note: first model does not involve a break.
## Models where all parameters are estimated (i.e. BD model):
  # 2 parameters for base model (no breakpoint) + 3 parameters (r, eps, breakpoint) for each subsequent model
  
  
# Determine number of piecewise models currently involved
	if (length(fit$par) < 3) # i.e. base model
	{
		num.models <- 1;
	} else {
		num.models <- length(fit$par[,1]);
	}
	
# Updated for more general models: check how many parameter values != NA
#	k <- 2 + (3 * (num.models - 1))
	k <- sum(!is.na(fit$par)) + (num.models - 1); # number of estimated parameters + number of breaks
	
	lnLik <- fit$lnLik;
	
	aic <- (-2 * lnLik) + (2*k);
	aicc <- aic + 2*k*(k+1)/(n-k-1);
	
	model.fit <- c(aic, aicc, k);
	return(model.fit);
}


## Get b and d values from r (b-d) and epsilson (d/b)
## Used in previous version of program; now in terms of r and epsilon
## Possibly of use to users wishing to translate results
get.b.d <- function (r, epsilon)
{
	b <- r/(1-epsilon);
	d <- b-r;   # Alternatively: d <- eps*r/(1-eps)
	return(list(b=b, d=d));
}

## Print out tree with ape-style node-numbering
## Possibly of interest for users to identify numbers of node(s) off interest
 ## If this is the case, make sure to pass in pruned tree
plotNN <- function (phy, time=TRUE, margin=TRUE, label.offset=0.5, cex=0.5, ...) 
{
	phy$node.label <- (length(phy$tip.label) + 1):max(phy$edge);
	plot.phylo(phy, show.node.label=TRUE, no.margin=!margin, label.offset=label.offset, cex=cex, ...);
	if (time && !margin) {cat("Cannot plot time axis without a margin.\n");}
	else if (time && margin) {axisPhylo(cex.axis=0.75)};
}