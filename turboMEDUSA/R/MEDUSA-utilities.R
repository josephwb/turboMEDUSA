## Only used for base model
medusaMLFitBase <- function (z, sp, model, fixPar, criterion)
{
	fit <- getOptimalModelFlavour(partition.id=1, z=z, sp=sp, model=model, fixPar=fixPar, criterion=criterion);
	model.fit <- calculateModelFit(fit=fit, z=z);
	
	return(list(par=matrix(fit$par, nrow=1, dimnames=list(NULL,c("r", "epsilon"))), lnLik.part=fit$lnLik, lnLik=fit$lnLik,
		split.at=min(z[,"anc"]), aic=round(model.fit[1], digits=7), aicc=round(model.fit[2], digits=7), num.par=model.fit[3],
		cut.at="node", model=fit$model));
}


## Split the edge matrix 'z' by adding a partition at node 'node'.
## The list 'desc' is a list of descendants (see makeCacheMedusa, above).
## Returns a list with elements:
##   z: new medusa matrix, with the new partition added
##   affected: indices of the partitions affected by the split (n == 2).
medusaSplitStem <- function (node, z, desc)
{
	part <- z[,"partition"];
	base <- min(part[z[,"anc"] == node | z[,"dec"] == node]);
	tag <- max(part) + 1;
	i <- desc[[node]];
	idx <- i[part[i] == base];
	z[idx,"partition"] <- tag;
	z[which(z["dec"] == node),"partition"] <- tag; # Possible to have several edges to consider
	return(list(z=z, affected=c(unique(part[idx]), tag)));
}

medusaSplitNode <- function (node, z, desc)
{
	part <- z[,"partition"];
	base <- min(part[z[,"anc"] == node]);
	tag <- max(part) + 1;
	i <- desc[[node]];
	idx <- i[part[i] == base];
	z[idx,"partition"] <- tag;
	z[which(z["dec"] == node),"partition"] <- tag; # Possible to have several edges to consider
	return(list(z=z, affected=c(unique(part[idx]), tag)));
}

medusaSplit <- function (node, z, desc, shiftCut)
{
	descendants <- desc$desc.stem; # will be correct half the time
	if (shiftCut == "node") {descendants <- desc$desc.node;}
	
	part <- z[,"partition"];
	base <- min(part[z[,1] == node | z[,2] == node]);
	tag <- max(part) + 1;
	i <- descendants[[node]];
	idx <- i[part[i] == base];
	z[idx,"partition"] <- tag;
	z[which(z["dec"] == node),"partition"] <- tag; # Possible to have several edges to consider
	return(list(z=z, affected=c(unique(part[idx]), tag)));
}


## Pre-fit values for pendant edges; DON'T recalculate later; should account for ~25% of all calculations
## Also cache values for virgin nodes; useful until subsetted.
## shiftCut can only be "stem" or "node" (not "both"), as both are evaluated separately
medusaMLPrefitNode <- function (node, z, desc, sp, model, fixPar, criterion)
{
	z.node <- medusaSplitNode(node=node, z=z, desc=desc)$z;
	fit <- getOptimalModelFlavour(partition.id=2, z=z.node, sp=sp, model=model, fixPar=fixPar, criterion=criterion);
	return(fit);
}


## Pre-fit values for pendant edges; DON'T recalculate later; should account for ~25% of all calculations
## Also cache values for virgin nodes; useful until subsetted.
## shiftCut can only be "stem" or "node" (not "both"), as both are evaluated separately
medusaMLPrefitStem <- function (node, z, desc, sp, model, fixPar, criterion)
{
	z.stem <- medusaSplitStem(node=node, z=z, desc=desc)$z;
	fit <- getOptimalModelFlavour(partition.id=2, z=z.stem, sp=sp, model=model, fixPar=fixPar, criterion=criterion);
	return(fit);
}


## When model == mixed, fit both and find optimal flavour
getOptimalModelFlavour <- function (partition.id, z, sp, model, fixPar, criterion)
{
	fit.bd <- NULL;
	fit.yule <- NULL;
	fit <- NULL;
	
	if (model == "yule" | model == "mixed")
	{
		fit.yule <- medusaMLFitPartition(partition.id, z, sp=sp, model="yule");
		fit.yule$model <- "yule";
	}
	if (model == "bd" | model == "mixed")
	{
		if (is.na(sp[2])) {sp[2] <- 0.5;}
		fit.bd <- medusaMLFitPartition(partition.id, z, sp=sp, model="bd");
		fit.bd$model <- "bd";
	}
	if (model != "mixed" && model != "bd" && model != "yule")
	{
		fit <- medusaMLFitPartition(partition.id, z, sp=sp, model=model, fixPar=fixPar);
		fit$model <- model;
		return(fit);
	}
	
## Figure out which model fits best
	if (is.null(fit.bd))
	{
		fit <- fit.yule;
	} else if (is.null(fit.yule)) {
		fit <- fit.bd;
	} else {
## Considering both places for a shift
		fit.bd.val <- calculateModelFit(fit=fit.bd, z=z);
		fit.yule.val <- calculateModelFit(fit=fit.yule, z=z);
		
		if (criterion == "aic") {element <- 1;} else {element <- 2;}
		
		if (fit.bd.val[[element]] < fit.yule.val[[element]])
		{
			fit <- fit.bd;
		} else {
			fit <- fit.yule;
		}
	}
	return(fit);
}


## 'fit' contains parameter values from previous model, used to initialize subsequent model.
## Pass in pre-fitted values for pendant edges and virgin nodes (in 'prefit'); DON'T recalculate.
## Need to consider the possibility of birth-death, yule, or mixed models.
## Need to consider where shft is placed (shiftCut). Placement affects not only new clade, but
## also the size of the split clade. Only relevant if shiftCut = "both".
## fit1 model is already logged; only need to record fit2 model, and only non-prefitted nodes

# *** Maybe just call the two functions above? ***
medusaMLUpdate <- function (node, z, desc, fit, prefit, num.tips, root.node, model, fixPar, criterion, shiftCut, preserveModelFlavour)
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
	
	if (shiftCut == "stem" | shiftCut == "both" | node < root.node) # can enter on "node" if a tip
	{
## First, diminshed clade
		obj.stem <- medusaSplitStem(node=node, z=z, desc=desc$desc.stem);
		z.stem <- obj.stem$z;
		aff <- obj.stem$affected;

## Check that neither partition is not empty; can occur with "node" or "both" cutting. If so, kill it.
		if (sum(z.stem[,"partition"] == aff[1]) == 0 | sum(z.stem[,"partition"] == aff[2]) == 0)
		{
			fit$lnLik <- -Inf;
			fit$aic <- Inf;
			fit$aicc <- Inf;
			return(fit);
		}
		
## Everything is cool; proceed.
		sp <- op[aff[1],]; # Use previously fit parameter values from clade that is currently being split
		
		if (model == "mixed")
		{
			if (preserveModelFlavour)  ## In mixed models, may want to conserve flavour of previously fit model
			{
				if (sum(!is.na(sp)) < 2) # yule
				{
					fit1.stem <- medusaMLFitPartition(partition.id=aff[1], z=z.stem, sp=sp, model="yule");
				} else {
					fit1.stem <- medusaMLFitPartition(partition.id=aff[1], z=z.stem, sp=sp, model="bd");
				}
			} else {
## consider both model flavours
				fit1.stem <- getOptimalModelFlavour(partition.id=aff[1], z=z.stem, sp=sp, model=model, fixPar=fixPar, criterion=criterion);
			}
		} else {
			fit1.stem <- medusaMLFitPartition(partition.id=aff[1], z=z.stem, sp=sp, model=model, fixPar=fixPar);
			fit1.stem$model <- model;
		}
## Second, new clade
		if (node < root.node) # tip, already calculated
		{
			fit2.stem <- prefit$tips[[node]];
		} else if (length(unique(z.stem[(z.stem[,"partition"] == aff[2] & z.stem[,"dec"] < root.node),"dec"])) == prefit$num.tips[[node]]) {
## vigin node, already calculated
			fit2.stem <- prefit$virgin.nodes$stem[[node - root.node]];
		} else {
## novel shift
			fit2.stem <- getOptimalModelFlavour (partition.id=aff[2], z=z.stem, sp=sp, model=model, fixPar=fixPar, criterion=criterion);
		}
	}
	if ((shiftCut == "node" || shiftCut == "both") && (node > root.node)) # never enter if tip
	{
## First, diminshed clade
		obj.node <- medusaSplitNode(node=node, z=z, desc=desc$desc.node);
		z.node <- obj.node$z;
		aff <- obj.node$affected;
		
## Need to check if cut is valid. May be inadmissable because of pattern of previous breaks (especially with shiftCut=both)
		if (is.na(aff[1]) | is.na(aff[2]))
		{
			cool <- FALSE;
		}		
## Check that partition is not empty; can occur with "node" or "both" cutting.
#		if (length(new.part.1) == 0 | length(new.part.2) == 0 | !cool)
		if (sum(z.node[,"partition"] == aff[1]) == 0 | sum(z.node[,"partition"] == aff[2]) == 0 | !cool)
		{
			fit$lnLik <- -Inf;
			fit$aic <- Inf;
			fit$aicc <- Inf;
			return(fit);
		}
		
## Everything is cool; proceed.
		sp <- op[aff[1],]; # Use previously fit parameter values from clade that is currently being split
		
		if (model == "mixed")
		{
			if (preserveModelFlavour)  ## In mixed models, may want to conserve flavour of previously fit model
			{
				if (sum(!is.na(sp)) < 2) # yule
				{
					fit1.node <- medusaMLFitPartition(partition.id=aff[1], z=z.node, sp=sp, model="yule"); # should this change? probably.
				} else {
					fit1.node <- medusaMLFitPartition(partition.id=aff[1], z=z.node, sp=sp, model="bd");
				}
			} else {
## consider both model flavours
				fit1.node <- getOptimalModelFlavour(partition.id=aff[1], z=z.node, sp=sp, model=model, fixPar=fixPar, criterion=criterion);
			}
		} else {
			fit1.node <- medusaMLFitPartition(partition.id=aff[1], z=z.node, sp=sp, model=model, fixPar=fixPar);
			fit1.node$model <- model;
		}
## Second, new clade
		if (length(unique(z.node[(z.node[,"partition"] == aff[2] & z.node[,"dec"] < root.node),"dec"])) == prefit$num.tips[[node]]) {
## vigin node, already calculated
			fit2.node <- prefit$virgin.nodes$node[[node - root.node]];
		} else {
## novel shift
			fit2.node <- getOptimalModelFlavour (partition.id=aff[2], z=z.node, sp=sp, model=model, fixPar=fixPar, criterion=criterion);
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
		stem.val <- list(lnLik=stem.lik, par=stem.par, model=model);
		stem.fit <- calculateModelFit(fit=stem.val, z=z);
		
		node.lik <- (fit1.node$lnLik + fit2.node$lnLik);
		node.par <- rbind(fit1.node$par, fit2.node$par)
		node.val <- list(lnLik=node.lik, par=node.par, model=model);
		node.fit <- calculateModelFit(fit=node.val, z=z);
		
		if (criterion == "aic") {element <- 1;} else {element <- 2;}
		
		if (stem.fit[[element]] < node.fit[[element]])
		{
			fit1 <- fit1.stem; # is this correct
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


## sp = initializing values for r & epsilon
## Default values should never be used (except for first model), as the values from the previous model are passed in
medusaMLFitPartition <- function (partition.id, z, sp=c(0.05, 0.5), model, fixPar=NULL)
{
	# SHOULD I PUT medusaSplit HERE? WOULD HAVE TO PASS IN Aff, shiftCut, model, node, OTHERS?!?
	# WOULD PUT ALL(?) MODEL COMPARISONS HERE
	# HUGE medusaMLUpdate WOULD HAVE TO GO HERE, BUT AVOIDS 10 CALLS TO THIS FUNCTION; 19 IN ALL.
	
	new.part <- z[z[,"partition"] == partition.id,,drop=FALSE];
# Check that partition is not empty; can occur with "node" or "both" cutting.
	if (length(new.part) == 0 || (length(new.part[,"partition"] == 1) && new.part[,"t.len"] == 1)) ## WHAT WAS I TRYING TO DO HERE WITH (new.part[,"t.len"] == 1)?!? SHOULD THIS BE (new.part[,"t.len"] == 0)???
	{
		cat("\nWhoops! Something went wrong...\n")
		par <- c(NA, NA);
		return(list(par=par, lnLik=-Inf));
	}
# Construct likelihood function:
	lik <- makePartitionLikelihood(partition=new.part, model=model, fixPar=fixPar);
	foo <- function (x) {-lik(pars=exp(x));} # work with parameters in log-space to preserve precision
	
	if (model == "bd")
	{
		fit <- optim(fn=foo, par=log(sp), method="N", control=list(maxit=5000));
		
		if(fit$convergence != 0) {cat("\nDidn't converge. Shit. Convergence = ", fit$convergence, "\n");}
		
		return(list(par=exp(fit$par), lnLik=-fit$value));
		
	} else if (model == "yule") {
## OLD: fit2 <- optimize(f=lik, interval=c(0, 1), maximum=TRUE)
		
		node.richness  <- sum(new.part[,"n.t"], na.rm=TRUE);
		depth <- max(new.part[,"t.0"]);
		maxVal <- (log(node.richness) / depth) * 3;
		
		if (node.richness == 1)
		{
			maxVal <- 1
		}
		
		# node.richness <- sum(new.part[,"n.t"], na.rm=TRUE);
		# sum.time <- sum(new.part[,"t.len"]);
		# best.guess <- (node.richness - 2 + (1 + 2^-50)) / (sum.time + (1 + 2^-50));
		# best.guess2 <- log(node.richness) / sum.time
		
		# fit <- nlm(log(sp[1]), f=foo, iterlim=5000)
		# par <- c(exp(fit$estimate), NA);
		# return(list(par=par, lnLik=-fit$minimum));
		
		# fit <- nlminb(log(sp[1]), objective=foo, lower=1e-50, upper=log(maxVal), control=list(maxiter=1000))
		# par <- c(exp(fit$par), NA);
		# return(list(par=par, lnLik=-fit$objective));
		
# use different intervals based on model flavour
		
#		fit <- optimize(f=foo, interval=c(-25, log(maxVal)));
	
		fit <- optimize(f=foo, interval=c(-25, 1));
		par <- c(exp(fit$minimum), NA);
	
		return(list(par=par, lnLik=-fit$objective));
	
	} else if (model == "fixedD") {
		
		fit <- optimize(f=foo, interval=c(log(fixPar), 2.5));
		b <- exp(fit$minimum);
		par <- c((b - fixPar), (fixPar / b));
		return(list(par=par, lnLik=-fit$objective));
		
	} else if (model == "fixedEpsilon") {
		
		fit <- optimize(f=foo, interval=c(-25, 0));
		par <- c(exp(fit$minimum), sp[2]);
		return(list(par=par, lnLik=-fit$objective));
		
	} else if (model == "fixedB") {
		
		fit <- optimize(f=foo, interval=c(-50, log(fixPar)));
		d <- exp(fit$minimum);
		par <- c((fixPar - d), (d / fixPar));
		return(list(par=par, lnLik=-fit$objective));
		
	} else if (model == "fixedR") {
		
		fit <- optimize(f=foo, interval=c(-50, 0));
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
makePartitionLikelihood <- function (partition, model, fixPar)
{
	i.int <- is.na(partition[,"n.t"])
	i.pend <- !(i.int)
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
				if (r <= 0 | epsilon <= 0 | epsilon >= 1) {return(-Inf);}
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
				if (r <= 0 | epsilon <= 0 | epsilon >= 1) {return(-Inf);}
				l.int <- n.int * log(r) - r * sum.int.t.len - sum(log(1 - (epsilon * exp(-r * int.t.0))));
				l.pend <- sum(log(1 - epsilon) - log(exp(r * pend.t.len) - epsilon));
				return(l.int + l.pend);
			}
			return(f)
		} else if (model == "fixedEpsilon") {
			epsilon <- fixPar;
			f <- function(pars) {
				r <- pars[1];
				if (r <= 0 | epsilon <= 0 | epsilon >= 1) {return(-Inf);}
				l.int <- n.int * log(r) - r * sum.int.t.len - sum(log(1 - (epsilon * exp(-r * int.t.0))));
				l.pend <- sum(log(1 - epsilon) - log(exp(r * pend.t.len) - epsilon));
				return(l.int + l.pend);
			}
			return(f)
		} else if (model == "fixedR") {
			r <- fixPar;
			f <- function(pars) {
				epsilon <- pars[1];
				if (r <= 0 | epsilon <= 0 | epsilon >= 1) {return(-Inf);}
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
				if (r <= 0 | epsilon <= 0 | epsilon >= 1) {return(-Inf);}
				l.int <- n.int * log(r) - r * sum.int.t.len - sum(log(1 - (epsilon * exp(-r * int.t.0))));
				l.pend <- sum(log(1 - epsilon) - log(exp(r * pend.t.len) - epsilon));
				return(l.int + l.pend);
			}
			return(f)
		}
	} else {
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
				if (r <= 0 | epsilon <= 0 | epsilon >= 1) {return(-Inf);}
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
				if (r <= 0 | epsilon <= 0 | epsilon >= 1) {return(-Inf);}
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
				if (r <= 0 | epsilon <= 0 | epsilon >= 1) {return(-Inf);}
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
				if (r <= 0 | epsilon <= 0 | epsilon >= 1) {return(-Inf);}
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
				if (r <= 0 | epsilon <= 0 | epsilon >= 1) {return(-Inf);}
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
calculateModelFit <- function (fit, z)
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
	model = fit$model[1];
	
	if (model == "fixedEpsilon" || model == "fixedR" || model == "fixedB" || model == "fixedD")
	{
		k <- 2 * num.models - 1;
	} else {
		k <- sum(!is.na(fit$par)) + num.models - 1; # number of estimated parameters + number of breaks
	}
	
	lnLik <- fit$lnLik;
	
	aic <- -2 * lnLik + 2*k;
	aicc <- aic + 2*k*(k+1)/(n-k-1);
	
	model.fit <- c(aic, aicc, k);
	return(model.fit);
}