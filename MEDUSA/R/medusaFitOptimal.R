## re-fit optimal node shift to get all required values; significant memory issue for large trees if done otherwise.
## shiftCut and node are fixed. model is as well, but it seems tricky...
medusaFitOptimal <- function (node, z, desc, fit, prefit, root.node, model, fixPar, criterion, shiftCut, preserveModelFlavour)
{
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
	
	if (shiftCut == "stem")
	{
## First, diminshed clade
		obj <- medusaSplitStem(node=node, z=z, desc=desc$stem);
		z <- obj$z;
		aff <- obj$affected;

## Ensure that neither partition is empty; can occur with "node" or "both" cutting. If so, kill it.
		if (sum(z[,"partition"] == aff[1]) == 0 || sum(z[,"partition"] == aff[2]) == 0)
		{
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
		if (length(x) == 1) # a cut-at-stem scenario
		{
			y <- as.numeric(dimClade[x, "dec"]);
			if (length(unique(dimClade[(dimClade[,"dec"] < root.node),"dec"])) == prefit$num.tips[[y]])
			{
				if (y < root.node)
				{
					fit1 <- prefit$tips[[y]];
				} else {
					fit1 <- prefit$virgin.nodes$stem[[y - root.node]];
				}
			}
		}
				
		if (is.null(fit1))
		{
			if (model == "mixed")
			{
				if (preserveModelFlavour)  ## In mixed models, may want to conserve flavour of previously fit model
				{
					if (sum(!is.na(sp)) < 2) # yule
					{
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
		if (node < root.node) # tip, already calculated
		{
			fit2 <- prefit$tips[[node]];
		} else if (length(unique(z[(z[,"partition"] == aff[2] & z[,"dec"] < root.node),"dec"])) == prefit$num.tips[[node]]) {
## vigin node, already calculated
			fit2 <- prefit$virgin.nodes$stem[[node - root.node]];
		} else {
## novel shift
			newClade <- z[z[,"partition"] == aff[2],,drop=FALSE];
			fit2 <- getOptimalModelFlavour(z=newClade, sp=sp, model=model, fixPar=fixPar, criterion=criterion);
		}
	} else if (shiftCut == "node") # never enter if tip
	{
## First, diminshed clade
		obj <- medusaSplitNode(node=node, z=z, desc=desc$node);
		z <- obj$z;
		aff <- obj$affected;
		
## Need to check if cut is valid. May be inadmissable because of pattern of previous breaks
		if (is.na(aff[1]) || is.na(aff[2]))
		{
			cool <- FALSE;
		}		
## Check that partition is not empty; can occur with "node" cutting.
		if (sum(z[,"partition"] == aff[1]) == 0 || sum(z[,"partition"] == aff[2]) == 0 || !cool)
		{
			fit$lnLik <- -Inf;
			fit$aic <- Inf;
			fit$aicc <- Inf;
			return(fit);
		}
		
## Everything is cool; proceed.
		sp <- op[aff[1],]; # Use previously fit parameter values from clade that is currently being split
		
## first, consider diminshed clade. may result in a clade that has been cached
		dimClade <- z[z[,"partition"] == aff[1],,drop=FALSE];
		
		if (model == "mixed")
		{
			if (preserveModelFlavour)  ## In mixed models, may want to conserve flavour of previously fit model
			{
				if (sum(!is.na(sp)) < 2) # yule
				{
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