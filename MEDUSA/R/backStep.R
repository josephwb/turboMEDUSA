backStep <- function (currentModel, z, step, model, fixPar, criterion)
{
## As a first step, only consider removing entire shifts. Later deal with individual parameters.
	
	z.opt <- z;
	bestModel <- currentModel;
	bestScore <- as.numeric(bestModel[criterion]);
	allDeletedShifts <- NULL;
	bestRemoved <- NULL;
	improve <- T;
	
	while (improve) # may be possible to remove > 1 previously fit shift
	{
		allDeletedShifts <- c(allDeletedShifts, bestRemoved);
		currentModel <- bestModel;
		z <- z.opt;
		cuts <- bestModel$cut.at;
		nodes <- bestModel$split.at;
		pars <- bestModel$par;
		numModels <- length(bestModel$par)/2;
		improve <- F;
		
		if (numModels > 2)
		{
			for (i in 2:(numModels - 1)) # don't waste time removing last shift
			{
				fitModel <- currentModel;
				obj <- dissolveSplit(z, cut=cuts[i], node=nodes[i], aff=i);
				aff <- obj$affected;
				z.temp <- obj$z[obj$z[,"partition"] == aff,,drop=FALSE];
				
			# set par to mean of 2 affected partitions. perhaps weight (where weight comes form number of edges)
				sp <- c(mean(pars[c(aff,i),1]), mean(pars[c(aff,i),2], na.rm=T));
				
				fit <- getOptimalModelFlavour(z=z.temp, sp=sp, model=model, fixPar=fixPar, criterion=criterion);
				
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
				
				if (fitModel[criterion] < bestScore)
				{
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


dissolveSplit <- function (z, cut, node, aff)
{
## Grab ancestral branch partition membership
	anc <- z[which(z[,"dec"] == node)];
	root <- min(z[,"anc"]);
	tag <- NULL;
	
	if (cut == "node")
	{
		tag <- as.numeric(z[which(z[,"dec"] == node),"partition"]);
	} else if (cut == "stem" && anc > root)
	{
		tag <- as.numeric(z[which(z[,"dec"] == anc),"partition"]);
	} else if (cut == "stem" && anc == root) { # need to take other side of root
		dec <- z[which(z[,"anc"] == root),"dec"];
		tag <- as.numeric(z[which(z[,"dec"] == dec[which(dec != node)]),"partition"]); # ug. li.
	}
	
	idx <- which(z[,"partition"] == aff);
	z[idx,"partition"] <- tag;
	
	return(list(z=z, affected=tag));
}


updateZ <- function (z, deletedPart)
{
	idx <- z[,"partition"] > deletedPart;
	z[idx,"partition"] <- z[idx,"partition"] - 1;
	return(z);
}


## Only print if model improves AIC score
printRemovedShifts <- function (remove)
{
	for (i in 1:length(remove))
	{
		cat("  Removing shift at node #", remove[i], "\n", sep="");
	}
}