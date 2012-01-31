multiMEDUSA <-
function(phy, richness=NULL, modelLimit=20, stop="modelLimit", model="bd",
	criterion="aicc", shiftCut="stem", initialR=0.05, initialE=0.5, plotFig=FALSE, nexus=FALSE, verbose=FALSE, mc=FALSE, numCores=NULL, ...)
{
	if (nexus) phy <- read.nexus(phy)
	res <- list()
	if (class(phy) == "multiPhylo")
	{
		num.trees <- length(phy)
		for (i in 1:num.trees)
		{
			cat("\nProcessing tree ", i, " (of ", num.trees, ")...\n\n", sep="")
			res[[i]] <- MEDUSA(phy=phy[[i]], richness=richness, modelLimit=modelLimit, stop=stop, model=model, criterion=criterion, shiftCut=shiftCut, initialR=initialR, initialE=initialE, plotFig=plotFig, nexus=nexus, verbose=verbose, mc=mc, numCores=numCores, ...)
		}
	} else {
		res <- MEDUSA(phy, richness=richness, modelLimit=modelLimit, stop=stop, model=model, criterion=criterion, shiftCut=shiftCut, initialR=initialR, initialE=initialE, plotFig=plotFig, nexus=nexus, verbose=verbose, mc=mc, numCores=numCores, ...)
	}
	return(res)
}

