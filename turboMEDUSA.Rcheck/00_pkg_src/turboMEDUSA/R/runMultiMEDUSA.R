runMultiMEDUSA <-
function(phy, richness=NULL, model.limit=20, stop="model.limit", model="bd",
	criterion="aicc", shiftCut="stem", initialR=0.05, initialE=0.5, plotFig=FALSE, nexus=FALSE, verbose=FALSE, mc=FALSE, num.cores=NULL, ...)
{
	if (nexus) phy <- read.nexus(phy)
	res <- list()
	if (class(phy) == "multiPhylo")
	{
		num.trees <- length(phy)
		for (i in 1:num.trees)
		{
			cat("\nProcessing tree ", i, " (of ", num.trees, ")...\n\n", sep="")
			res[[i]] <- runTurboMEDUSA(phy=phy[[i]], richness=richness, model.limit=model.limit, stop=stop, model=model, criterion=criterion, shiftCut=shiftCut, initialR=initialR, initialE=initialE, plotFig=plotFig, nexus=nexus, verbose=verbose, mc=mc, num.cores=num.cores, ...)
		}
	} else {
		res <- runTurboMEDUSA(phy, richness=richness, model.limit=model.limit, stop=stop, model=model, criterion=criterion, shiftCut=shiftCut, initialR=initialR, initialE=initialE, plotFig=plotFig, nexus=nexus, verbose=verbose, mc=mc, num.cores=num.cores, ...)
	}
	return(res)
}

