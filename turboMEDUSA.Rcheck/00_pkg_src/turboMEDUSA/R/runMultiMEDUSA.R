runMultiMEDUSA <-
function(phy, richness=NULL, model.limit=20, stop="model.limit",
	criterion="aicc", initial.r=0.05, initial.e=0.5, plotFig=FALSE, nexus=FALSE, verbose=FALSE, mc=FALSE, num.cores=NULL, ...)
{
	if (nexus) phy <- read.nexus(phy)
	res <- list()
	if (class(phy) == "multiPhylo")
	{
		num.trees <- length(phy)
		for (i in 1:num.trees)
		{
			cat("Processing tree ", i, " (of ", num.trees, ")...\n\n", sep="")
			res[[i]] <- runTurboMEDUSA(phy=phy[[i]], richness=richness, model.limit=model.limit, stop=stop, criterion=criterion, initial.r=initial.r, initial.e=initial.e, plotFig=plotFig, nexus=nexus, verbose=verbose, mc=mc, num.cores=num.cores, ...)
		}
	} else {
		res <- runTurboMEDUSA(phy, richness=richness, model.limit=model.limit, stop=stop, criterion=criterion, initial.r=initial.r, initial.e=initial.e, plotFig=plotFig, nexus=nexus, verbose=verbose, mc=mc, num.cores=num.cores, ...)
	}
	return(res)
}

