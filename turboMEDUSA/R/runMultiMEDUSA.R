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
			foo <- runTurboMEDUSA(phy[i], richness, model.limit, stop, criterion, initial.r, initial.e, plotFig, nexus, verbose, mc, num.cores, ...)
		}
	} else {
		res <- runTurboMEDUSA(phy, richness, model.limit, stop, criterion, initial.r, initial.e, plotFig, nexus, verbose, mc, num.cores, ...)
	}
	return(res)
}

