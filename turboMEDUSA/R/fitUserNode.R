medusaFitClade <- function (phy=phy, richness=NULL, userNode=NULL, model="mixed",
	fixedEpsilon=NULL, criterion="aicc", shiftCut="both", extractClade=FALSE, initialR=0.05, initialE=0.5,
	plotFig=FALSE, verbose=TRUE, mc=FALSE, numCores=NULL, ...)
{
# Used for testing specific hypotheses, whether between clades, or against background
	phyData <- prepareData(phy=phy, richness=richness, verbose=verbose);
	phy <- phyData$phy;
	richness <- phyData$richness;
	
## Keep track of all nodes, internal and pendant (for keeping track of breakpoints)
	pend.nodes <- seq_len(length(phy$tip.label));   # Calculate pendant splits just once, keep track through various models
	int.nodes <- (length(phy$tip.label)+2):max(phy$edge); # Omit root node
	root.node <- length(phy$tip.label) + 1;
	all.nodes <- c(pend.nodes, root.node, int.nodes);
	
	## The important bits. Set up z, get descendants and number of tips per node
	obj <- makeCacheMedusa(phy=phy, richness=richness, all.nodes=all.nodes, mc=mc, numCores=numCores);
	
	desc <- list(desc.stem=obj$desc.stem, desc.node=obj$desc.node);
	z <- obj$z;
	z.orig <- z; # Save for summarizing models later

	fit <- medusaMLBaseModel(z=z, extractClade=extractClade, nodes=userNode, model=model,
		shiftCut=shiftCut, desc=desc, criterion=criterion, initialR=initialR, initialE=initialE);
	
	models <- list(fit);	
	
	for (i in 1:length(userNode))
	{
		fit <- list(medusaFitUserNode(userNode[i], z, desc, model, root.node, criterion, shiftCut, initialR, initialE));
		models <- c(models, fit);
	}
	results <- list(z=z.orig, desc=desc, models=models, phy=phy);
	
# *** Need new summary function ***
	
	
	
	
	
	
	
#	return(models);
}

medusaFitUserNode <- function (node, z, desc, model, root.node, criterion, shiftCut, initialR, initialE)
{
## various combinations possible
	fit.stem <- NULL;
	fit.node <- NULL;
	cut.at <- NULL;
	
	if (shiftCut == "stem" | shiftCut == "both")
	{
		obj.stem <- medusaSplitStem(node=node, z=z, desc=desc$desc.stem);
		z.stem <- obj.stem$z;
		
		fit.stem.bd <- NULL;
		fit.stem.yule <- NULL;
		
		if (model == "yule" | model == "mixed")
		{
			fit.stem.yule <- medusaMLFitPartition(partition.id=2, z.stem, sp=c(initialR, initialE), model="yule");
			fit.stem.yule$model <- "yule";
		}
		if (model == "bd" | model == "mixed")
		{
			fit.stem.bd <- medusaMLFitPartition(partition.id=2, z.stem, sp=c(initialR, initialE), model="bd");
			fit.stem.bd$model <- "bd";
		}
## Figure out which model fits best
		if (is.null(fit.stem.bd))
		{
			fit.stem <- fit.stem.yule;
		} else if (is.null(fit.stem.yule)) {
			fit.stem <- fit.stem.bd;
		} else {
## Considering both places for a shift
			fit.stem.bd.val <- calculateModelFit(fit=fit.stem.bd, z=z);
			fit.stem.yule.val <- calculateModelFit(fit=fit.stem.yule, z=z);
			if (criterion == "aic") {element <- 1;} else {element <- 2;}
			if (fit.stem.bd.val[[element]] < fit.stem.yule.val[[element]])
			{
				fit.stem <- fit.stem.bd;
			} else {
				fit.stem <- fit.stem.yule;
			}
		}
	}
	if ((shiftCut == "node" || shiftCut == "both") && (node > root.node))
	{
		obj.node <- medusaSplitNode(node=node, z=z, desc=desc$desc.node);
		z.node <- obj.node$z;
		
		fit.node.bd <- NULL;
		fit.node.yule <- NULL;
		
		if (model == "yule" | model == "mixed")
		{
			fit.node.yule <- medusaMLFitPartition(partition.id=2, z.node, sp=c(initialR, initialE), model="yule");
			fit.node.yule$model <- "yule";
		}
		if (model == "bd" | model == "mixed")
		{
			fit.node.bd <- medusaMLFitPartition(partition.id=2, z.node, sp=c(initialR, initialE), model="bd");
			fit.node.bd$model <- "bd";
		}
## Figure out which model fits best
		if (is.null(fit.node.bd))
		{
			fit.node <- fit.node.yule;
		} else if (is.null(fit.node.yule)) {
			fit.node <- fit.node.bd;
		} else {
## Considering both places for a shift
			fit.node.bd.val <- calculateModelFit(fit=fit.node.bd, z=z);
			fit.node.yule.val <- calculateModelFit(fit=fit.node.yule, z=z);
			if (criterion == "aic") {element <- 1;} else {element <- 2;}
			if (fit.node.bd.val[[element]] < fit.node.yule.val[[element]])
			{
				fit.node <- fit.node.bd;
			} else {
				fit.node <- fit.node.yule;
			}
		}	
	}
## Now, figure out which shift position is optimal	
	if (is.null(fit.node))
	{
		fit <- fit.stem;
		cut.at <- "stem";
	} else if (is.null(fit.stem)) {
		fit <- fit.node;
		fit <- fit.node;
		cut.at <- "node";
	} else {
## Considering both places for a shift
		stem.lik <- (fit.stem$lnLik + fit.stem$lnLik);
		stem.par <- rbind(fit.stem$par, fit.stem$par)
		stem.val <- list(lnLik=stem.lik, par=stem.par);
		stem.fit <- calculateModelFit(fit=stem.val, z=z);
		
		node.lik <- (fit.node$lnLik + fit.node$lnLik);
		node.par <- rbind(fit.node$par, fit.node$par)
		node.val <- list(lnLik=node.lik, par=node.par);
		node.fit <- calculateModelFit(fit=node.val, z=z);
		
		if (criterion == "aic") {element <- 1;} else {element <- 2;}
		
		if (stem.fit[[element]] < node.fit[[element]])
		{
			fit <- fit.stem; # is this correct
			fit <- fit.stem;
			cut.at <- "stem";
		} else {
			fit <- fit.node;
			fit <- fit.node;
			cut.at <- "node";
		}
	}
	
	fit$split.at <- node;
	model.fit <- calculateModelFit(fit=fit, z=z);
	fit$aic <- model.fit[1];
	fit$aicc <- model.fit[2];
	fit$num.par <- model.fit[3];
	fit$cut.at <- cut.at;
	
	return(list(par=fit$par, lnLik=fit$lnLik, split.at=fit$split.at, aic=round(model.fit[1], digits=7),
		aicc=round(model.fit[2], digits=7), num.par=model.fit[3], cut.at=fit$cut.at, model=fit$model,
		richness=sum(z[z[,"partition"] == 1,"n.t"], na.rm=TRUE)));
}

## Only used for base model
medusaMLBaseModel <- function (z, extractClade, nodes, model, shiftCut, desc, criterion, initialR, initialE)
{
	rootnode <- min(z[,"anc"]);
	fit <- NULL;
	model.id <- NULL;
	partition.id <- 1;
	
# exclude nodes of interest from estimate of base rate
	if (extractClade == TRUE) {
		for (i in 1:length(nodes)) {
			z <- medusaSplit(node=nodes[i], z=z, desc=desc, shiftCut=shiftCut)$z;
		}
	}
		
	if (model == "bd") {
		
		fit <- medusaMLFitPartition(partition.id=1, z=z, sp=c(initialR, initialE), model="bd");
		model.fit <- calculateModelFit(fit=fitted, z=z);
		model.id <- "bd";
	} else if (model == "yule") {
		
		fit <- medusaMLFitPartition(partition.id=1, z=z, sp=c(initialR, initialE), model="yule");
		model.fit <- calculateModelFit(fit=fitted, z=z);
		model.id <- "yule";
	} else {
		
		fitted.yule <- medusaMLFitPartition(partition.id=1, z=z, sp=c(initialR, initialE), model="yule");
		fitted.bd <- medusaMLFitPartition(partition.id=1, z=z, sp=c(initialR, initialE), model="bd");
		
## Dealing with a 'mixed' model here; need to consider number of parameters.
		bd.model.fit <- calculateModelFit(fit=fitted.bd, z=z);
		yule.model.fit <- calculateModelFit(fit=fitted.yule, z=z);
		
		if (criterion == "aic") {element <- 1;} else {element <- 2;}
		
		if (bd.model.fit[[element]] < yule.model.fit[[element]])
		{
			fit <- fitted.bd;
			model.fit <- bd.model.fit;
			model.id <- "bd";
		} else {
			fit <- fitted.yule;
			model.fit <- yule.model.fit;
			model.id <- "yule";
		}
	}
	return(list(par=matrix(fit$par, nrow=1, dimnames=list(NULL,c("r", "epsilon"))), lnLik=fit$lnLik,
		split.at=rootnode, aic=round(model.fit[1], digits=7), aicc=round(model.fit[2], digits=7), num.par=model.fit[3],
		cut.at="node", model=model.id, richness=sum(z[z[,"partition"] == 1,"n.t"], na.rm=TRUE)));
}