summarizeTurboMEDUSA <-
function(results, modelNum=NULL, cutoff="threshold", criterion="aicc", plotTree=TRUE, time=TRUE, node.labels=TRUE, cex=0.5,
	plotSurface=FALSE, printTitle=TRUE, n.points=100, ...)
{
# Desirables:
#  1. table listing parameter values of selected model
#  2. list parameters of base model
#  3. tree printed with colour-coded edges, node labels to indicate split position(s)
#  4. plot likelihood surface
	
# Extract constituent components from results
	fit <- results$models;
	phy <- results$phy;
	z <- results$z;
	desc <- results$desc;
	modelSummary <- results$modelSummary;
	threshold <- results$threshold;
	cutAt <- as.character(modelSummary$Cut.at);
	cutAt[1] <- "NA";
	
# First, determine which model is desired
	model.id <- 0;
	if (!is.null(modelNum))
	{
		model.id <- modelNum;
	} else {   # Find best model using some criterion (threshold or user-defined)
		if (cutoff != "threshold") {threshold <- cutoff}
		else {cat("\nSelecting model based on corrected threshold (improvement in information theoretic score of ", threshold, " units).\n", sep="");}
		model.id <- 1;
		while (1)
		{
			if ((model.id + 1) > length(fit)) break;
			if ((unlist(fit[[model.id]][criterion]) - unlist(fit[[model.id+1]][criterion])) < threshold) break;
			model.id <- model.id + 1;
		}
	}
	
	break.pts <- fit[[model.id]]$split.at;
	cuts <- cutAt[1: model.id];
		
	opt.model <- cbind(N.Models=seq(1:length(fit[[model.id]]$split.at)), Shift.Node=fit[[model.id]]$split.at, fit[[model.id]]$par, LnLik.part=fit[[model.id]]$lnLik.part)
	opt.model <- as.data.frame(opt.model);
	opt.model <- cbind(opt.model[,c(1:2)], Cut.at=cuts, opt.model[,c(3:5)])
	opt.model[1,2] <- NA # root node for base model
	
	base.model <- as.data.frame(fit[[1]]$par);
	
	cat("\nEstimated parameter values for model #", model.id, ":\n\n", sep="");
	print.data.frame(opt.model, digits=5);
	
	opt.fit <- 0;
	base.fit <- 0;
	if (criterion == "aicc")
	{
		opt.fit <- modelSummary$aicc[model.id];
		base.fit <- modelSummary$aicc[1];
	} else { # aic used
		opt.fit <- modelSummary$aic[model.id];
		base.fit <- modelSummary$aic[1];
	}
	cat("\nModel fit summary for model #", model.id, ":\n\n", sep="");
	cat("\tLog-likelihood = ", as.numeric(results$models[[model.id]]["lnLik"]), "\n", sep="");
	cat("\t", criterion, " = ", opt.fit, "\n", sep="");
	
	if (model.id != 1)
	{
		if (modelSummary$N.Param[1] == 1) {model <- "Yule";} else {model <- "BD";}
		cat("\nFor comparison, estimated values for the base (single homogeneous-", model, ") model are:\n\n", sep="");
		print.data.frame(base.model, digits=5, row.names=FALSE);
		cat("\nModel fit summary for base model:\n\n", sep="");
		cat("\tLog-likelihood = ", as.numeric(results$models[[1]]["lnLik"]), "\n", sep="");
		cat("\t", criterion, " = ", base.fit, "\n", sep="");
	}
	
# Get desired tree-model conformation
# Check on cutAtStem option to accurately recover model and paint tree
	labels <- break.pts;
	if (length(break.pts) > 1)
	{
		labels[1] <- break.pts[1];
		for (i in 2:length(break.pts))
		{
			tmp <- medusa.split(node=break.pts[i], z=z, desc=desc, shiftCut=cutAt[i]);
			if (cutAt[i] == "stem")
			{
				labels[i] <- z[which(z[,"dec"] == break.pts[i]),"anc"]
			}
			z <- tmp$z;
		}
	}
	
	mm <- match(phy$edge[,2], z[,"dec"]);
# Plot tree with purdy colours and labelled nodes (to better map between tree and table)
	if (plotTree)
	{
		dev.new();
		margin <- FALSE;
		
		if (time) {margin=TRUE;}
		plot.phylo(phy, edge.color=z[mm,"partition"], no.margin=!margin, cex=cex, ...);
		if (time)
		{
			axisPhylo(cex.axis=0.75);
			mtext("Divergence Time (MYA)", at=(max(get("last_plot.phylo", envir = .PlotPhyloEnv)$xx)*0.5), side = 1, line = 2, cex=0.75);
		}
		if (node.labels) # label with stem/node breaks in mind
		{
			for (i in  1:length(break.pts))
			{
				nodelabels(i, node=labels[i], frame = "c", font = 1, cex=0.5);
			}
		}
	}
	
## Due to my current ineptitude in plotting in R, plotting surface is currently [0,1] for
#  both r and epsilon, even if values of interest are clustered in some subregion.
	if (plotSurface)
	{
		n.pieces <- length(opt.model[,1]);
		for (k in 1: n.pieces)
		{	
	## *** NEED TO MAKE CHECKS IF CURRENT MODEL IS YULE, AND THEN FIGURE HOW TO PLOT IT ***
		#	if Yule, need different plot type
			if (fit[[model.id]]$model[k] == "bd")
			{
				lik <- make.lik.medusa.part(z[z[,"partition"] == k,,drop=FALSE], model="bd");
				
## center these around the MLEs
				opt.r <- as.numeric(opt.model[k,"r"]);
				opt.epsilon <- as.numeric(opt.model[k,"epsilon"]);
				
				lik.vals <- matrix(nrow=n.points, ncol=n.points);
				r.vals <- seq(from=(opt.r/10), to=(opt.r*2), length.out=n.points);
				eps.vals <- seq(from=1e-10, to=1.0, length.out=n.points);
				
				for (i in 1:length(r.vals))
				{
					for (j in 1:(length(eps.vals))) {lik.vals[i,j] <- lik(pars=c(r.vals[i], eps.vals[j]));}
				}
				if (n.pieces > 1) {cat("Completed computing surface for piecewise model #", k, "\n", sep="");}
				
# Contour plot
				plot.new()
				
				max.lik <- as.numeric(opt.model[k,"LnLik.part"])  # MLE
				
				lines <- c(max.lik-0.5, max.lik-1, max.lik-2, max.lik-3, max.lik-4, max.lik-5, max.lik-10, max.lik-50, max.lik-100)
				contour(lik.vals, levels=lines, labels=c(0.5, 1, 2, 3, 4, 5, 10, 50, 100), axes=FALSE, xlab="r (b-d)", ylab="epsilon (d/b)", method="flattest")
				tics<-floor(c(1, n.points/4, n.points/2, n.points*3/4, n.points));
				axis(1, at=c(0, 0.25, 0.5, 0.75, 1), labels=round(r.vals[tics], 3));
				axis(2, at=c(0, 0.25, 0.5, 0.75, 1), labels=round(eps.vals[tics], 3));
				points(x=opt.r*(1/(min(r.vals)+max(r.vals))), y=opt.epsilon*(1/(min(eps.vals)+max(eps.vals))), pch=16, col="black"); # inferred ML values
				if (n.pieces > 1 && printTitle) {title(main=paste("Piecewise Birth-Death Model #", k, sep=""));}
			} else {
				lik <- make.lik.medusa.part(z[z[,"partition"] == k,,drop=FALSE], model="yule");
				
				opt.r <- as.numeric(opt.model[k,"r"]);
				lik.vals <- numeric(n.points);
				r.vals <- seq(from=0, to=0.05, length.out=n.points);
#				r.vals <- seq(from=(opt.r/10), to=(opt.r*2), length.out=n.points);
				
				for (i in 1:length(r.vals))
				{
					lik.vals[i] <- lik(pars=r.vals[i]);
				}
				if (n.pieces > 1) {cat("Completed computing surface for piecewise model #", k, "\n", sep="");}
				
				max.lik <- as.numeric(opt.model[k,"LnLik.part"]);  # MLE
				
				plot(r.vals, lik.vals, xlab="r (b-d)", ylab="Log-Likelihood", type="l");
				points(x=opt.r, y=max.lik, pch=16, col="black");
				if (n.pieces > 1 && printTitle) {title(main=paste("Piecewise Yule Model #", k, sep=""));}
			}
		}
	}
	treeParameters <- list(z=z, edge.colour=mm, break.pts=break.pts, phy=phy, labels=labels);
}

