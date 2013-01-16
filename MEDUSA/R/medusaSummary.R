## Need to update this given new output format (i.e. only 1 optimal model)

medusaSummary <- function(results, plotTree=TRUE, time=TRUE, node.labels=TRUE, colourOffset=0, 
	cex=0.5, plotSurface=FALSE, printTitle=TRUE, n.points=100, ...) {
	
# Extract components from results
	fit <- results$optModel;
	phy <- results$phy;
	modelSummary <- results$modelSummary;
	threshold <- results$threshold;
	cutAt <- as.character(modelSummary$Cut.At);
	cutAt[1] <- "NA";
	fixPar <- results$fixPar;
	criterion <- results$criterion;
	
	break.pts <- fit$split.at;
	modelFlavour <- fit$model;
	z <- fit$z;
	
	# Get desired tree-model conformation
# Check on cutAtStem option to accurately recover model and paint tree; made more difficult by stepBack ***
	labels <- break.pts;
	if (length(break.pts) > 1) {
		labels[1] <- break.pts[1];
		for (i in 2:length(break.pts)) {
			if (cutAt[i] == "stem") {
				labels[i] <- z[which(z[,"dec"] == break.pts[i]),"anc"];
			}
		}
	}
	
# Determine parameter confidence intervals, print in table
	profLikes <- getProfileLikelihoods(z=z, parm=fit$par, models=modelFlavour, fixPar=fixPar);
		
	opt.model <- cbind(N.Models=seq(1:length(fit$split.at)), Shift.Node=fit$split.at,
		round(fit$par, digits=7), LnLik.part=fit$lnLik.part)
	opt.model <- as.data.frame(opt.model);
	if (all(modelFlavour == "yule")) {
		opt.model <- cbind(opt.model[,c(1:2)], Cut.at=cutAt, Model=modelFlavour, r=opt.model[,3], profLikes[,c(1:2)],
			opt.model[,c(4:5)])
	} else {
		opt.model <- cbind(opt.model[,c(1:2)], Cut.at=cutAt, Model=modelFlavour, r=opt.model[,3], profLikes[,c(1:2)],
			epsilon=opt.model[,4], profLikes[,c(3:4)], LnLik.part=opt.model[,5])
	}
	
	cat("\nEstimated parameter values for optimal MEDUSA model:\n\n", sep="");
	print.data.frame(opt.model, digits=5);
	cat("\n95% confidence intervals calculated from profile likelihoods\n");
	
	cat("\nModel fit summary for optimal MEDUSA model:\n\n", sep="");
	cat("\tAIC threshold = ", threshold, "\n", sep="");
	cat("\tLog-likelihood = ", fit$lnLik, "\n", sep="");
	cat("\t", criterion, " = ", as.numeric(fit[criterion]), "\n", sep="");
	cat("\tNumber of parameters = ", fit$num.par, "\n", sep="");
	if (!is.null(results$fixPar)) {
		cat("\tModel constrained to: ", fit$model[1], " with parameter fixed at ", results$fixPar, "\n", sep="");
	}
			
	mm <- match(phy$edge[,2], z[,"dec"]);
	edge.color <- z[mm, "partition"] + colourOffset; # ugh. make this nicer
	
# Plot tree with purdy colours and labelled nodes (to better map between tree and table)
	if (plotTree) {
		dev.new();
		margin <- FALSE;
		
		if (time) margin <- TRUE;
		plot.phylo(phy, edge.color=edge.color, no.margin=!margin, cex=cex, ...);
		if (time) {
			axisPhylo(cex.axis=0.75);
			mtext("Divergence Time (MYA)", at=(max(get("last_plot.phylo", envir = .PlotPhyloEnv)$xx)*0.5),
				side = 1, line = 2, cex=0.75);
		}
		if (node.labels) { # label with stem/node breaks in mind
			for (i in  1:length(break.pts)) {
				nodelabels(i, node = labels[i], frame = "c", font = 1, cex = 0.5);
			}
		}
	}
	
## *** BRING THIS OUTSIDE TO FORM ITS OWN FUNCTION ***
	if (plotSurface) {
		n.pieces <- length(opt.model[,1]);
		cat("\n");
		for (k in 1:n.pieces) {	
			if (fit$model[k] == "bd") {
				lik <- makePartitionLikelihood(z[z[,"partition"] == k,,drop=FALSE], model="bd");
## center these around the MLEs
				opt.r <- as.numeric(opt.model[k,"r"]);
				opt.epsilon <- as.numeric(opt.model[k,"epsilon"]);
				
				lik.vals <- matrix(nrow=n.points, ncol=n.points);
				r.vals <- seq(from=(opt.r/10), to=(opt.r*2), length.out=n.points);
				eps.vals <- seq(from=1e-10, to=1.0, length.out=n.points);
				
				for (i in 1:length(r.vals)) {
					for (j in 1:(length(eps.vals))) {lik.vals[i,j] <- lik(pars=c(r.vals[i], eps.vals[j]));}
				}
				if (n.pieces > 1) {cat("Completed computing surface for piecewise model #", k, "\n", sep="");}
				
# Contour plot
				dev.new();
				
				max.lik <- as.numeric(opt.model[k,"LnLik.part"])  # MLE
				
				lines <- c(max.lik-0.5, max.lik-1, max.lik-2, max.lik-3, max.lik-4, max.lik-5,
					max.lik-10, max.lik-50, max.lik-100)
				contour(lik.vals, levels=lines, labels=c(0.5, 1, 2, 3, 4, 5, 10, 50, 100), axes=FALSE,
					xlab="r (b-d)", ylab="epsilon (d/b)", method="flattest")
				tics<-floor(c(1, n.points/4, n.points/2, n.points*3/4, n.points));
				axis(1, at=c(0, 0.25, 0.5, 0.75, 1), labels=round(r.vals[tics], 3));
				axis(2, at=c(0, 0.25, 0.5, 0.75, 1), labels=round(eps.vals[tics], 3));
				points(x=opt.r*(1/(min(r.vals)+max(r.vals))), y=opt.epsilon*(1/(min(eps.vals)+max(eps.vals))),
					pch=16, col="black"); # inferred ML values
				if (n.pieces > 1 && printTitle) {title(main=paste("Piecewise Birth-Death Model #", k, sep=""));}
			} else if (fit$model[k] == "yule" || fit$model[k] == "fixedEpsilon" || fit$model[k] == "fixedR") {
				lik <- makePartitionLikelihood(z[z[,"partition"] == k,,drop=FALSE], model=fit$model[k], fixPar=fixPar);
				
				opt.val <- as.numeric(opt.model[k,"r"]);
				lik.vals <- numeric(n.points);
				
				if (opt.val > 0) {
					if (fit$model[k] == "fixedR") {
						par.vals <- seq(from=(profLikes[k,3]/2), to=(profLikes[k,4] * 2), length.out=n.points);
					} else {
						par.vals <- seq(from=(profLikes[k,1]/2), to=(profLikes[k,2] * 2), length.out=n.points);	
					}
				} else {
					par.vals <- seq(from=0, to=1e-5, length.out=n.points);
				}
				
				for (i in 1:length(par.vals)) {
					lik.vals[i] <- lik(pars=par.vals[i]);
				}
				if (n.pieces > 1) {cat("Completed computing surface for piecewise model #", k, "\n", sep="");}
				
				dev.new();
				
				max.lik <- as.numeric(opt.model[k,"LnLik.part"]);  # MLE
				
				xLabel <- "r (b-d)";
				if (fit$model[k] == "fixedR") {xLabel <- "epsilon (d/b)";}
				plot(par.vals, lik.vals, xlab=xLabel, ylab="Log-Likelihood", type="l");
				points(x=opt.val, y=max.lik, pch=16, col="black");
				if (n.pieces > 1 && printTitle) {title(main=paste("Piecewise ", fit$model[k], " Model #", k, sep=""));}
			}
		}
	}
	treeParameters <- list(z=z, edge.colour=edge.color, break.pts=break.pts, phy=phy, labels=labels,
		par=fit$par, aicc=as.numeric(fit["aicc"]));
	class(treeParameters) <- "medusa.summary";
	invisible(treeParameters);
}


## Prints out a table of likelihoods, parameters, and aic scores
## because a threshold is used, aic weights are meaningless.
optModelSummary <- function (optModel) {
	modelSize <- length(optModel$split.at);
	
	summ <- data.frame(cbind(seq(1:modelSize), optModel$split.at, optModel$cut.at, optModel$model,
		signif(optModel$lnLik.part, digits=7)), signif(optModel$par, digits=6));
	colnames(summ) <- c("Model.ID", "Shift.Node", "Cut.At", "Model", "Ln.Lik.part", "r", "epsilon");
	
	return(summ);
}


# calculateModelFitSummary <- function (models, threshold, ...) {
	# tmp <- matrix(nrow=(length(models)), ncol=6);
	# colnames(tmp) <- c("N.Models", "Shift.Node", "N.Param", "Ln.Lik", "aic", "aicc");
	
	# w.aic <- numeric(length(models));
	# w.aicc <- numeric(length(models));
	# cut.at <- character(length(models));
	# model <- character(length(models));
	
	# for (i in 1:length(tmp[,1])) {
		# tmp[i,] <- c(length(models[[i]]$split.at), tail(models[[i]]$split.at,1), models[[i]]$num.par, models[[i]]$lnLik,
			# models[[i]]$aic, models[[i]]$aicc);
		# cut.at[i] <- tail(models[[i]]$cut.at,1);
		# model[i] <- tail(models[[i]]$model,1);
	# }
	# cut.at[1] <- "NA";
	
	# all.res <- as.data.frame(tmp);
	
	# if (threshold == 0) {
		# w.aic <- round(calculateModelWeights(all.res$aic), digits=5);
		# w.aicc <- round(calculateModelWeights(all.res$aicc), digits=5);
		# all.res <- cbind(all.res[,c(1:2)], Cut.at=cut.at, Model=model, all.res[,c(3:5)], w.aic=w.aic$w,
			# aicc=all.res$aicc, w.aicc=w.aicc$w);
	# } else {
		# all.res <- cbind(all.res[,c(1:2)], Cut.at=cut.at, Model=model, all.res[,c(3:6)]);
	# }
	
	# all.res[1,2] <- NA # root node for base model

	# return(all.res);
# }


## 'fit' is a single vector of AIC scores across all models
## These are meaningless when using a threshold criterion
calculateModelWeights <- function (fit) {
	best <- min(fit);
	delta <- fit-best;
	sumDelta <- sum(exp(-0.5 * delta));
	w <- (exp(-0.5 * delta)/sumDelta);
	
	results <- data.frame(fit=fit, delta=delta, w=w);
	
	return(results);
}


## Create a plot of model-fit vs. model-size. deprecated. easily generated.
plotModelFit <- function (all.res) {
	ylim <- c(min(all.res[,"aic"],all.res[,"aicc"]), max(all.res[,"aic"],all.res[,"aicc"]));
	plot(all.res[,"N.Models"],all.res[,"aicc"], xlab="Number of Piecewise Models", ylab="Model Fit",
		ylim=ylim, type="l", col="blue");
	points(all.res[,"N.Models"],all.res[,"aicc"], col="blue", pch=21, bg="white");
	points(all.res[,"N.Models"],all.res[,"aic"], col="black", type="l");
	points(all.res[,"N.Models"],all.res[,"aic"], col="black", pch=21, bg="white");
	
	legend("topleft", c("aicc","aic"), pch=21, pt.bg="white", lty=1, col=c("blue", "black"),
		inset=0.05, cex=0.75, bty="n"); # 'bottomright' also works
}


## passed-in function fun will have form: lik <- makePartitionLikelihood(partition=new.part, model=model);
## passed-in parameters parm are MLEs stored in a matrix
getProfileLikelihoods <- function (z, parm, models, fixPar, crit=1.92) {
	res <- matrix(nrow=length(parm[,1]), ncol=4);
	colnames(res) <- c("r.low", "r.high", "eps.low", "eps.high")
	inc <- 0.05;
	
	for (i in 1:length(parm[,1])) {
#		cat("Model",i,"\n")
		model <- models[i]
		sp <- as.numeric(parm[i,]);
		new.part <- z[z[,"partition"] == i,,drop=FALSE];
		lik <- makePartitionLikelihood(partition=new.part, model=model, fixPar=fixPar);
		
		if (model == "yule") {
			par <- sp[1];
			maxLik <- lik(par); if (maxLik == -Inf) maxLik <- 0; # correct for -Inf at boundary lambda == 0
			
			threshold <- function (x) lik(x) - maxLik + crit; # find roots on either side of maxLik
			
	## need intelligent bounds
			if (par != 0) {
				low.bound <- par - par/2;
				up.bound <- par + par/2;
			} else {
				low.bound <- par;
				up.bound <- par + inc/2;
			}
			
			if (low.bound != 0) {
				while (threshold(low.bound) > 0) {
					low.bound <- low.bound - inc;
				}
			}
			while (threshold(up.bound) > 0) {
				up.bound <- up.bound + inc;
			}
			
			if (low.bound <= 0) {
				res[i,1] <- 0;
			} else {
				res[i,1] <- uniroot(threshold, lower=low.bound, upper=par)$root;
			}
			if (par == 0) par <- 1e-10; # avoid -Inf at boundary
			res[i,2] <- uniroot(threshold, lower=par, upper=up.bound)$root;
			
		} else if (model == "bd") {
			par1 <- sp[1]; par2 <- sp[2];
			maxLik <- lik(sp);
			
	## first, r
			thresholdR <- function (x) lik(c(x, par2)) - maxLik + crit;
			
			low.bound <- par1 - par1/2;
			up.bound <- par1 + par1/2;
			
			while (thresholdR(low.bound) > 0) {
				low.bound <- low.bound - inc;
			}
			while (thresholdR(up.bound) > 0) {
				up.bound <- up.bound + inc;
			}
			if (low.bound <= 0) low.bound <- 0;
			
			res[i,1] <- uniroot(thresholdR, lower=low.bound, upper=par1)$root;
			res[i,2] <- uniroot(thresholdR, lower=par1, upper=up.bound)$root;
			
	## now, epsilon
			thresholdE <- function (x) lik(c(par1, x)) - maxLik + crit;
			
			low.bound <- par2 - par2/2;
			up.bound <- par2 + par2/2;
			
			while (thresholdE(low.bound) > 0) {
				low.bound <- low.bound - inc;
			}
			while (thresholdE(up.bound) > 0) {
				up.bound <- up.bound + inc;
			}
			
			if (low.bound < 0) low.bound <- 0;
			if (up.bound > 1) up.bound <- 1;
			
			if (low.bound == 0) {
				res[i,3] <- 0;
			} else {
				res[i,3] <- uniroot(thresholdE, lower=0, upper=par2)$root;
			}
			if (up.bound == 1) {
				res[i,4] <- 1;
			} else {
				res[i,4] <- uniroot(thresholdE, lower=par2, upper=up.bound)$root;
			}
			
			if (res[i,3] == par2) res[i,3] <- 0; # precision problem?!? check optimization ***
			
		} else if (model == "fixedR" || model == "fixedEpsilon") {
			par <- sp[1]; # for model == "fixedEpsilon"
			if (model == "fixedR") {
				par <- sp[2];
			}
			
			maxLik <- lik(par);
			threshold <- function (x) lik(x) - maxLik + crit; # find roots on either side of maxLik
			
			low.bound <- par - par/2;
			up.bound <- par + par/2;
			
			while (threshold(low.bound) > 0) {
				low.bound <- low.bound - inc;
			}
			while (threshold(up.bound) > 0) {
				up.bound <- up.bound + inc;
			}
			
			if (low.bound < 0) low.bound <- 0;
			if (up.bound > 1 && model == "fixedR") up.bound <- 1;
			
			if (model == "fixedEpsilon") {
				res[i,1] <- uniroot(threshold, lower=low.bound, upper=par)$root;
				res[i,2] <- uniroot(threshold, lower=par, upper=up.bound)$root;
			} else if (model == "fixedR") {
				if (par < 1e-10 || low.bound == 0) {
					res[i,3] <- 0;
				} else {
					res[i,3] <- uniroot(threshold, lower=low.bound, upper=par)$root;
				}
				if (up.bound == 1) {
					res[i,4] <- 1;
				} else {
					res[i,4] <- uniroot(threshold, lower=par, upper=up.bound)$root;
				}
			}
			
		} else if (model == "fixedB" || model == "fixedD") {
			par <- NULL;
			if (model == "fixedB") {
				par <- BD(sp)$d;
			} else {
				par <- BD(sp)$b;
			}
			maxLik <- lik(par);
			threshold <- function (x) lik(x) - maxLik + crit; # find roots on either side of maxLik
			
			low.bound <- par - par/2;
			up.bound <- par + par/2;
			
			while (threshold(low.bound) > 0) {
				low.bound <- low.bound - inc;
			}
			while (threshold(up.bound) > 0) {
				up.bound <- up.bound + inc;
			}
			
			if (low.bound <= 0) low.bound <- 0;
			if (model == "fixedD") low.bound <- max(low.bound, fixPar);
			if (model == "fixedB") up.bound <- min(up.bound, fixPar);
			
			down <- suppressWarnings(uniroot(threshold, lower=low.bound, upper=par)$root);
			up <- suppressWarnings(uniroot(threshold, lower=par, upper=up.bound)$root);
			
			if (model == "fixedB") {
				res[i,1] <- fixPar - up;
				res[i,2] <- fixPar - down;
				res[i,3] <- down / fixPar;
				res[i,4] <- up / fixPar;
			} else { # fixedD
				res[i,1] <- down - fixPar;
				res[i,2] <- up - fixPar;
				res[i,3] <- fixPar / up;
				res[i,4] <- fixPar / down;
			}
		}
	}
	res <- as.data.frame(round(res, digits=7));
	return(res);
}


## Get b and d values from r (b-d) and epsilson (d/b)
## Used in previous version of program; now in terms of r and epsilon
## Possibly of use to users wishing to translate results
## More general version of the function getBD below
BD <- function (par1, par2=NULL) {
	if (is.null(par2)) {
		r <- par1[1];
		epsilon <- par1[2];
	} else {
		r <- par1;
		epsilon <- par2;
	}
	
	if (is.na(epsilon)) {epsilon <- 0;}
	
	b <- r/(1-epsilon);
	d <- b-r;   # Alternatively: d <- eps*r/(1-eps)
	return(list(b=b, d=d));
}


## Get b and d values from r (b-d) and epsilson (d/b)
## Used in previous version of program; now in terms of r and epsilon
## Possibly of use to users wishing to translate results
getBD <- function (r, epsilon) {
	if (is.na(epsilon)) {epsilon <- 0;}
	
	b <- r/(1-epsilon);
	d <- b-r;   # Alternatively: d <- eps*r/(1-eps)
	return(list(b=b, d=d));
}

## Print out tree with ape-style node-numbering
## Possibly of interest for users to identify node number(s) of interest
plotNN <- function (phy, time=TRUE, margin=TRUE, label.offset=0.5, cex=0.5, ...)  {
	phy$node.label <- (length(phy$tip.label) + 1):max(phy$edge);
	plot.phylo(phy, show.node.label=TRUE, no.margin=!margin, label.offset=label.offset, cex=cex, ...);
#	if (time && !margin) {cat("Cannot plot time axis without a margin.\n");}
	if (time && margin) {axisPhylo(cex.axis=0.75)};
}

print.medusa <- function(x, ...) {
	cat("\n");
	print(x$modelSummary);
	cat("\n");
}

print.multiMedusa <- function(x, ...) {
	cat("\n");
	cat("MEDUSA results for ", length(x$results), " trees.\n\n", sep="");
}

print.multiMedusaSummary <- function(x, ...) {
	cat("\n");
	cat("multiMEDUSA summary results for ", x$num.trees, " trees.\n", sep="")
	cat("\n");
	
	cat("\tmedusaVersion: ", x$medusaVersion, sep="", "\n");
	cat("\tmodel.sizes: vector of optimal model sizes across all ", x$num.trees, " trees.\n", sep="");
	cat("\tsummary.tree: tree annotated with average rates across all ", x$num.trees, " trees.\n", sep="");
	cat("\tshift.summary: summary statistics for most frequent shift positions (below).\n");
	cat("\n");
	print(x$shift.summary);
}