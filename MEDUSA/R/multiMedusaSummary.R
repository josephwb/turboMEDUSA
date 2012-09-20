# summarize MEDUSA results across a distribution of trees
# requires:
# 1. results (of class multiMEDUSA)
# 2. consensus tree on which to summarize results
# 3. richness information (if consensus tree does not already reflect taxonomic sampling); now a component of results
# returns: labelled tree, table with split frequencies, etc.
# cutOff pertains to displaying shift positions i.e. ignore those below cutOff

# TODO:
# 1. ladderize tree prior to summary - DONE
# 2. rate legend - DONE
# 3. shift proportion legend - DONE

# get rid of all stem vs. node stuff

multiMedusaSummary <- function (res, conTree, cutOff=0.05, plotModelSizes=TRUE,
	plotTree=TRUE, cex=0.5, resolveConTree=FALSE, ...) {
	richness <- res$richness;
	results <- res$results;
	
# prune consensus tree with richness information (if necessary)
# an issue here is that tip label ordering in conTree may differ from those from the results (which probably used a translation table)
# NEED to fix this, as it may be a general problem.
	conTree <- prepareData(phy=conTree, richness=richness, verbose=FALSE, resolveTree=resolveConTree)$phy;
	
# reorder tip.labels in conTree to correspond to those in the multiMedusa analyses
	conTree <- ape:::.compressTipLabel(c(results[[1]]$phy, conTree))[[2]];
	
	num.trees <- length(results);
	cat("Summarizing MEDUSA results across ", num.trees, " trees.\n\n", sep="");
	
# check if all tip labels are in the same order (including consensus tree); makes everything easier.
# tree distribution. This should be checked earlier i.e. prior to analysis.
	tipLabels <- lapply(lapply(results, FUN="[[", "phy"), FUN="[[", "tip.label")
	
	if (length(unique(tipLabels)) == 1 && identical(tipLabels[[1]], conTree$tip.label)) {
		cat("All translation tables identical. Summary straightforward.\n\n");
	} else {
		stop("Not all translation tables identical. Functionality not yet implemented.\n\n");
	}

	rm(tipLabels);
	
# ladderize for plotting purposes
	conTree <- ladderize(conTree);
	
# number of tips/edges should be same for all trees
	n.tips <- length(conTree$tip.label);
	num.edges <- length(conTree$edge[,1]);
	root.node <- n.tips + 1;
	
	obj <- makeCacheMedusa(phy=conTree, richness=richness, all.nodes=seq_len((2 * n.tips) -1), shiftCut="both", verbose=FALSE, mc=F);
	con.desc <- list(stem=obj$desc.stem, node=obj$desc.node);
	con.z <- obj$z;
	
# store tip descedants for each edge of conTree, ordered as in con.z
	con.edge.tip.desc <- lapply(con.z[,"dec"], FUN=getTips, z=con.z, desc=con.desc$stem, n.tips=n.tips);
	
# this will change if I alter MEDUSA output format to keep only final model
# will anyone ever really want an intermediate model? may be a lot of data to store for very large trees
	model.sizes <- numeric(num.trees);
	for (i in 1:length(results)) {
		model.sizes[i] <- length(results[[i]]$optModel$split.at);
	}
	
# not terribly useful, but perhaps interesultsting (maybe)
	mean.n.models <- mean(model.sizes);
	sd.n.models   <- sd(model.sizes);
	min.n.models  <- min(model.sizes);
	max.n.models  <- max(model.sizes);
	if (plotModelSizes) hist(model.sizes, main=NULL, xlab="Number of Piecewise Models");
	
# keep only optimal (last) model for each tree; other bits (e.g. phy, desc, etc.) remain
	cleaned.results <- results; # this is deprecated, as only best model is now saved
	
# for each edge in conTree, store associated estimated parameters from replicate MEDUSA results
	est.pars <- matrix(ncol=(2 * num.trees), nrow=num.edges);
	colnames(est.pars) <- rep(c("r", "epsilon"), num.trees);
	est.splits <- NULL; # possible for some not to map to consensus tree (i.e. incompatible)
	est.cuts <- NULL; # i.e. stem vs. node
	est.shift.magnitudes <- NULL; # store magnitude of shift changes
	
	for (i in 1:num.trees) {
		i.z <- cleaned.results[[i]]$optModel$z;
		i.desc <- cleaned.results[[i]]$desc;
# get tips descended from each edge in i.z
		i.edge.tip.desc <- lapply(i.z[,"dec"], FUN=getTips, z=i.z, desc=i.desc$stem, n.tips=n.tips);
		i.par <- cleaned.results[[i]]$optModel$par;
		i.splits <- cleaned.results[[i]]$optModel$split.at[-1]; # need to map these
		i.cuts <- cleaned.results[[i]]$optModel$cut.at[-1]; # -1 gets rid of root 'shift'
		
# use this to map edges (rows) between consensus tree and replicate trees. some may be NA.
		idx.conToRep <- match(con.edge.tip.desc, i.edge.tip.desc);
		idx.repToCon <- match(i.edge.tip.desc, con.edge.tip.desc);
		
		temp.par <- matrix(ncol=2, nrow=num.edges); colnames(temp.par) <- c("r", "epsilon");
		
		for (j in 1:length(idx.conToRep)) {
			part <- i.z[idx.conToRep[j], "partition"]; # put in order of conTree edges
			
			#cat("Model membership is: ", part, "\n", sep="");
			temp.par[j,] <- i.par[part,];
		}
		
		est.pars[, c(((2 * i) - 1), (2 * i))] <- temp.par;		
		
# now, process shift positions
# map shifted node to node in consensus tree using tip complements
		mapped.splits <- rep(NA, length(i.splits));
		
# all tips (i.e. < root.node) will match identical rows
		for (d in 1:length(i.splits)) {
			# shiftNode <- i.splits[d]; # good.
			# shiftNodePosition <- which(i.z[,"dec"] == i.splits[d]); # good.
			# mappedShiftEdge <- idx.repToCon[shiftNodePosition]; # now, good.
			# mappedShiftNode <- con.z[mappedShiftEdge,]; # got it.
			
			if (i.cuts[d] == "node") {
				mapped.splits[d] <- as.integer(con.z[idx.repToCon[which(i.z[,"dec"] == i.splits[d])],"dec"]);
			} else { # stem shift
				mapped.splits[d] <- as.integer(con.z[idx.repToCon[which(i.z[,"dec"] == i.splits[d])],"anc"]);
			}
		}
		
		est.splits <- c(est.splits, mapped.splits);
		est.cuts <- c(est.cuts, i.cuts);
		
		shift.magnitudes <- rep(NA, length(i.splits));
		for (k in 1:length(i.splits)) {
			if (!is.na(mapped.splits[k])) {
				
				parent.class <- NULL;
				descendant.class <- NULL;
				
				if (i.cuts[k] == "stem") { # grab rate one node up
					parent.node <- as.integer(i.z[which(i.z[,"dec"] == i.splits[k]), "anc"]);
					parent.class <- as.integer(i.z[which(i.z[,"dec"] == parent.node), "partition"]);
					descendant.class <- as.integer(i.z[which(i.z[,"dec"] == i.splits[k]), "partition"]);
					#descendant.class <- k; # by definition. hmm, maybe not if a stepback occurred...
		# check above is always true. should be, even if deletion occurs
				} else { # node cut
					parent.class <- as.integer(i.z[which(i.z[,"dec"] == i.splits[k]), "partition"]);
					descendant.class <- k; # + 1; # the +1 is because the initial 'shift' at the root is removed
				}
			# got partition classes. store differences.
				shift.magnitudes[k] <- as.numeric(i.par[descendant.class, "r"] - i.par[parent.class, "r"]);
			#	cat("Tree #", i, ": comparing rate ", as.numeric(i.par[descendant.class, "r"]), 
			#		" and ", as.numeric(i.par[parent.class, "r"]), "\n", sep="");
			}
		}
		est.shift.magnitudes <- c(est.shift.magnitudes, shift.magnitudes);
	}
	
# summarize edge-specific rates across trees
	rates <- matrix(ncol=7, nrow=num.edges);
	colnames(rates) <- c("r.mean", "r.median", "r.sd", "eps.mean", "eps.median", "eps.sd", "freq");
	
	for (i in 1:num.edges) {
		i.r <- as.numeric(est.pars[i, seq(from=1, to=(num.trees * 2), by=2)]);
		idx.valid <- !is.na(i.r);
		i.eps <- as.numeric(est.pars[i, seq(from=2, to=(num.trees * 2), by=2)]);
		i.eps <- i.eps[idx.valid]; i.eps[which(is.na(i.eps))] <- 0;
		
		r.mean <- mean(i.r, na.rm=TRUE);
		r.median <- median(i.r, na.rm=TRUE);
		r.sd <- sd(i.r, na.rm=TRUE);
		eps.mean <- mean(i.eps);
		eps.median <- median(i.eps);
		eps.sd <- sd(i.eps);
		freq <- sum(idx.valid) / num.trees; # how often the particular edge occurs across trees
		
		rates[i,] <- c(r.mean, r.median, r.sd, eps.mean, eps.median, eps.sd, freq);
	}
	
# map con.z edges to conTree edges, annotate tree with edge-specific rates
	mapping <- match(conTree$edge[, 2], con.z[, 2]); # hmm. is this correct?!? yes, should be.
	#all(conTree$edge[,] == as.integer(con.z[mapping, c(1:2)])); # check
	conTree$rates <- rates[mapping,];
	
# summarize shift positions, mapped to consensus tree
# get rid of stem vs. node shifts
# seems to have a problem if the are zero 'node' shifts (same for the opposite?)
	idx.valid <- which(!is.na(est.splits));
	shift.pos <- as.data.frame(cbind(est.splits[idx.valid], est.cuts[idx.valid])); # get rid of shifts that cannot map to consensus tree
	shift.summary <- data.frame(cbind(shift.node=as.integer(rownames(table(shift.pos))), table(shift.pos)/num.trees));
	
# problem here when either stem or node is never observed
	if (length(shift.summary[1,]) < 3) {
		if (is.null(shift.summary$node)) {
			shift.summary <- cbind(shift.node=shift.summary[,1], node=rep(0, length(shift.summary[,1])),
				stem=shift.summary$stem);
		} else if (is.null(shift.summary$stem)) {
			shift.summary <- cbind(shift.summary[,1:2], stem=rep(0, length(shift.summary[,1])));
		} else {
			stop("\nUm, I don't know what is wrong here.\n")
		}
	}
	
	colnames(shift.summary)[2:3] <- c("cut.at.node", "cut.at.stem");
	
	unique.shifts <- shift.summary[,"shift.node"];
	num.unique.shifts <- length(unique.shifts);
	
	mean.shift <- rep(NA, num.unique.shifts);
	median.shift <- rep(NA, num.unique.shifts);
	min.shift <- rep(NA, num.unique.shifts);
	max.shift <- rep(NA, num.unique.shifts);
	sd.shift <- rep(NA, num.unique.shifts);
	
	for (i in 1:length(unique.shifts)) { # first will always be the root, which is not a shift
		idx.shift <- which(est.splits == unique.shifts[i])
		cur.shift.mag <- est.shift.magnitudes[idx.shift];
		
		mean.shift[i] <- mean(cur.shift.mag, na.rm=TRUE);
		median.shift[i] <- median(cur.shift.mag, na.rm=TRUE);
		min.shift[i] <- min(cur.shift.mag, na.rm=TRUE);
		max.shift[i] <- max(cur.shift.mag, na.rm=TRUE);
		sd.shift[i] <- sd(cur.shift.mag, na.rm=TRUE);
	}
	
	shift.summary <- cbind(shift.node=shift.summary[,1], sum.prop=(shift.summary[,"cut.at.node"] + shift.summary[,"cut.at.stem"]),
		mean.shift=mean.shift, median.shift=median.shift, min.shift=min.shift, max.shift=max.shift, sd.shift=sd.shift);
	shift.summary  <- shift.summary[order(shift.summary[,"sum.prop"], decreasing=TRUE),]; # reorder by frequency
	
# if desired, only keep shifts presultsent above cutOff thresultshold
	shift.summary <- shift.summary[which(shift.summary[,"sum.prop"] >= cutOff),];
	rownames(shift.summary) <- NULL;
	
	if (cutOff > 0) {
		cat("Mapped rate shift positions present in at least ", cutOff * 100, "% (of ", num.trees, " total) trees:\n\n", sep="");
	} else {
		cat("Mapped rate shift positions across all ", num.trees, " trees:\n\n", sep="");
	}
	
	print(shift.summary);
	
	summary <- list(model.sizes=model.sizes, num.trees=num.trees, shift.summary=shift.summary,
		summary.tree=conTree, richness=richness);
	class(summary) <- "multiMedusaSummary";
	
	if (plotTree) plotMultiMedusa(summary, ...);
	
	invisible(summary);
}


## Plot tree with summarized results.

# TODO:
# 1. make flexible legend for rate shifts
# 2. make item placement more general/robust

plotMultiMedusa <- function (summary, treeRearrange="down", annotateShift=TRUE, annotateRate="r.median",
	plotRichnesses=TRUE, richPlot="log", time=TRUE, tip.cex=0.3, shiftScale=1, label.offset=0.5,
	font=3, shift.leg.pos="left", power=1.5, ...) {
	
	dev.new(); # make a new plotting window
	conTree <- summary$summary.tree;
	shift.summary <- summary$shift.summary;
	rates <- summary$summary.tree$rates;
	richness <- summary$richness; # use for pplotting species richnesses (if desired)
	
# discretize rates into some set number
	rates <- rates[,annotateRate];
	rates[which(is.nan(rates))] <- NA; # this occurs if edge in conTree occurs in no other tree
	
	rateColours <- diverge_hcl(20, power=power); # might change number of colours here
	rateSeq <- seq(min(rates, na.rm=T), max(rates, na.rm=T), length=20);
	
# suppressWarnings is used in case some edges have no rate estimates
	edgeColours <- suppressWarnings(rateColours[unname(sapply(rates, function(x) min(which(abs(rateSeq-x) == min(abs(rateSeq-x))), na.rm=T)))]);
	edgeColours[which(is.na(edgeColours))] <- "#000000"; # set to black those without estimates. shouldn't happen with a decent tree.
	
	minMax <- c(min(rateSeq), max(rateSeq));
	
# shift positions (with label size proportional to frequency)
	margin <- FALSE; if (time) margin <- TRUE;
	
	plot.phylo(conTree, edge.color=edgeColours, no.margin=!margin, cex=tip.cex, label.offset=label.offset, font=font, ...);
	
# store parameters from tree plotting
	lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv);
	maxX <- lastPP$x.lim[2];
	
	if (time) {
		axisPhylo(cex.axis=0.75);
		mtext("Divergence Time (MYA)", at=(max(lastPP$xx) * 0.5), side=1, line=2, cex=0.75);
	}
	
# plot barplot of extant tip richnesses. need to fix the mapping here.
	if (plotRichnesses) {
		richness2Plot <- NULL;
		
	# reorder richnesses according to tip.labels
		richness <- richness[match(conTree$tip.label, richness[,1]),];
		
		if (richPlot == "log") {
			richness2Plot <- log(richness[,2] + 1);
			names(richness2Plot) <- richness[,1];
		} else {
			richness2Plot <- richness[,2];
			names(richness2Plot) <- richness[,1];
		}
		maxVal <- max(as.numeric(richness2Plot[conTree$tip.label]));
		fontSize <- lastPP$cex;
		longestName <- max(str_length(conTree$tip.label));
		nTips <- length(conTree$tip.label);
		
	# reorder (again) according to edge ordering. results from laddering conTree.
		richness2Plot <- as.numeric(richness2Plot[conTree$edge[which(conTree$edge[,2] <= length(conTree$tip.label)),2]]);
		
	# seems to work for large and small trees
		startPos <- max(lastPP$xx) + (longestName * fontSize) + lastPP$label.offset + max(lastPP$xx)/20;
		
	# adjest segment lengths
		richMultiplier <- ((maxX - startPos) / maxVal) * 0.75;
		
		segments(rep(startPos, nTips), 1:nTips, rep(startPos, nTips) + richMultiplier * richness2Plot,
			1:nTips, lwd=(tip.cex * 2), col="blue");
		
	# get max and min to find best labels
		prettyVals <- pretty(0:maxVal);
		plotAt <- startPos + (prettyVals * richMultiplier);
		axisPlacement <- mean(plotAt);
		
		axis(1, at=plotAt, labels=prettyVals, cex.axis=0.75);
		mtext("ln(species count + 1)", at=axisPlacement, side = 1, line = 2, cex=0.75);
	}
	
	if (annotateShift) {
		plotcolor <- rgb(red=255, green=0, blue=0, alpha=150, maxColorValue=255);
		for (i in  1:length(shift.summary[,"shift.node"])) {
			nodelabels(node=shift.summary[,"shift.node"][i], pch=21, cex=(shift.summary[,"sum.prop"][i]) * shiftScale,
				bg=plotcolor);
		}
		
		legend(x=shift.leg.pos, c("1.00", "0.75", "0.50", "0.25"), pch=21, pt.bg=plotcolor,
			pt.cex=(shiftScale * c(1, 0.75, 0.5, 0.25)), inset=0.05, cex=0.5, bty="n", title="Shift Frequency");
	}
	
# the weird position information used here fucks up subsequent positioning
	colorlegend(posy=c(0.30, 0.55), posx=c(0.05, 0.075), col=rateColours, zlim=minMax,
		zval=seq(min(rateSeq), max(rateSeq), length=10), dz=0.5, digit=3, cex=0.5, zlevels=NULL,
		main.cex=0.5, main="Rate");
}

# this function returns the phy$tip.label indices of tips decended from each edge in z
# compare this against consensus tree
getTips <- function (node, z, desc, n.tips) {
	x <- as.integer(z[desc[[node]], "dec"]); # gives descendant node(s) of a given node
	return(x[x <= n.tips]); # only return tip indices
}