# summarize MEDUSA results across a distribution of trees
# requires:
# 1. results (of class multiMEDUSA)
# 2. consensus tree on which to summarize results
# 3. richness information (if consensus tree does not already reflect taxonomic sampling); now a component of results
# returns: labelled tree, table with split frequencies, etc.
# cutOff pertains to displaying shift positions i.e. ignore those below cutOff

# TODO: ladderize tree prior to summary - DONE


multiMedusaSummary <- function (res, conTree, cutOff=0.05, plotModelSizes=FALSE,
	cex=0.5, ...) {
	richness <- res$richness;
	results <- res$results;
	
# prune consensus tree with richness information (if possible)
	conTree <- MEDUSA:::prepareData(phy=conTree, richness=richness, verbose=FALSE)$phy;
	
# number of tips/edges should be same for all trees
	n.tips <- length(conTree$tip.label);
	num.edges <- length(conTree$edge[,1]);
	root.node <- n.tips + 1;
	
	obj <-MEDUSA::: makeCacheMedusa(phy=conTree, richness=richness, all.nodes=seq_len((2 * n.tips) -1), shiftCut="both", verbose=FALSE, mc=F);
	con.desc <- list(stem=obj$desc.stem, node=obj$desc.node);
	con.z <- obj$z;
	
# store tip descedants for each edge of conTree, ordered as in con.z
	con.edge.tip.desc <- lapply(con.z[,"dec"], FUN=getTips, z=con.z, desc=con.desc$stem, n.tips=n.tips);
	
	num.trees <- length(results);
	
# check if all tip labels are in the same order (including consensus tree); makes everything easier
	pruned.trees <- lapply(results, FUN="[[", "phy"); names(pruned.trees) <- NULL;
	tipLabels <- lapply(pruned.trees, FUN="[[", "tip.label");
	if (length(unique(tipLabels)) == 1 && identical(tipLabels[[1]], conTree$tip.label)) {
		cat("All translation tables identical. Summary straighforward.\n\n");
	} else {
		stop("Not all translation tables identical. Functionality not yet implemented.\n\n");
	}

# this will change if I alter MEDUSA output format to keep only final model
# will anyone ever really want an intermediate model? may be a lot of data to store for very large trees
	model.sizes <- numeric(num.trees);
	for (i in 1:length(results)) {
		model.sizes[i] <- length(results[[i]]$models);
	}
	
# not terribly useful, but perhaps interesultsting (maybe)
	mean.n.models <- mean(model.sizes);
	sd.n.models   <- sd(model.sizes);
	min.n.models  <- min(model.sizes);
	max.n.models  <- max(model.sizes);
	if (plotModelSizes) hist(model.sizes, main=NULL, xlab="Number of Piecewise Models");
	
# keep only optimal (last) model for each tree; other bits (e.g. phy, desc, etc.) remain
	cleaned.results <- results;
	for (i in 1:num.trees) {
		cleaned.results[[i]]$models <- cleaned.results[[i]]$models[[model.sizes[i]]];
	}

# for each edge in conTree, store associated estimated parameters from replicate MEDUSA results
	est.pars <- matrix(ncol=(2 * num.trees), nrow=num.edges); colnames(est.pars) <- rep(c("r", "epsilon"), num.trees)
	est.splits <- NULL; # possible for some not to map to consensus tree (i.e. incompatible)
	est.cuts <- NULL; # i.e. stem vs. node
	
	for (i in 1:num.trees) {
		i.z <- cleaned.results[[i]]$models$z;
		i.desc <- cleaned.results[[i]]$desc;
# get tips descended from each edge in i.z
		i.edge.tip.desc <- lapply(i.z[,"dec"], FUN=getTips, z=i.z, desc=i.desc$stem, n.tips=n.tips);
		i.par <- cleaned.results[[i]]$models$par;
		i.splits <- cleaned.results[[i]]$models$split.at;
		i.cuts <- cleaned.results[[i]]$models$cut.at;
		
# use this to map edges between replicate trees and consensus tree. some may be NA.
		idx <- match(con.edge.tip.desc, i.edge.tip.desc);
		
		temp.par <- matrix(ncol=2, nrow=num.edges); colnames(temp.par) <- c("r", "epsilon");
		
		for (j in 1:length(idx)) {
			part <- i.z[idx[j], "partition"];
			
			#cat("Model membership is: ", part, "\n", sep="");
			temp.par[j,] <- i.par[part,];
		}
		
		est.pars[, c(((2 * i) - 1), (2 * i))] <- temp.par;
		
# now, process shift positions
# map shifted node to node in consensus tree using tip complements
		mapped.splits <- rep(NA, length(i.splits));
		idx.mapped <- which(i.splits <= root.node);
		mapped.splits[idx.mapped] <- i.splits[idx.mapped];
		unmapped <- which(is.na(mapped.splits));
		
		for(j in 1:length(unmapped)) {
			x <- which(i.z[,"dec"] == i.splits[unmapped[j]]); # row/edge index
			if (!length(which(idx == x)) < 1) {
				mapped.splits[unmapped[j]] <- as.integer(con.z[which(idx == x), "dec"]);
			}
		}
		est.splits <- c(est.splits, mapped.splits);
		est.cuts <- c(est.cuts, i.cuts);
	}
	
# summarize edge-specific rates across trees
	rates <- matrix(ncol=7, nrow=num.edges); colnames(rates) <- c("r.mean", "r.median", "r.sd", "eps.mean", "eps.median", "eps.sd", "freq");
	
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
	mapping <- match(conTree$edge[, 2], con.z[, 2]);
	#all(conTree$edge[,] == as.integer(con.z[mapping, c(1:2)])); # check
	conTree$rates <- rates[mapping,];
	
# summarize shift positions, mapped to consensus tree
	shift.pos <- as.data.frame(cbind(est.splits[!is.na(est.splits)], est.cuts[!is.na(est.splits)])); # get rid of shifts that cannot map to consensus tree
	shift.summary <- data.frame(cbind(shift.node=as.integer(rownames(table(shift.pos))), table(shift.pos)/num.trees));
	colnames(shift.summary)[2:3] <- c("cut.at.node", "cut.at.stem");
	
	#sum.prop <- shift.summary[,"cut.at.node"] + shift.summary[,"cut.at.stem"];
	shift.summary <- cbind(shift.summary, sum.prop=(shift.summary[,"cut.at.node"] + shift.summary[,"cut.at.stem"]));
	shift.summary  <- shift.summary[order(shift.summary[,"sum.prop"], decreasing=TRUE),]; # reorder by frequency
	
# if desired, only keep shifts presultsent above cutOff thresultshold
	shift.summary <- shift.summary[which(shift.summary["sum.prop"] >= cutOff),];
	rownames(shift.summary) <- NULL;
	
	if (cutOff > 0) {
		cat("Mapped rate shift positions present in at least ", cutOff * 100, "% (of ", num.trees, " total) trees:\n\n", sep="");
	} else {
		cat("Mapped rate shift positions across all ", num.trees, " trees:\n\n", sep="");
	}
	
	print(shift.summary);
	
	summary <- list(model.sizes=model.sizes, num.trees=num.trees, shift.summary=shift.summary, summary.tree=conTree);
	class(summary) <- "multiMedusaSummary";
	
	invisible(summary);
}

## Plot tree with summarized summaryults.
plot.multiMedusaSummary <- function (summary, treeRearrange="up", annotateShift=TRUE, annotateRate="r.mean", plotRate=FALSE,
	time=TRUE, cex=0.5, shiftScale=1, label.offset=0.5, font=3, shift.leg.pos="bottomleft", power=1, ...) {
	conTree <- summary$summary.tree;
	shift.summary <- summary$shift.summary;
	rates <- summary$summary.tree$rates;
	
# rearrange tree for better viewing and placement of legend(s)
	if (!is.null(treeRearrange)) {
		if (treeRearrange == "up") {
			conTree <- ladderize(conTree, FALSE);
		} else {
			conTree <- ladderize(conTree);
		}
	}
	
# discretize rates into some set number
	rates <- rates[,annotateRate];
	
	rateColours <- diverge_hcl(20, power=power); # might change number of colours here
	rateSeq <- seq(min(rates), max(rates), length=20);
	edgeColours <- rateColours[unname(sapply(rates, function(x) min(which(abs(rateSeq-x) == min(abs(rateSeq-x))))))]
	
# shift positions (with label size proportional to frequency)
	#dev.new();
	margin <- FALSE; if (time) margin <- TRUE;
	
	plot.phylo(conTree, edge.color=edgeColours, no.margin=!margin, cex=cex, label.offset=label.offset, font=font, ...);
	if (time) {
		axisPhylo(cex.axis=0.75);
		mtext("Divergence Time (MYA)", at=(max(get("last_plot.phylo", envir = .PlotPhyloEnv)$xx)*0.5),
			side = 1, line = 2, cex=0.75);
	}
	
	if (annotateShift) {
		plotcolor <- rgb(red=255, green=0, blue=0, alpha=150, max=255);
		for (i in  1:length(shift.summary[,"shift.node"])) {
			nodelabels(node = shift.summary[,"shift.node"][i], pch=21, cex=(shift.summary[,"sum.prop"][i]) * shiftScale, bg=plotcolor);
		}
		
		legend(x=shift.leg.pos, c("1.0", "0.5", "0.1"), pch=21, pt.bg=plotcolor, pt.cex=(shiftScale * c(1, 0.5, 0.1)),
			inset=0.05, cex=(cex * shiftScale), bty="n", title="Shifts");
	}
}

# this function returns the phy$tip.label indices of tips decended from each edge in z
# compare this against consensus tree
getTips <- function (node, z, desc, n.tips) {
	x <- as.integer(z[desc[[node]], "dec"]) # gives descendant node(s) of a given node
	return(x[x <= n.tips]); # only return tip indices
}