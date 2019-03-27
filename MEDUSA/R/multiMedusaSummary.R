# summarize MEDUSA results across a distribution of trees
# requires:
# 1. results (of class multiMEDUSA)
# 2. consensus tree on which to summarize results
# 3. richness information (if consensus tree does not already reflect taxonomic sampling); now a component of results
# returns: labelled tree, table with split frequencies, etc.
# cutOff pertains to displaying shift positions i.e. ignore those below cutOff

multiMedusaSummary <- function (res, conTree, cutOff=0.05, plotModelSizes=TRUE,
    plotTree=TRUE, cex=0.5, resolveConTree=FALSE, mc=FALSE, numCores=NULL, ...) {
    richness <- res$richness;
    results <- res$results;
    medusaVersion <- res$medusaVersion;
    
    if (mc) {
        if (is.null(numCores)) {
            numCores <- parallel::detectCores();
        }
        cat("Using", numCores, "cores.\n\n");
    }
    
    if (is.null(medusaVersion)) medusaVersion <- "< 0.93.4.19"; # tag MEDUSA version
    
    richness <- formatRichness(richness); # for older data sets that may have used different colnames
    
# prune consensus tree with richness information (if necessary)
    conTree <- prepareData(phy=conTree, richness=richness, verbose=FALSE, resolveTree=resolveConTree)$phy;
    
# reorder tip.labels in conTree to correspond to those in the multiMedusa analyses
    conTree <- manageTipLabels(c(results[[1]]$phy, conTree), mc=mc, numCores=numCores)[[2]];
    
# ladderize for plotting purposes
    conTree <- ladderize(conTree);
    
    num.trees <- length(results);
    cat("Summarizing MEDUSA results across ", num.trees, " trees.\n\n", sep="");
    
# check if all tip labels are in the same order (including consensus tree); makes everything easier.
    if (length(unique(lapply(lapply(results, FUN="[[", "phy"), FUN="[[", "tip.label"))) == 1 &&
        identical(results[[1]]$phy$tip.label, conTree$tip.label)) {
        cat("All translation tables identical. Summary straightforward.\n\n");
    } else {
        stop("Not all translation tables identical. Functionality not yet implemented.\n\n");
    }
        
# number of tips/edges should be same for all trees
    n.tips <- length(conTree$tip.label);
    num.edges <- length(conTree$edge[,1]);
    root.node <- n.tips + 1;
    
    cat("Processing conTree... ");
    obj <- makeCacheMedusa(phy=conTree, richness=richness, all.nodes=seq_len((2 * n.tips) -1), shiftCut="both",
                           verbose=FALSE, mc=mc, numCores=numCores);
    con.desc <- list(stem=obj$desc.stem, node=obj$desc.node);
    con.z <- obj$z;
    cat("done.\n");
    
# store tip descedants for each edge of conTree, ordered as in con.z
    
    cat("Storing tip sets... ");
    con.edge.tip.desc <- NULL;
    if (mc) {
        con.edge.tip.desc <- parallel::mclapply(con.z[,"dec"], FUN=getTips, z=con.z, desc=con.desc$stem, n.tips=n.tips,
                                      mc.cores=numCores);
    } else {
        con.edge.tip.desc <- lapply(con.z[,"dec"], FUN=getTips, z=con.z, desc=con.desc$stem, n.tips=n.tips);
    }
    cat("done");
    
    model.sizes <- numeric(num.trees);
    for (i in 1:length(results)) { # fasdt, but inelegant
        model.sizes[i] <- length(results[[i]]$optModel$split.at);
    }
    
# not terribly useful, but perhaps interesting (maybe)
    mean.n.models <- mean(model.sizes);
    sd.n.models   <- sd(model.sizes);
    min.n.models  <- min(model.sizes);
    max.n.models  <- max(model.sizes);
    if (plotModelSizes) {
        hist(model.sizes, main = NULL, xlab = "Number of Piecewise Models", prob = TRUE);
        lines(density(model.sizes, adjust = 2), lty = "longdash", col = "red");
    }
    
# for each edge in conTree, store associated estimated parameters from replicate MEDUSA results
    est.pars <- matrix(ncol=(2 * num.trees), nrow=num.edges); # important to preallocate size
    colnames(est.pars) <- rep(c("r", "epsilon"), num.trees);
    
    
# for each model size, there will be n-1 shifts.
    n.shifts <- sum(model.sizes) - num.trees;
    est.splits <- rep(NA, n.shifts); # possible for some not to map to consensus tree (i.e. incompatible)
    est.cuts <- rep(NA, n.shifts); # i.e. stem vs. node
    
# store magnitude of shift changes. consider r, b, d, and epsilon. should this be a data.frame instead?
    est.shift.magnitudes.r <- rep(NA, n.shifts);
    est.shift.magnitudes.b <- rep(NA, n.shifts);
    est.shift.magnitudes.d <- rep(NA, n.shifts);
    est.shift.magnitudes.eps <- rep(NA, n.shifts);
    
    indx.pos <- 1;
    
    for (i in 1:num.trees) {
        i.z <- results[[i]]$optModel$z;
# get tips descended from each edge in i.z
        
        i.edge.tip.desc <- NULL;
        
        if (mc) {
            i.edge.tip.desc <- parallel::mclapply(i.z[,"dec"], FUN=getTips, z=i.z, desc=results[[i]]$desc$stem, n.tips=n.tips,
                                      mc.cores=numCores);
        } else {
            i.edge.tip.desc <- lapply(i.z[,"dec"], FUN=getTips, z=i.z, desc=results[[i]]$desc$stem, n.tips=n.tips);
        }
        
        i.par <- results[[i]]$optModel$par;
    # -1 gets rid of root 'shift'
        i.splits <- results[[i]]$optModel$split.at[-1]; # need to map these
        i.cuts <- results[[i]]$optModel$cut.at[-1];
        
# use this to map edges (rows) between consensus tree and replicate trees. some will be NA.
        idx.conToRep <- match(con.edge.tip.desc, i.edge.tip.desc);
        idx.repToCon <- match(i.edge.tip.desc, con.edge.tip.desc);
        
        est.pars[, c(((2 * i) - 1), (2 * i))] <- i.par[i.z[idx.conToRep, "partition"],];
        
# now, process shift positions. may NOT exist!
# map shifted node to node in consensus tree using tip complements
        if (length(i.splits) > 0) {
            mapped.splits <- rep(NA, length(i.splits));
# all tips (i.e. < root.node) will match identical rows. could make this more concise...
            for (d in 1:length(i.splits)) {
                if (i.cuts[d] == "node") {
                    mapped.splits[d] <- as.integer(con.z[idx.repToCon[which(i.z[,"dec"] == i.splits[d])],"dec"]);
                } else { # stem shift
                    mapped.splits[d] <- as.integer(con.z[idx.repToCon[which(i.z[,"dec"] == i.splits[d])],"anc"]);
                }
            }
            
            shift.magnitudes.r <- rep(NA, length(i.splits));
            shift.magnitudes.b <- rep(NA, length(i.splits));
            shift.magnitudes.d <- rep(NA, length(i.splits));
            shift.magnitudes.eps <- rep(NA, length(i.splits));
            
            mappable.magnitude <- TRUE;
            for (k in 1:length(i.splits)) {
                if (!is.na(mapped.splits[k])) {
                    parent.class <- NULL;
                    descendant.class <- NULL;
                    mappable.magnitude <- TRUE;
                    if (i.cuts[k] == "stem" && mapped.splits[k] != root.node) { # grab rate one edge up
                        parent.node <- as.integer(i.z[which(i.z[,"dec"] == i.splits[k]), "anc"]);
                        if (parent.node != root.node) {
                            parent.class <- as.integer(i.z[which(i.z[,"dec"] == parent.node), "partition"]);
                            descendant.class <- as.integer(i.z[which(i.z[,"dec"] == i.splits[k]), "partition"]);
                        } else {
                            mappable.magnitude <- FALSE;
                        }
            # check above is always true. should be, even if deletion occurs
                    } else { # node cut
                        parent.class <- as.integer(i.z[which(i.z[,"dec"] == i.splits[k]), "partition"]);
                        descendant.class <- k + 1; # the +1 is because the initial 'shift' at the root is removed
                    }
                # got partition classes. store differences. r is easy, as always present. eps may not be.
                    if (mappable.magnitude) {
                        parms <- i.par[c(parent.class, descendant.class),]; # parent is first!
                        #shift.magnitudes.r[k] <- as.numeric(i.par[descendant.class, "r"] - i.par[parent.class, "r"]);
                        shift.magnitudes.r[k] <- as.numeric(parms[2, "r"] - parms[1, "r"]);
                        #if (any(!is.na(i.par[c(descendant.class, parent.class), "epsilon"]))) {
                        if (any(!is.na(parms[, "epsilon"]))) { # only enter if extinction is involved somewhere. otherwise, NAs.
                            parms[which(is.na(parms[,2])), 2] <- 0;
                            shift.magnitudes.eps[k] <- as.numeric(parms[2, "epsilon"] - parms[1, "epsilon"]);
                            parent.bd <- getBD(as.numeric(i.par[parent.class, 1]), as.numeric(i.par[parent.class, 2]));
                            descendant.db <- getBD(as.numeric(i.par[descendant.class, 1]), as.numeric(i.par[descendant.class, 2]));
                            
                            shift.magnitudes.b[k] <- as.numeric(descendant.db$b - parent.bd$b);
                            shift.magnitudes.d[k] <- as.numeric(descendant.db$d - parent.bd$d);
                        }
                    }
                }
            }
            
            # preallocate all of these vectors! done.
            for (j in 1:length(i.splits)) {
                est.shift.magnitudes.r[indx.pos] <- shift.magnitudes.r[j];
                
                est.shift.magnitudes.eps[indx.pos] <- shift.magnitudes.eps[j];
                est.shift.magnitudes.b[indx.pos] <- shift.magnitudes.b[j];
                est.shift.magnitudes.d[indx.pos] <- shift.magnitudes.d[j];
                
                est.splits[indx.pos] <- mapped.splits[j]; # maybe preallocate est.splits vector? probably faster!
                est.cuts[indx.pos] <- i.cuts[j];
                indx.pos <- indx.pos + 1;
            }
        }
    }
        
# summarize edge-specific rates across trees. these look okay, but shift magnitudes are wrong.
    rates <- matrix(ncol=7, nrow=num.edges);
    colnames(rates) <- c("r.mean", "r.median", "r.sd", "eps.mean", "eps.median", "eps.sd", "freq");
    
    # vectorize?
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
    conTree$rates <- rates[mapping,];
    
# summarize shift positions, mapped to consensus tree
# get rid of stem vs. node shifts
    idx.valid <- which(!is.na(est.splits));
    shift.pos <- as.data.frame(cbind(est.splits[idx.valid], est.cuts[idx.valid])); # get rid of shifts that cannot map to consensus tree
    shift.summary <- data.frame(cbind(shift.node=as.integer(rownames(table(shift.pos))), table(shift.pos)/num.trees));
    
    
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
    
    for (i in 1:length(unique.shifts)) {
        idx.shift <- which(est.splits == unique.shifts[i]);
        cur.shift.mag <- est.shift.magnitudes.r[idx.shift];
        
        mean.shift[i] <- mean(cur.shift.mag, na.rm=TRUE);
        median.shift[i] <- median(cur.shift.mag, na.rm=TRUE);
        min.shift[i] <- min(cur.shift.mag, na.rm=TRUE);
        max.shift[i] <- max(cur.shift.mag, na.rm=TRUE);
        sd.shift[i] <- sd(cur.shift.mag, na.rm=TRUE);
    }
    
    shift.summary <- cbind(shift.node=shift.summary[,1], sum.prop=(shift.summary[,"cut.at.node"] + shift.summary[,"cut.at.stem"]),
        mean.shift=mean.shift, median.shift=median.shift, min.shift=min.shift, max.shift=max.shift, sd.shift=sd.shift);
    shift.summary  <- shift.summary[order(shift.summary[,"sum.prop"], decreasing=TRUE),]; # reorder by frequency
    
# if desired, only keep shifts presultsent above cutOff threshold
    # update this to include shift magnitudes too
    shift.summary <- shift.summary[which(shift.summary[,"sum.prop"] >= cutOff),,drop=FALSE];
    rownames(shift.summary) <- NULL;
    
    
    
    ## TEMP STUFF!!!! ##
    
    # this is ugly, but for the moment just want something that works
    
# make a dataframe containing all shifts; enables plotting densities
    # maybe some other format is better than a data.frame (as nodes will have different number of shifts)
        # list?
    shifts <- as.data.frame(cbind(node=est.splits, shift.r=est.shift.magnitudes.r, shift.eps=est.shift.magnitudes.eps,
        shift.b=est.shift.magnitudes.b, shift.d=est.shift.magnitudes.d));
    shifts <- shifts[-which(is.na(shifts$node)),];
    
    idx <- shift.summary[,1];
    
    foo <- function (id) {
        tmp <- shifts[which(shifts$node == id),];
        r   <- tmp$shift.r;
        eps <- tmp$shift.eps;
        b   <- tmp$shift.b;
        d   <- tmp$shift.d;
        node <- id;
        return(list(node=node, r=r, eps=eps, b=b, d=d));
    }
    
    shifts <- lapply(idx, foo);
    
    if (length(shift.summary) == 0) {
        cat("WARNING: no node shifts occur above cutoff of ", cutOff, ". Try setting lower cutoff.\n", sep="");
    }
    
    if (cutOff > 0) {
        cat("Mapped rate shift positions present in at least ", cutOff * 100, "% (of ", num.trees, " total) trees:\n\n", sep="");
    } else {
        cat("Mapped rate shift positions across all ", num.trees, " trees:\n\n", sep="");
    }
    
    print(shift.summary);
    
    summary <- list(model.sizes=model.sizes, num.trees=num.trees, shift.summary=shift.summary,
        summary.tree=conTree, richness=richness, shift.magnitudes=shifts, medusaVersion=medusaVersion);
    class(summary) <- "multiMedusaSummary";
    
    if (plotTree) plotMultiMedusa(summary, ...);
    
    invisible(summary);
}





## Plot shift magnitude. hmm, expects summaries for all parameters. currently only exports 1 above
## not currently exported
plotShiftMagnitude <- function (summary, nodeID, par="r") {
    # map from nodeID from shift.summary
    idx <- which(summary$shift.summary[,"shift.node"] == nodeID);
    tmp <- summary$shift.magnitudes[[idx]];
    
    rr <- c("r", "eps", "b", "d");
    rl <- c("Net Diversification", "Relative Extinction", "Birth", "Death");
    
    if (length(par) == 1) {
        plot(density(tmp[[par]], na.rm=TRUE), type = "n", main="", xlab = paste("Magnitude of ", 
            rl[which(rr == par)], " Rate Shifts (Node #", nodeID, ")", sep=""));
        polygon(density(tmp[[par]], na.rm=TRUE), col = "snow3");
        rug(tmp[[par]]);
    } else { # plot multiple densities
        dd <- NULL;
        xRange <- NULL;
        yRange <- NULL;
        for (i in 1:length(par)) {
            di <- density(tmp[[par[i]]], bw="sj", na.rm=TRUE);
            dd <- c(dd, list(di));
            xRange <- c(xRange, di$x);
            yRange <- c(yRange, di$y);
        }
        
        xlim <- range(xRange, na.rm=TRUE);
        ylim <- range(0, yRange, na.rm=T);
        
        #colours <- c(rgb(0,1,0,0.5), rgb(1,0,0,0.5), rgb(0,0,1,0.5), rgb(1,0,1,0.5));
        colours <- c(rgb(0,0,1,0.5), rgb(1,0,0,0.5), rgb(0,1,0,0.5), rgb(1,0,1,0.5));
        
        plot(dd[[1]], xlim = xlim, ylim = ylim, xlab = paste("Rate Shift Magnitudes (node #", nodeID, ")", sep=""), main = "", type= "n");
        
        for (i in 1:length(par)) {
            polygon(dd[[i]], col = colours[i]);
        }
        legend('topleft', par, fill=colours[1:length(par)], bty="n");
    }
}


## TODO: use digits for rates/shifts
writeFigtree <- function (summary, file="", digits=10) {
    if (!(inherits(summary, "multiMedusaSummary"))) {
        stop("object \"summary\" must be of class \"multiMedusaSummary\"");
    }
    phy <- createAnnotatedTree(summary);
    
    res <- .write.tree3(phy, digits=digits);
    
    if (file == "") {
        return(res);
    } else {
            cat(res, file=file, sep="\n");
    }
}




# straight-up copy of ape's .write.tree2
# added bits to make Nexus. otherwise figtree chokes on Nexus annottaions in a newick tree
.write.tree3 <- function(phy, digits=10, tree.prefix="")  {
    brl <- !is.null(phy$edge.length);
    nodelab <- !is.null(phy$node.label);
    edgelab <- !is.null(phy$edge.label);
    phy$tip.label <- checkLabel(phy$tip.label); # this should be okay
    #if (nodelab) {
    #    # this is to get newick-compliant labels. not appropriate for Nexus annotation
    #    phy$node.label <- checkLabel(phy$node.label);
    #}
    f.d <- paste("%.", digits, "g", sep="");
    cp <- function(x) {
        STRING[k] <<- x;
        k <<- k + 1;
    }
    add.internal <- function (i) {
        cp("(");
        desc <- kids[[i]];
        for (j in desc) {
            if (j > n) {
                add.internal(j);
            } else {
                    add.terminal(ind[j]);
            }
            if (j != desc[length(desc)]) {
                cp(",");
            }
        }
        cp(")");
        if (nodelab && i > n) {
            cp(phy$node.label[i - n]);
        }
        ## *** add edge annotations here ***
        if (brl) {
            cp(":");
            cp(sprintf(f.d, phy$edge.length[ind[i]])); # what is this 'ind' thing?
            cp(phy$edge.label[ind[i]]);
        }
    }
    add.terminal <- function (i) {
        cp(phy$tip.label[phy$edge[i, 2]])
        ## *** add edge annotations here ***
        if (brl) {
            cp(":"); # why separate this?
            cp(sprintf(f.d, phy$edge.length[i]));
            cp(phy$edge.label[i]);
        }
    }
    n <- length(phy$tip.label);
    parent <- phy$edge[, 1];
    children <- phy$edge[, 2];
    kids <- vector("list", n + phy$Nnode);
    for (i in 1:length(parent)) {
            kids[[parent[i]]] <- c(kids[[parent[i]]], children[i]);
    }
    ind <- match(1:max(phy$edge), phy$edge[, 2]);
    LS <- 4 * n + 5;
    if (brl) {
        LS <- LS + 4 * n;
    }
    if (nodelab) {
        LS <- LS + n;
    }
    if (edgelab) {
        LS <- LS + length(phy$edge.label);
    }
    
    STRING <- character(LS);
    k <- 1;
    cp(tree.prefix);
    cp("(");
    getRoot <- function(phy) phy$edge[, 1][!match(phy$edge[, 1], phy$edge[, 2], 0)][1];
    root <- getRoot(phy);
    desc <- kids[[root]];
    for (j in desc) {
        if (j > n) {
            add.internal(j);
        } else {
                add.terminal(ind[j]);
        }
        if (j != desc[length(desc)]) 
            cp(",")
    }
    if (is.null(phy$root.edge)) {
        cp(")");
        if (nodelab) {
            cp(phy$node.label[1]);
        }
        cp(";");
    } else {
        cp(")");
        if (nodelab) {
            cp(phy$node.label[1])
        }
        cp(":");
        cp(sprintf(f.d, phy$root.edge));
        cp(";");
    }
    tstring <- paste(STRING, collapse="");
    tstring <- paste0("#NEXUS\nBegin trees;\ntree MEDUSA_TREE = [&R] ", tstring, "\nEnd;");
    return(tstring);
}




createAnnotatedTree <- function (medusa.summary) {
    phy <- medusa.summary$summary.tree;
    phy$node.label <- createNodeLabels(medusa.summary);
    phy$edge.label <- createEdgeLabels(medusa.summary);
    return (phy);
}

# combine these
createNodeLabels <- function (medusa.summary) {
    shift.summary <- medusa.summary$shift.summary;
    nodelabs <- character(medusa.summary$summary.tree$Nnode);
    ntips <- length(medusa.summary$summary.tree$tip.label);
    cols <- colnames(shift.summary);
    pos <- NULL;
    lab <- NULL;
    for (i in 1:length(shift.summary[,1])) {
        pos <- as.numeric(shift.summary[i,1]) - ntips;
        lab <- paste0("[&node_id=", shift.summary[i,1]);
        for (j in 2:length(shift.summary[i,])) {
            lab <- paste0(lab, ",", cols[j], "=", shift.summary[i,j]);
        }
        lab <- paste0(lab, "]");
        nodelabs[pos] <- lab;
    }
    return (nodelabs);
}

createEdgeLabels <- function (medusa.summary) {
    rates <- medusa.summary$summary.tree$rates;
    edgelabs <- character(length(rates[,1]));
    ntips <- length(medusa.summary$summary.tree$tip.label);
    cols <- colnames(rates);
    lab <- NULL;
    for (i in 1:length(rates[,1])) {
        lab <- "[&";
        first <- TRUE;
        for (j in 1:length(rates[i,])) {
            if (!first) {
                lab <- paste0(lab, ",", cols[j], "=", rates[i,j]);
            } else {
                lab <- paste0(lab, cols[j], "=", rates[i,j]);
                first <- FALSE;
            }
        }
        lab <- paste0(lab, "]");
        edgelabs[i] <- lab;
    }
    return (edgelabs);
}








## Plot tree with summarized results.

# TODO:
# 1. make flexible legend for rate shifts
# 2. make item placement more general/robust

plotMultiMedusa <- function (summary, treeRearrange="down", annotateShift=TRUE, annotateRate="r.median",
    plotRichnesses=TRUE, richPlot="log", time=TRUE, tip.cex=0.3, shiftScale=1, label.offset=0.5,
    font=3, shift.leg.pos="left", power=1.5, pal=1, revColours=FALSE, ...) {
    
    #dev.new(); # make a new plotting window
    conTree <- summary$summary.tree;
    shift.summary <- summary$shift.summary;
    rates <- summary$summary.tree$rates;
    richness <- summary$richness; # use for plotting species richnesses (if desired)
    
# discretize rates into some set number
    rates <- rates[,annotateRate];
    rates[which(is.nan(rates))] <- NA; # this occurs if edge in conTree occurs in no other tree
    
    #rateSeq <- seq(min(rates, na.rm=T), max(rates, na.rm=T), length=nColours);
    rateSeq <- sort(unique(rates));
    
    nColours = length(rateSeq);
    #rates <- exp(rates);
    if (pal == 1) {
        rateColours <- diverge_hcl(n=nColours, power=power);
    } else if (pal == 2) {
        rateColours <- diverge_hsv(n=nColours, power=power);
    } else if (pal == 3) {
        rateColours <- terrain_hcl(n=nColours, power=c(1/10, 1));
    } else if (pal == 4) {
        rateColours <- rainbow_hcl(n=nColours);
    } else if (pal == 5) {
        rateColours <- heat_hcl(n=nColours, power=power);
    } else if (pal == 6) {
        rateColours <- sequential_hcl(n=nColours, power=power);
    }
    
    if (revColours) {
        rateColours <- rev(rateColours); 
    }
    
# suppressWarnings is used in case some edges have no rate estimates
    edgeColours <- suppressWarnings(rateColours[unname(sapply(rates, function(x) min(which(abs(rateSeq-x) == min(abs(rateSeq-x))), na.rm=TRUE)))]);
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
    
    if (all(richness$n.taxa == 1)) plotRichnesses <- FALSE;
# plot barplot of extant tip richnesses. need to fix the mapping here.
    if (plotRichnesses) {
        richness2Plot <- NULL;
    # reorder richnesses according to tip.labels
        richness <- richness[match(conTree$tip.label, richness$taxon),];
        
        if (richPlot == "log") {
            richness2Plot <- log(richness$n.taxa + 1);
        } else {
            richness2Plot <- richness$n.taxa;
        }
        names(richness2Plot) <- richness$taxon;
        maxVal <- max(as.numeric(richness2Plot[conTree$tip.label]));
        fontSize <- lastPP$cex;
        longestName <- max(nchar(conTree$tip.label));
        nTips <- length(conTree$tip.label);
        
    # reorder (again) according to edge ordering. results from laddering conTree.
        richness2Plot <- as.numeric(richness2Plot[conTree$edge[which(conTree$edge[,2] <= length(conTree$tip.label)),2]]);
    # seems to work for large and small trees
        startPos <- max(lastPP$xx) + (longestName * fontSize) + lastPP$label.offset + max(lastPP$xx)/20;
    # adjust segment lengths
        richMultiplier <- ((maxX - startPos) / maxVal) * 0.75;
        
        segments(rep(startPos, nTips), 1:nTips, rep(startPos, nTips) + richMultiplier * richness2Plot,
            1:nTips, lwd=(tip.cex * 2), col="blue");
        
    # get max and min to find best labels
        prettyVals <- pretty(0:maxVal);
        plotAt <- startPos + (prettyVals * richMultiplier);
        axisPlacement <- mean(plotAt);
        
        axis(1, at=plotAt, labels=prettyVals, cex.axis=0.75);
        if (richPlot == "log") {
            mtext("ln(species count + 1)", at=axisPlacement, side = 1, line = 2, cex=0.75);
        } else {
            mtext("species count", at=axisPlacement, side = 1, line = 2, cex=0.75);
        }
    }
    
    if (annotateShift && (length(shift.summary) > 0)) {
        plotcolor <- rgb(red=0, green=0, blue=255, alpha=150, maxColorValue=255);
        for (i in  1:length(shift.summary[,"shift.node"])) {
            nodelabels(node=shift.summary[,"shift.node"][i], pch=21, cex=(shift.summary[,"sum.prop"][i]) * shiftScale,
                bg=plotcolor);
        }
        legend(x=shift.leg.pos, c("1.00", "0.75", "0.50", "0.25"), pch=21, pt.bg=plotcolor,
            pt.cex=(shiftScale * c(1, 0.75, 0.5, 0.25)), inset=0.05, cex=0.5, bty="n", title="Shift Frequency");
    }
    
# the weird position information used here fucks up subsequent positioning
    colorlegend(posy=c(0.30, 0.55), posx=c(0.05, 0.075), col=rateColours, zlim=minMax,
        zval=rateSeq, dz=0.5, digit=3, cex=0.25, zlevels=NULL,
        main.cex=0.5, main=paste0("Rate (", annotateRate, ")"));
}

# this function returns the phy$tip.label indices of tips descended from each edge in z
# compare this against consensus tree
getTips <- function (node, z, desc, n.tips) {
    x <- as.integer(z[desc[[node]], "dec"]); # gives descendant node(s) of a given node
    return(x[x <= n.tips]); # only return tip indices
}

print.multiMedusaSummary <- function(x, ...) {
    cat("\n");
    cat("multiMEDUSA summary results for ", x$num.trees, " trees.\n", sep="")
    cat("\n");
    
    cat("\tmedusaVersion: ", as.character(x$medusaVersion), sep="", "\n");
    cat("\tmodel.sizes: vector of optimal model sizes across all ", x$num.trees, " trees.\n", sep="");
    cat("\tsummary.tree: tree annotated with average rates across all ", x$num.trees, " trees.\n", sep="");
    cat("\tshift.summary: summary statistics for most frequent shift positions (below).\n");
    cat("\n");
    print(x$shift.summary);
}
