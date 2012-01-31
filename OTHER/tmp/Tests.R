require(geiger);
load("/Users/josephwb/Projects/R_working/turboMEDUSA/OTHER/tmp.rda");
phy <- tmp$phy;
richness <- tmp$richness;
richness=NULL; modelLimit=5; stop="modelLimit"; model="mixed"; fixEpsilon=NULL; fixR=NULL; criterion="aicc"; shiftCut="both"; userNode=NULL; initialR=0.05; initialE=0.5; plotFig=FALSE; verbose=TRUE; mc=FALSE; numCores=NULL; epsilon=NULL; r=NULL; b=NULL; d=NULL; preserveModelFlavour=T; modelNum=NULL; cutoff="threshold"; crit=1.92; stepBack=TRUE; modelNum=NULL; cutoff="threshold"; criterion="aicc"; plotTree=TRUE; time=TRUE; node.labels=TRUE; cex=0.5; plotSurface=FALSE; printTitle=TRUE; n.points=100; fixThreshold=NULL; stepBack=T;

check.multicore<-function () 
{
    tmp = rownames(installed.packages())
    if ("multicore" %in% tmp) {
        require(multicore)
        return(TRUE)
    }
    else {
        return(FALSE)
    }
}

likelihoods.tips <- NULL;
par.tips <- NULL;
for (i in 1:length(tips))
{
	cat("Node=", i, "; lnLik=", tips[[i]]$lnLik, "\n", sep="");
	likelihoods.tips[i] <- tips[[i]]$lnLik;
	par.tips[i] <- tips[[i]]$par[1];
}

plot(likelihoods.tips)
plot(par.tips)

likelihoods.res <- NULL;
for (i in 1:length(res))
{
	cat("Node=", i, "; lnLik=", res[[i]]$lnLik, "\n", sep="");
	likelihoods.res[i] <- res[[i]]$lnLik;
}

dev.new()
plot(likelihoods.res)

likelihoods.res.tips <- NULL;
for (i in 1:length(tips))
{
	cat("Node=", i, "; lnLik=", res[[i]]$lnLik, "\n", sep="");
	likelihoods.res.tips[i] <- res[[i]]$lnLik;
}

dev.new()
plot(likelihoods.res.tips)

likelihoods.virgin.stem <- NULL;
likelihoods.virgin.node <- NULL;
likelihoods.virgin.diff <- NULL;
for (i in 1:length(virgin.node))
{
	cat("Node=", i, "; lnLik=", virgin.nodes$node[[i]]$lnLik, "\n", sep="");
	cat("Node=", i, "; lnLik=", virgin.nodes$stem[[i]]$lnLik, "\n", sep="");
	likelihoods.virgin.node[i] <- virgin.nodes$node[[i]]$lnLik;
	likelihoods.virgin.stem[i] <- virgin.nodes$stem[[i]]$lnLik;
	likelihoods.virgin.diff[i] <- virgin.nodes$node[[i]]$lnLik - virgin.nodes$stem[[i]]$lnLik;
}

dev.new()
plot(likelihoods.virgin.stem)
dev.new()
plot(likelihoods.virgin.node)
dev.new()
plot(likelihoods.virgin.diff)
plot(likelihoods.virgin.stem, likelihoods.virgin.node)

medusaSplitStem <- function (node, z, desc)
{
	part <- z[,"partition"];
	base <- min(part[z[,"anc"] == node | z[,"dec"] == node]);
	tag <- max(part) + 1;
	i <- desc[[node]];
	idx <- i[part[i] == base];
	z[idx,"partition"] <- tag;
	z[which(z["dec"] == node),"partition"] <- tag;
	return(list(z=z, affected=c(unique(part[idx]), tag)));
}

z.orig
partition.id
