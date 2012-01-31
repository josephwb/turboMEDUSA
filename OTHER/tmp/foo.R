require(geiger)
load("/Users/josephwb/Projects/R_working/turboMEDUSA/DATA/tmp.rda")
phy <- tmp$phy;
plotNN(phy, margin=F)

richness=NULL; modelLimit=20; stop="modelLimit"; model="mixed"; fixedEpsilon=NULL;
	criterion="aicc"; shiftCut="both"; initialR=0.05; initialE=0.5; plotFig=FALSE; nexus=FALSE;
	verbose=TRUE; mc=FALSE; numCores=NULL;


runTurboMEDUSA(tmp$phy, modelLimit=10, model.limit=10, shiftCut="stem", model="yule") -> res.yule.stem
runTurboMEDUSA(tmp$phy, modelLimit=10, model.limit=10, shiftCut="node", model="yule") -> res.yule.node
runTurboMEDUSA(tmp$phy, modelLimit=10, model.limit=10, shiftCut="both", model="yule") -> res.yule.both
runTurboMEDUSA(tmp$phy, modelLimit=10, model.limit=10, shiftCut="stem", model="bd") -> res.bd.stem
runTurboMEDUSA(tmp$phy, modelLimit=10, model.limit=10, shiftCut="node", model="bd") -> res.bd.node
runTurboMEDUSA(tmp$phy, modelLimit=10, model.limit=10, shiftCut="both", model="bd") -> res.bd.both
runTurboMEDUSA(tmp$phy, modelLimit=10, model.limit=10, shiftCut="stem", model="mixed") -> res.mixed
runTurboMEDUSA(tmp$phy, modelLimit=10, model.limit=10, shiftCut="node", model="mixed") -> res.mixed
runTurboMEDUSA(tmp$phy, modelLimit=10, model.limit=10, shiftCut="both", model="mixed") -> res.mixed
runTurboMEDUSA(tmp$phy, tmp$richness, modelLimit=10, model.limit=10, shiftCut="stem", model="bd") -> res.bd.rich
runTurboMEDUSA(tmp$phy, tmp$richness, modelLimit=10, model.limit=10, shiftCut="stem", model="yule") -> res.yule.rich
runTurboMEDUSA(tmp$phy, tmp$richness, modelLimit=10, model.limit=10, shiftCut="node", model="bd") -> res.bd.rich
runTurboMEDUSA(tmp$phy, tmp$richness, modelLimit=10, model.limit=10, shiftCut="node", model="yule") -> res.yule.rich
runTurboMEDUSA(tmp$phy, tmp$richness, modelLimit=10, model.limit=10, shiftCut="both", model="yule") -> res.mixed.rich
runTurboMEDUSA(tmp$phy, tmp$richness, modelLimit=10, model.limit=10, shiftCut="stem", model="mixed") -> res.mixed.rich
runTurboMEDUSA(tmp$phy, tmp$richness, modelLimit=10, model.limit=10, shiftCut="both", model="mixed") -> res.mixed.rich
runTurboMEDUSA(tmp$phy, tmp$richness, modelLimit=10, model.limit=10, shiftCut="both", model="mixed") -> res.mixed.rich



# fastYuleCalculation <- function (node=node, numNodes=numNodes, desc=desc, z=z, model=model)
# {
	# res <- list();
	# res <- vf(node, numNodes, desc, z, model)
# }

# numNodes <- fNv(node, phy)

# fN <- function (node, phy) {length(node.leaves(phy,node));}
# fNv <- Vectorize(fN, "node");

# fF <- function (node, numNodes, desc, z, model)
# {
	# fit <- NULL;
	# fitted.yule <- list();
	# if(model == "yule" || model == "mixed")
	# {
		# if (numNodes <= 2) { # presently only shiftCut == "node"
			# obj <- medusaSplitNode(node, z, desc=desc$desc.node);
			# res <- medusaMLFitPartition(2, z=obj$z, model="yule");
			# fitted.yule$par <- res$par;
			# fitted.yule$lnLik <- res$lnLik;
			# fitted.yule$model <- "yule";
			# fitted.yule$cut.at <- "node";
		# } else {
			# x <- numNodes - 2;
			# sum.t <- sum(z[desc$desc.node[[node]], "t.len"])
			# r <- x/sum.t;
			# fitted.yule$par <- r;
			# fitted.yule$lnLik <- x * log(r) - r * sum.t;
			# fitted.yule$model <- "yule";
			# fitted.yule$cut.at <- "node";
		# }
	# }
	# return(fitted.yule);
# }
# vf <- Vectorize(fF, c("node", "numNodes"), SIMPLIFY=FALSE);

gah <- function()
{
	for (i in 1:length(virgin.node.yule))
	{
		diff <- virgin.node.yule[[i]]$lnLik - virgin.node.bd[[i]]$lnLik;
		print(diff);
	}
}


