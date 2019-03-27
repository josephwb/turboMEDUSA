medusaFitClade <- function (phy=phy, richness=NULL, userNode=NULL, model="mixed",
    epsilon=NULL, r=NULL, b=NULL, d=NULL, criterion="aicc", shiftCut="both",
    extractClade=FALSE, initialR=0.05, initialE=0.5, verbose=TRUE, mc=FALSE, numCores=NULL, ...) {
# Used for testing specific hypotheses, whether between clades, or against background
    phyData <- prepareData(phy=phy, richness=richness, verbose=verbose);
    phy <- phyData$phy;
    richness <- phyData$richness;
    
    conf <- configureModel(model=model, epsilon=epsilon, r=r, b=b, d=d,
        initialR=initialR, initialE=initialE);
    sp <- conf$sp;
    model <- conf$model;
    fixPar <- conf$fixPar;
    
## Keep track of all nodes, internal and pendant (for keeping track of breakpoints)
    pend.nodes <- seq_len(length(phy$tip.label));   # Calculate pendant splits just once, keep track through various models
    int.nodes <- (length(phy$tip.label)+2):max(phy$edge); # Omit root node
    root.node <- length(phy$tip.label) + 1;
    all.nodes <- c(pend.nodes, root.node, int.nodes);
    
    ## The important bits. Set up z, get descendants and number of tips per node
    obj <- makeCacheMedusa(phy=phy, richness=richness, all.nodes=all.nodes, mc=mc, numCores=numCores);
    
    desc <- list(stem=obj$desc.stem, node=obj$desc.node);
    z <- obj$z;

    fit <- medusaMLBaseModel(z=z, extractClade=extractClade, nodes=userNode, model=model,
        shiftCut=shiftCut, desc=desc, criterion=criterion, sp=sp, fixPar=fixPar);
    
    models <- list(fit);    
    
    for (i in 1:length(userNode)) { # ugh. get rid of loop and use *apply instead
        fit <- list(medusaFitUserNode(userNode[i], z, desc, model, root.node, criterion, shiftCut, sp=sp, fixPar));
        models <- c(models, fit);
    }
    results <- list(desc=desc, models=models, phy=phy);
    
# *** Need new summary function ***
    
    
    
    
    
    
    
#    return(models);
}

medusaFitUserNode <- function (node, z, desc, model, root.node, criterion, shiftCut, sp, fixPar) {
## various combinations possible
    fit.stem <- NULL;
    fit.node <- NULL;
    cut.at <- NULL;
    
    if (shiftCut == "stem" | shiftCut == "both") {
        obj.stem <- medusaSplitStem(node=node, z=z, desc=desc$stem);
        z.stem <- obj.stem$z[obj.stem$z[,"partition"] == 2,,drop=FALSE];
        
        fit.stem <- getOptimalModelFlavour(z=z.stem, sp=sp, model=model, fixPar=fixPar, criterion=criterion);
    }
    if ((shiftCut == "node" || shiftCut == "both") && (node > root.node)) {
        obj.node <- medusaSplitNode(node=node, z=z, desc=desc$node);
        z.node <- obj.node$z[obj.node$z[,"partition"] == 2,,drop=FALSE];
        
        fit.node <- getOptimalModelFlavour(z=z.node, sp=sp, model=model, fixPar=fixPar, criterion=criterion);
    }
## Now, figure out which shift position is optimal    
    if (is.null(fit.node)) {
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
        
        if (stem.fit[[element]] < node.fit[[element]]) {
            fit <- fit.stem;
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
        aicc=round(model.fit[2], digits=7), num.par=model.fit[3], cut.at=fit$cut.at, model=fit$model));
}

## Only used for base model
medusaMLBaseModel <- function (z, extractClade, nodes, model, shiftCut, desc, criterion, sp, fixPar) {
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
    
    z.base <- z[z[,"partition"] == 1,,drop=FALSE];
    fit <- getOptimalModelFlavour(z=z.base, sp=sp, model=model, fixPar=fixPar, criterion=criterion);
    model.fit <- calculateModelFit(fit=fit, z=z.base);
    
    return(list(par=matrix(fit$par, nrow=1, dimnames=list(NULL,c("r", "epsilon"))), lnLik=fit$lnLik,
        split.at=rootnode, aic=round(model.fit[1], digits=7), aicc=round(model.fit[2], digits=7), num.par=model.fit[3],
        cut.at="node", model=model.id));
}