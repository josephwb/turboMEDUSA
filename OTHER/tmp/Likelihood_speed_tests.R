system.time(runTurboMEDUSA(tmp$phy, modelLimit=100, shiftCut="both", model="mixed") -> res.mixed) -> time.mixed
system.time(runTurboMEDUSA(tmp$phy, modelLimit=100, shiftCut="both", model="bd") -> res.mixed) -> time.bd
system.time(runTurboMEDUSA(tmp$phy, modelLimit=100, shiftCut="both", model="yule") -> res.mixed) -> time.yule


medusa.fun <- function()
{
	obj <- make.cache.medusa(phy, richness, all.nodes, mc, num.cores)
	z <- obj$z
	new.part <- z[z[,"partition"] == 1,,drop=FALSE];
	for (i in 1:100000)
	{
		lik <- make.lik.medusa.part(partition=new.part, model=model);
		foo <- function (x) {-lik(pars=exp(x));}
		optimize(f=foo, interval=c(-25, 1));
	}
}

yule.fun <- function()
{
	for (i in 1:100000)
	{
		yule(phy);
	}
}

pureBirth.fun <- function()
{
	a <- branching.times(phy);
	for (i in 1:100000)
	{
		pureBirth(a);
	}
}

system.time(medusa.fun());
system.time(yule.fun());
system.time(pureBirth.fun());

> system.time(yule.fun());
   user  system elapsed 
 87.793   0.809  95.259 
 
 > system.time(medusa.fun());
   user  system elapsed 
 11.358   0.327  11.751 

> system.time(pureBirth.fun());
   user  system elapsed 
 54.278   0.352  54.400 


########################################################################

old MEDUSA PB Likelihood

    isInt <- is.na(z[,"n.t"])
    isTerm <- !isInt
    int<-z[isInt,]
    term<-z[isTerm,]

    nint <- sum(isInt)
    nterm <- sum(isTerm)
    
    lbetaPB <- function(b, aa, tt) {
    	res<-log((exp(b*tt)-1)/(exp(b*tt)))
        res
    }
	
	Lfunc_comb_pb <- function(r) {
    	if(r<0) return(-Inf)

    	b<-r
    	d<-0
    	eps<-0
    	
            branchLengths<-int[,"t.len"]
    		lint<-nint * log(r) - r * sum(branchLengths) - sum(log(1 - (eps * exp(-r * int[, 3]))))
    #		a<-term[,"startR"] # Eek! Same name as initial value of r
    #		n<-term[,"endR"]
    		
    		a<-term[,"n.0"]
    		n<-term[,"n.t"]
    		
    		timeInterval<-term[,"t.len"]
    		endj<-pmin(a, n)
    		sum<-0
    		lnl<-numeric(nterm)
    		for(i in 1:nterm) {
    			lnxx<-numeric(length=endj[i])
    			for(j in 1:endj[i]) { #using likelihoods from Foote et al. 1999, a correction of Raup
    				logAlpha<- -Inf
    				logBeta<-lbetaPB(b, a[i], timeInterval[i])
    				s1<-lchoose(a[i],j)+lchoose(n[i]-1,j-1)
    				
    				if(logAlpha==-Inf) s2<-0 else s2<-(a[i]-j)*logAlpha
    				s3<-log(((1-exp(logAlpha))*(1-exp(logBeta)))^j)
    				s4<-(n[i]-1)*logBeta
    				s5<-log(1-exp(logAlpha)) # Conditioning on survival to the present
       				lnxx[j]<-s1+s2+s3+s4-s5
    			}
    			lnl[i]<-logspace_sum(lnxx)
    		}
    		lterm<-sum(lnl)
    	return(lint+lterm)
	}
	
    res <- list()


    foo.old.medusa <- function(x) 
    		- Lfunc_comb_pb(r=exp(x[1]))
	sp<-log(startR)
	o<-optimize(foo, interval=c(-25, 1))
	res$LH <- -o$objective
   	res$par <- exp(o$minimum)

tester <- function (func) {
	for (i in 1:10000) {
		optimize(func, interval=c(-25, 1))
	}
}

o<-optimize(foo.old.medusa, interval=c(-25, 1));
res <- optimize(f= foo.turbomedusa, interval=c(-25, 1));
res <- optimize(f= foo.turbomedusa.2, interval=c(-25, 1));

system.time(tester(foo.old.medusa));
system.time(tester(foo.turbomedusa));
system.time(tester(foo.turbomedusa.2));
system.time(tester(foo.turbomedusa.3));


> system.time(tester(foo.old.medusa));
   user  system elapsed 
154.919   0.383 154.591 
> system.time(tester(foo.turbomedusa));
   user  system elapsed 
  0.639   0.002   0.666 
> system.time(tester(foo.turbomedusa.2));
   user  system elapsed 
  1.974   0.110   2.101 
> system.time(tester(foo.turbomedusa.3));
   user  system elapsed 
  2.830   0.117   2.941 


tip.tester <- function (model)
{
	for (i in 1:100)
	{
		tips <- NULL;
		tips <- lapply(pend.nodes, medusaMLPrefitStem, z=z, desc=desc$desc.stem, initialR=initialR, initialE=initialE,
				model=model, criterion=criterion);
	}
}

system.time(tip.tester(model="yule")) -> yule.tip.time; yule.tip.time;
system.time(tip.tester(model="bd")) -> bd.tip.time; bd.tip.time;
system.time(tip.tester(model="mixed")) -> mixed.tip.time; mixed.tip.time;


internal.tester <- function (model)
{
	for (i in 1:100)
	{
		virgin.stem <- list(); virgin.node <- list();
		if (shiftCut == "stem" || shiftCut == "both") {
			virgin.stem <- lapply(int.nodes, medusaMLPrefitStem, z=z, desc=desc$desc.stem, initialR=initialR, initialE=initialE,
				model=model, criterion=criterion);
		}
		if (shiftCut == "node" || shiftCut == "both") {
			virgin.node <- lapply(int.nodes, medusaMLPrefitNode, z=z, desc=desc$desc.node, initialR=initialR, initialE=initialE,
				model=model, criterion=criterion);
		}
	}
}

system.time(internal.tester(model="yule")) -> yule.int.time; yule.int.time;
system.time(internal.tester(model="bd")) -> bd.int.time; bd.int.time;
system.time(internal.tester(model="mixed")) -> mixed.int.time; mixed.int.time;



prefit.tester <- function (model)
{
	for (i in 1:100)
	{
		tips <- NULL;
		tips <- lapply(pend.nodes, medusaMLPrefitStem, z=z, desc=desc$desc.stem, initialR=initialR, initialE=initialE,
			model=model, criterion=criterion);
		
		virgin.stem <- list(); virgin.node <- list();
		if (shiftCut == "stem" || shiftCut == "both") {
			virgin.stem <- lapply(int.nodes, medusaMLPrefitStem, z=z, desc=desc$desc.stem, initialR=initialR, initialE=initialE,
				model=model, criterion=criterion);
		}
		if (shiftCut == "node" || shiftCut == "both") {
			virgin.node <- lapply(int.nodes, medusaMLPrefitNode, z=z, desc=desc$desc.node, initialR=initialR, initialE=initialE,
				model=model, criterion=criterion);
		}
	}
}

system.time(prefit.tester(model="yule")) -> yule.prefit.time; yule.prefit.time;
system.time(prefit.tester(model="bd")) -> bd.prefit.time; bd.prefit.time;
system.time(prefit.tester(model="mixed")) -> mixed.prefit.time; mixed.prefit.time;


base.tester <- function (model)
{
	for (i in 1:10000)
	{
		fit <- list();
		fit <- medusaMLFitBase(z=z, initialR=initialR, initialE=initialE, model=model, criterion=criterion);
		models <- list(fit);
	}
}

system.time(base.tester(model="yule")) -> yule.base.time; yule.base.time;
system.time(base.tester(model="bd")) -> bd.base.time; bd.base.time;
system.time(base.tester(model="mixed")) -> mixed.base.time; mixed.base.time;

node.list <- (1:(root.node - 1)) # just tips
node.list <- ((root.node + 1):max(all.nodes)) # just int
node.list <- all.nodes # all


update.tester <- function (model)
{
	fit <- list();
	fit <- medusaMLFitBase(z=z, initialR=initialR, initialE=initialE, model="bd", criterion=criterion);
	
	for (i in 1:10)
	{
		res <- lapply(node.list, medusaMLUpdate, z=z, desc=desc, fit=fit, prefit=prefit, num.tips=num.tips,
			root.node=root.node, model=model, criterion=criterion, shiftCut=shiftCut);
	}
}

gc();
system.time(update.tester(model="yule")) -> yule.update.time; yule.update.time;
system.time(update.tester(model="bd")) -> bd.update.time; bd.update.time;
system.time(update.tester(model="mixed")) -> mixed.update.time; mixed.update.time;

