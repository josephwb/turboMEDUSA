rm(list=ls())

## WORKING with turboMEDUSA 0.92-3

## PLOTTING FUNCTIONS ##
add.transparency=function (col, alpha) 
{
	tmp <- col2rgb(col)/255
	rgb(tmp[1, ], tmp[2, ], tmp[3, ], alpha = alpha)
}

plot.turboMEDUSA=function (results, model = NULL, cutoff = "threshold", criterion = "aicc", nd.cex=0.5, edge.transp=0.8, nd.transp=1, cex=0.5, param="r", legend=TRUE, sigdig=2, ...) 
{
#	model = NULL; cutoff = "threshold"; criterion = "aicc"; nd.cex=0.5; cex=0.5; param="r"; legend=TRUE; sigdig=2; edge.transp=0.8; nd.transp=0.9
	require(colorspace)
	plotTree=TRUE; time = FALSE; node.labels = TRUE;     
    fit <- results$models
    phy <- results$phy
    z <- results$z
    anc <- results$desc
    model.summary <- results$modelSummary
    threshold <- results$threshold
    model.id <- 0
    if (!is.null(model)) {
        model.id <- model
    } else {
        if (cutoff != "threshold") {
            threshold <- cutoff
        } else {
            cat("\nSelecting model based on corrected threshold (decrease in information theoretic score of ", 
                -threshold, " units).\n", sep = "")
        }
        model.id <- 1
        while (1) {
            if ((model.id + 1) > length(fit)) 
			break
            if ((unlist(fit[[model.id]][criterion]) - unlist(fit[[model.id + 
															 1]][criterion])) < threshold) 
			break
            model.id <- model.id + 1
        }
    }
    break.pts <- fit[[model.id]]$split.at
    opt.model <- as.data.frame(cbind(Split.node = break.pts, 
									 fit[[model.id]]$par, LnLik.part = fit[[model.id]]$lnLik.part))
    base.model <- as.data.frame(fit[[1]]$par)
    cat("\nEstimated parameter values for model #", model.id, 
        ":\n\n", sep = "")
    print.data.frame(opt.model, digits = 5)
    opt.weight <- 0
    opt.fit <- 0
    base.weight <- 0
    base.fit <- 0
    if (criterion == "aicc") {
        opt.weight <- model.summary$w.aicc[model.id]
        opt.fit <- model.summary$aicc[model.id]
        base.weight <- model.summary$w.aicc[1]
        base.fit <- model.summary$aicc[1]
    } else {
        opt.weight <- model.summary$w.aic[model.id]
        opt.fit <- model.summary$aic[model.id]
        base.weight <- model.summary$w.aic[1]
        base.fit <- model.summary$aic[1]
    }
    cat("\nModel fit summary for model #", model.id, ":\n\n", 
        sep = "")
    cat("\tLog-likelihood = ", as.numeric(results$models[[model.id]]["lnLik"]), 
        "\n", sep = "")
    cat("\t", criterion, " = ", opt.fit, "\n", sep = "")
    cat("\t", criterion, " weight = ", opt.weight, "\n\n", sep = "")
    if (model.id != 1) {
        cat("\nFor comparison, estimated values for the base (homogeneous-BD) model are:\n\n")
        print.data.frame(base.model, digits = 5, row.names = FALSE)
        cat("\nModel fit summary for base model:\n\n", sep = "")
        cat("\tLog-likelihood = ", as.numeric(results$models[[1]]["lnLik"]), 
            "\n", sep = "")
        cat("\t", criterion, " = ", base.fit, "\n", sep = "")
        cat("\t", criterion, " weight = ", base.weight, "\n\n", 
            sep = "")
    }
	z[,"partition"][]=1
    if (length(break.pts) > 1) {
        for (i in 2:length(break.pts)) {
            tmp <- medusa.split(break.pts[i], z, anc, "stem")
            z <- tmp$z
        }
    }
    if (plotTree) {
        col<-col.orig<-diverge_hcl(nrow(opt.model))
        
		if(param=="r") {
			rr=sort(opt.model$r)
			col=col[match(opt.model$r, rr)]
		} else {
			tmp=opt.model$epsilon
			tmp[is.na(tmp)]=0
			rr=sort(tmp)
			col=col[match(tmp, rr)]
		} 
        margin <- FALSE
        mm <- match(phy$edge[, 2], z[, "dec"])
        if (time) {
            margin = TRUE
        }
		
		shift.colors=diverge_hcl(3)
		shift.colors[2]="black"
		breakpt.colors=character(length(break.pts))
		
		ests=results$models[[model.id]]$par
		ests[is.na(ests)]=0
		for(i in 1:length(break.pts)){
			
			if(i==1){
				breakpt.colors[i]=shift.colors[2]
			} else {
				cur_part=as.numeric(z[z[,"dec"]==break.pts[i],"partition"])
				
				cur_parm=ests[cur_part,param]
				anc=get.ancestor.of.node(break.pts[i], phy)
				anc_part=unique(z[z[,"dec"]==anc,"partition"])
				anc_parm=ests[anc_part,param]
				if(cur_parm>anc_parm) {
					breakpt.colors[i]=shift.colors[3] 
				} else if(cur_parm<anc_parm) {
					breakpt.colors[i]=shift.colors[1]
				} else if(cur_parm==anc_parm){
					breakpt.colors[i]=shift.colors[2]
				} else {
					stop()
				}
				
			}
		}
				
		## TREE PLOTTING 
		diversitree:::plot.clade.tree(phy, edge.color = add.transparency(col[z[mm, "partition"]], edge.transp), no.margin = !margin, cex = cex, type="fan", ...)
		
		if(legend){
			sub.rr=pretty(range(rr), n=7)
			col.rr=col.orig[seq(1,length(col.orig), length=length(sub.rr))]
			legend(x="topleft", legend=sprintf(paste("%.",sigdig,"f",sep=""),sub.rr), bty="n", pch=22, pt.bg=col.rr, col=col.rr)
		}
        if (time) {
            axisPhylo(cex.axis = 0.75)
            mtext("Divergence Time (MYA)", at = (max(get("last_plot.phylo", 
														 envir = .PlotPhyloEnv)$xx) * 0.5), side = 1, 
				  line = 2, cex = 0.75)
        }
        if (node.labels) {
			if(length(break.pts)>1){
				for (i in 2:length(break.pts)) {
					node.scale=sqrt((length(break.pts)-i)/(length(break.pts)))*nd.cex+nd.cex
					if(break.pts[i]<=Ntip(phy)){
## FIXME: deal better with edgelabels on fan plot
#						sel=which(phy$edge[,2]==break.pts[i])
#						edgelabels.auteur(text=NULL, edge = sel, frame = "c", font = 1, cex = node.scale, col=add.transparency(shift.colors[2], nd.transp), bg=add.transparency(breakpt.colors[i], nd.transp), pch=21)
#						edgelabels.auteur(text=i, edge = sel, frame = "n", font = 1, cex = node.scale/4, col=add.transparency("white",nd.transp))
						nodelabels(text=NULL, node = break.pts[i], frame = "c", font = 1, cex = node.scale, col=add.transparency(shift.colors[2], nd.transp), bg=add.transparency(breakpt.colors[i], nd.transp), pch=21)
						nodelabels(text=i, node = break.pts[i], frame = "n", font = 1, cex = node.scale/4, col=add.transparency("white",nd.transp))
					} else {
						nodelabels(text=NULL, node = break.pts[i], frame = "c", font = 1, cex = node.scale, col=add.transparency(shift.colors[2], nd.transp), bg=add.transparency(breakpt.colors[i], nd.transp), pch=21)
						nodelabels(text=i, node = break.pts[i], frame = "n", font = 1, cex = node.scale/4, col=add.transparency("white",nd.transp))
					}
				}
			}
		}
    }
	return(results$models[[model.id]])
}

summaryplot.turboMEDUSA=function (results, shifts, model = NULL, cutoff = "threshold", criterion = "aicc", nd.cex=0.5, edge.transp=0.8, nd.transp=0.7, cex=0.5, param="r", legend=TRUE, sigdig=2, ...) 
{
#	model = NULL; cutoff = "threshold"; criterion = "aicc"; nd.cex=0.5; cex=0.5; param="r"; legend=TRUE; sigdig=2;
#   shifts: a named vector of shift directions (-1 to 1), summarized across a distribution of trees
#   results: summary from a single analysis (e.g., the MLE tree) 
	require(colorspace)
	plotTree=TRUE; time = FALSE; node.labels = TRUE;     
    fit <- results$models
    phy <- results$phy
    z <- results$z
    anc <- results$desc
    model.summary <- results$model.summary
    threshold <- results$threshold
    model.id <- 0
    if (!is.null(model)) {
        model.id <- model
    } else {
        if (cutoff != "threshold") {
            threshold <- cutoff
        } 
        model.id <- 1
        while (1) {
            if ((model.id + 1) > length(fit)) 
			break
            if ((unlist(fit[[model.id]][criterion]) - unlist(fit[[model.id + 
															 1]][criterion])) < threshold) 
			break
            model.id <- model.id + 1
        }
    }
    break.pts <- fit[[model.id]]$split.at
    opt.model <- as.data.frame(cbind(Split.node = break.pts, 
									 fit[[model.id]]$par, LnLik.part = fit[[model.id]]$lnLik.part))
    base.model <- as.data.frame(fit[[1]]$par)
    print.data.frame(opt.model, digits = 5)
    opt.weight <- 0
    opt.fit <- 0
    base.weight <- 0
    base.fit <- 0
    if (criterion == "aicc") {
        opt.weight <- model.summary$w.aicc[model.id]
        opt.fit <- model.summary$aicc[model.id]
        base.weight <- model.summary$w.aicc[1]
        base.fit <- model.summary$aicc[1]
    } else {
        opt.weight <- model.summary$w.aic[model.id]
        opt.fit <- model.summary$aic[model.id]
        base.weight <- model.summary$w.aic[1]
        base.fit <- model.summary$aic[1]
    }
	z[,"partition"][]=1
    if (length(break.pts) > 1) {
        for (i in 2:length(break.pts)) {
            tmp <- medusa.split(break.pts[i], z, anc, "stem")
            z <- tmp$z
        }
    }
    if (plotTree) {
        col<-col.orig<-diverge_hcl(nrow(opt.model))
        
		ests=opt.model
		ests[is.na(ests)]=0

		if(param=="r") {
			rr=sort(ests$r)
			col=col[match(ests$r, rr)]
		} else {
			rr=sort(ests$epsilon)
			col=col[match(ests$epsilon, rr)]
		} 
        margin <- FALSE
        mm <- match(phy$edge[, 2], z[, "dec"])
        if (time) {
            margin = TRUE
        }
		
		shift.colors=diverge_hcl(3)
		shift.colors[2]="black"
		breakpt.colors=character(length(break.pts))
		
		## TREE PLOTTING 
		diversitree:::plot.clade.tree(phy, edge.color = add.transparency(col[z[mm, "partition"]], edge.transp), no.margin = !margin, cex = cex, type="fan", ...)
		
		if(legend){
			sub.rr=pretty(range(rr), n=7)
			col.rr=col.orig[seq(1,length(col.orig), length=length(sub.rr))]
			legend(x="topleft", legend=sprintf(paste("%.",sigdig,"f",sep=""),sub.rr), bty="n", pch=22, pt.bg=col.rr, col=col.rr)
			legend(x="center", legend=sprintf(paste("%.",sigdig,"f",sep=""),seq(0,1,length=3)), pt.cex=seq(0,nd.cex,length=3), bty="n", pch=21, pt.bg=add.transparency("white",0.2), col=1)
		}
        if (time) {
            axisPhylo(cex.axis = 0.75)
            mtext("Divergence Time (MYA)", at = (max(get("last_plot.phylo", 
														 envir = .PlotPhyloEnv)$xx) * 0.5), side = 1, 
				  line = 2, cex = 0.75)
        }
        if (node.labels) {
			all_labels=c(phy$tip.label, phy$node.label)
			for (i in 1:length(shifts)) {
				nd=which(all_labels==names(shifts[i]))
				nd_col=ifelse(shifts[i]==0, shift.colors[2], ifelse(shifts[i]>0, shift.colors[3], ifelse(shifts[i]<0, shift.colors[1], shift.colors[2])))
				nodelabels(text=NULL, node = nd, frame = "c", font = 1, cex = nd.cex*abs(shifts[i]), col=add.transparency(shift.colors[2], nd.transp), bg=add.transparency(nd_col, nd.transp), pch=21)
			}
		}
    }
	return(results$models[[model.id]])
}

## END PLOTTING FUNCTIONS ##

## turboMEDUSA SUPPLEMENTARY FUNCTIONS ##

results.turboMEDUSA=function (results, model = NULL, cutoff = "threshold", criterion = "aicc", nd.cex=0.5, nd.transp=1, cex=0.5, sigdig=2, param="r", ...) 
{
#	model = NULL; cutoff = "threshold"; criterion = "aicc"; nd.cex=0.5; nd.transp=1; cex=0.5; sigdig=2; param="r";
#	results: come from runTurboMEDUSA
	require(colorspace)
	time = FALSE; node.labels = TRUE;     
    fit <- results$models
    phy <- results$phy
    z <- results$z
    anc <- results$desc
    model.summary <- results$modelSummary
    threshold <- results$threshold
    model.id <- 0
    if (!is.null(model)) {
        model.id <- model
    } else {
        if (cutoff != "threshold") {
            threshold <- cutoff
        } 
        model.id <- 1
        while (1) {
            if ((model.id + 1) > length(fit)) 
                break
            if ((unlist(fit[[model.id]][criterion]) - unlist(fit[[model.id + 
                1]][criterion])) < threshold) 
                break
            model.id <- model.id + 1
        }
    }
    break.pts <- fit[[model.id]]$split.at
    opt.model <- as.data.frame(cbind(Split.node = break.pts, 
        fit[[model.id]]$par, LnLik.part = fit[[model.id]]$lnLik.part))
    base.model <- as.data.frame(fit[[1]]$par)
    opt.weight <- 0
    opt.fit <- 0
    base.weight <- 0
    base.fit <- 0
    if (criterion == "aicc") {
        opt.weight <- model.summary$w.aicc[model.id]
        opt.fit <- model.summary$aicc[model.id]
        base.weight <- model.summary$w.aicc[1]
        base.fit <- model.summary$aicc[1]
    } else {
        opt.weight <- model.summary$w.aic[model.id]
        opt.fit <- model.summary$aic[model.id]
        base.weight <- model.summary$w.aic[1]
        base.fit <- model.summary$aic[1]
    }
	z[,"partition"][]=1
    if (length(break.pts) > 1) {
        for (i in 2:length(break.pts)) {
            tmp <- medusa.split(break.pts[i], z, anc, model.summary$Cut.at[i])
            z <- tmp$z
        }
    }
	
	shift.colors=diverge_hcl(3)
	breakpt.colors=character(length(break.pts))
	breakpt.trend=numeric(length(break.pts))
	ancestral_r=numeric(length(break.pts))
	ancestral_epsilon=numeric(length(break.pts))
	
	ests=results$models[[model.id]]$par
	ests[is.na(ests)]=0
	r_est=ests[z[,"partition"],"r"]
	epsilon_est=ests[z[,"partition"],"epsilon"]

	newz=cbind(z, r=r_est, epsilon=epsilon_est)
	
	for(param in c("r","epsilon")){
		for(i in 1:length(break.pts)){
			
			if(i==1){
				breakpt.colors[i]=shift.colors[2]
				breakpt.trend[i]=0
				if(param=="r") ancestral_r[i]=NA else ancestral_epsilon[i]=NA
				
				
			} else {
				cur_part=as.numeric(z[z[,"dec"]==break.pts[i],"partition"])
				
				cur_parm=ests[cur_part,param]
				anc=get.ancestor.of.node(break.pts[i], phy)
				anc_part=unique(z[z[,"dec"]==anc,"partition"])
				anc_parm=ests[anc_part,param]
				if(param=="r") ancestral_r[i]=anc_parm else ancestral_epsilon[i]=anc_parm
				
				if(cur_parm>anc_parm) {
					breakpt.colors[i]=shift.colors[3] 
					breakpt.trend[i]=1
				} else if(cur_parm<anc_parm){
					breakpt.colors[i]=shift.colors[1]
					breakpt.trend[i]=-1
				} else if(cur_parm==anc_parm){
					breakpt.colors[i]=shift.colors[2]
					breakpt.trend[i]=0
				} else {
					stop()
				}
			}
		}
		
		tmp=results$modelSummary
		tmp$trend=breakpt.trend
		names(tmp)[ncol(tmp)]=paste(param, "trend", sep="_")
		results$modelSummary=tmp		
	}
	
	results$modelSummary$ancestral_r=ancestral_r
	results$modelSummary$ancestral_epsilon=ancestral_epsilon

	
	fix_splits=function(phy, split.at){
		split.at[is.na(split.at)]=Ntip(phy)+1
		xx=split.at
		all=c(phy$tip.label, phy$node.label)
		names(xx)=all[xx]
		df=data.frame(name=names(xx),node=xx,stringsAsFactors=FALSE)
		rownames(df)=1:nrow(df)
		return(df)
	}
	
	results$model.summary=cbind(results$modelSummary, fix_splits(phy, results$modelSummary$Break.Node))
	results$model.summary=cbind(results$model.summary, results$models[[model.id]]$par)

	collect_times=function(phy, results, z){
		N=Ntip(phy)
		res=apply(results, 1, function(x){
			 cur=as.numeric(x[["Break.Node"]])
			 type=x[["Cut.at"]]
			 if(is.na(cur)) return(max(z[,"t.0"]))
			  if(cur<=N){
				return(z[which(z[,"dec"]==cur),"t.0"])
			  } else {
				if(type=="node"){
					return(z[which(z[,"dec"]==cur),"t.1"])
				} else {
					return(z[which(z[,"dec"]==cur),"t.0"])
				}
			  }
		})
		return(as.numeric(res))
		
	}
	
	results$model.summary$time=collect_times(phy, results$model.summary, z)
	
	return(list(z=newz, summary=results$model.summary))
}

fix_hash=function(phy, edges){
	tmp=store_edges(phy, edges)$label
	hash=attributes(edges)$hash
	nn=(Ntip(phy)+1):(Ntip(phy)+Nnode(phy))
	nd=hash[tmp[nn]]
	phy$hash=nd
	phy
}

turbomedusa.rank=function(phy, taxonomies=list(), edges=NULL, rank=c("genus","family","order"), ...){
# phy: ultrametric tree (tips should be at least genus level)
# richness: named vector of richnesses (names are tips in phy)
# taxonomies: 'phy' and 'richness' linkage between tips of and higher clades (e.g., family, order, etc)
	
	require(phylo)
	require(turboMEDUSA)
	require(diversitree)
	
	tax=get(data(spermatophyta_APGIII_taxonomy))
	fam=get(data(spermatophyta_APGIII_families))
	ord=get(data(spermatophyta_APGIII_orders))
	tmp=unlist(sapply(1:length(ord), function(idx) {y=ord[[idx]]; names(y)=rep(names(ord)[idx], length(y)); return(y)}))
	ord=data.frame(family=tmp, order=names(tmp), stringsAsFactors=FALSE)
	
	# prune tree
	phy_tax=taxonomies$phy
	unplc=rownames(phy_tax)[is.na(phy_tax[,rank])]
	pruned=drop.tip(phy, unplc)
	tre=prune_to_rank.phylo(pruned, phy_tax[rank], rank=rank)
	
	## RICHNESSES 
	richness_tax=taxonomies$richness
	tmp=sapply(split(richness_tax$richness, richness_tax[,rank]),sum)
	
	## RUN MEDUSA
	dat=data.frame(taxon=names(tmp), n.taxa=tmp, stringsAsFactors=FALSE)
	tmp=prune.tree.merge.data(tre, dat, FALSE)
	prn.phy=tmp$phy
	prn.dat=tmp$richness
	res=runTurboMEDUSA(prn.phy, prn.dat, stop="threshold", ...)
	
	# build cladetree to display results
	clades=lapply(prn.phy$tip.label, function(x) {
				  ntax=prn.dat$n.taxa[which(prn.dat$taxon==x)]
				  paste(x, 1:ntax, sep="_")
				  })
	names(clades)=prn.phy$tip.label
	
	cladetree=make.clade.tree(prn.phy, clades)
	
	# build lookup table for clade tree
	cld=get(data(spermatophyta_APGIII_clades))
	tmp=as.data.frame(tax[match(prn.phy$tip.label,tax[,rank]),which(names(tax)==rank):ncol(tax)])
	if(ncol(tmp)==1) prn.tax=NULL else prn.tax=tmp
	lkp=build_lookup.phylo(cladetree, prn.tax, cld)
	cladetree=nodelabel.phylo(cladetree, lkp, rank=NULL)
	
	if(is.null(edges)){
		edges=hash_edges(cladetree)
	}
	
	hash_phy=fix_hash(cladetree, edges)
	cladetree$node.label[cladetree$node.label==""]=hash_phy$hash[cladetree$node.label==""]
	
	res$phy=ladderize(cladetree, right=FALSE)
	res$edges=edges
	res
}

shifts_timeline=function(dat, root, dt=0.1, abs=TRUE){
	tt=seq(0,root,by=dt)
	yy=numeric(length(tt))
	cat(paste("using ", names(dat)[2], " as shift variable\n", sep=""))
	for(i in 1:nrow(dat)){
		sft=dat[i,2]
		time=dat$time[i]
		ll=abs(tt-time)
		ww=min(which(ll==min(ll)))
		if(abs) yy[ww]=yy[ww]+abs(sft)	else yy[ww]=yy[ww]+sft
	}
	return(cbind(tt,yy))	
}

## END turboMEDUSA SUPPLEMENTARY FUNCTIONS ##

## PROCESSING ## 
batch=TRUE
require(phylo)
require(multicore)
require(diversitree)
require(turboMEDUSA)
base="spermatophyta_AToL_639_PL_MEDUSA_BATCH"
cores=12

## FIND MLE ESTIMATE
# collect taxonomic information
pp=read.csv(paste(base,"MLE_0.92.phy-taxonomy.csv",sep="."), row=1, stringsAsFactors=FALSE)
rr=read.csv(paste(base,"MLE_0.92.richness-taxonomy.csv",sep="."), row=1, stringsAsFactors=FALSE)
taxonomies=list(phy=pp, richness=rr)

# collect MLE tree
ml=read.tree("out_dates.tre")
drop_genera=list(Escalloniaceae="Platyspermation", Olacaceae="Ximenia", Phytolaccaceae="Rivina", Gesneriaceae="Saintpaulia", Aristolochiaceae="Saruma")
mle=drop.tip(ml, unlist(drop_genera))

mle_result<-mle_res<-turbomedusa.rank(mle, taxonomies, edges=NULL, rank="family")
tmp=results.turboMEDUSA(mle_result)
mle_res$z=tmp$z
mle_res$result=tmp$summary
save(mle_res, file=paste(base, "familial_MLE.rda", sep="."))

pdf(paste(base, "familial_MLE.pdf", sep="."), height=60, width=60)

	hashes=attributes(mle_res$edges)$hash

	plot.turboMEDUSA(mle_res, nd.transp=0.5, nd.cex=10, transform=sqrt, param="r", cex=2)
	nodelabels(mle_res$phy$node.label, frame="n", cex=ifelse(mle_res$phy$node.label%in%hashes, 0.01, 3.5), col=add.transparency(1,0.5))
	mtext("net diversification")

	plot.turboMEDUSA(mle_res, nd.transp=0.5, nd.cex=10, transform=sqrt, param="epsilon", cex=2)
	nodelabels(mle_res$phy$node.label, frame="n", cex=ifelse(mle_res$phy$node.label%in%hashes, 0.01, 3.5), col=add.transparency(1,0.5))
	mtext("relative extinction")

	plot.turboMEDUSA(mle_res, nd.transp=0.5, nd.cex=10, transform=sqrt, param="r", cex=2, show.tip=FALSE)
	nodelabels(mle_res$phy$node.label, frame="n", cex=ifelse(mle_res$phy$node.label%in%hashes, 0.01, 3.5), col=add.transparency(1,0.5))
	mtext("net diversification")

	plot.turboMEDUSA(mle_res, nd.transp=0.5, nd.cex=10, transform=sqrt, param="epsilon", cex=2, show.tip=FALSE)
	nodelabels(mle_res$phy$node.label, frame="n", cex=ifelse(mle_res$phy$node.label%in%hashes, 0.01, 3.5), col=add.transparency(1,0.5))
	mtext("relative extinction")

dev.off()



if(batch){
	## RUN PL REPLICATES
	resdir="MEDUSA_batch"
	dir.create(resdir)
	phy=read.tree("bootstrap.clocking.tre")
	edges=mle_res$edges
	all_labs=c(mle_res$phy$tip.label, mle_res$phy$node.label)
	curbase=paste(base,"familial_bs",sep="_")
	ss=seq(1,length(phy),by=cores)
	for(b in ss){
		reps=b:(min(b+cores-1, length(phy)))
		x=mclapply(reps, function(idx){ 
				   tt=drop.tip(phy[[idx]], unlist(drop_genera))
				   if(!is.ultrametric(tt)) tt=ultrametricize(tt, 1, "mean")
				   cur=turbomedusa.rank(tt, taxonomies, edges=mle_res$edges, rank="family")
				   res=results.turboMEDUSA(cur)
				   if(!all(res$summary$name%in%all_labs)) cat(paste("Node labels not recognized in replicate: ", idx, sep=""))
				   save(res, file=paste(resdir, paste(idx, "rda", sep="."), sep="/"))
				   return(res)
				   }, mc.cores=cores)
		x=NULL
	}
	
	
	## COMPILE PL REPLICATES
	require(phylo)
	require(turboMEDUSA)
	resdir="MEDUSA_batch"
	base="spermatophyta_AToL_639_PL_MEDUSA_BATCH"
	curbase=paste(base,"familial_MLE",sep=".")
	mle_res=get(load(paste(curbase,"rda",sep=".")))
	mle=mle_res$phy
	
	resfiles=dir(resdir, full.names=TRUE)
	
	r.res<-e.res<-t.res<-matrix(0, nrow=Ntip(mle)+Nnode(mle), ncol=length(resfiles))
	rownames(r.res)<-rownames(e.res)<-rownames(t.res)<-rownm<-c(mle$tip.label, mle$node.label)
	colnames(r.res)<-colnames(e.res)<-colnames(t.res)<-basename.noext(resfiles)
	
	## compile results from PL replicates
	for(i in 1:length(resfiles)){
		
		cat("\n\n",i,"\n\n")
		
		tmp=get(load(resfiles[i]))
		dat=tmp$summary
		dat$epsilon[is.na(dat$epsilon)]=0
		z=tmp$z
		
		new=dat[-which(is.na(dat$Break.Node)),]
		new$r_shift=new$r-new$ancestral_r
		new$epsilon_shift=new$epsilon-new$ancestral_epsilon
		
		mm=match(new$name,rownm)
		t.res[mm,i]=new$time
		r.res[mm,i]=new$r_shift
		e.res[mm,i]=new$epsilon_shift
	}
	
	batch.res=list(e=e.res, r=r.res, time=t.res)
	save(batch.res, file=paste(base, "effectsize", "rda", sep="."))
	
	# compile PL replicate results for mean shift time, effect sizes, and weight of evidence
	shift_pct=apply(t.res, 1, function(x) sum(x>0))/ncol(t.res)
	tmean=sapply(1:nrow(t.res), function(idx) {x=t.res[idx,]; y=mean(x[x>0]); if(is.na(y)) return(0) else return(y)})
	rmean=sapply(1:nrow(r.res), function(idx) {x=r.res[idx,]; y=mean(x[x!=0]); if(is.na(y)) return(0) else return(y)})
	emean=sapply(1:nrow(e.res), function(idx) {x=e.res[idx,]; y=mean(x[x!=0]); if(is.na(y)) return(0) else return(y)})
	mean_res=as.data.frame(cbind(pct=shift_pct, time=tmean, r=rmean, epsilon=emean), stringAsFactors=FALSE)
	
	cutoff=0.75
	root=max(mle_res$z[,"t.0"])
	
	rs=shifts_timeline(mean_res[mean_res$pct>=cutoff,c("time", "r")], root, dt=1, abs=FALSE)
	es=shifts_timeline(mean_res[mean_res$pct>=cutoff,c("time", "epsilon")], root, dt=1, abs=FALSE)
	
	pdf("spermatophyta_AToL_639_PL_MEDUSA_BATCH.shifteffects_BATCH.pdf", width=24, height=12)
		# PLOT MAGNITUDE of shifts for shifts occuring in at least cutoff percent of tree replicates
		# contrast is shifted parameter value - direct ancestral value (positive values are increases in diversification parameter) 
		mx=max(c(abs(es[,2]),abs(rs[,2])))
		plot(rs[,1],rs[,2],type="l", xlim=c(root,0), xlab="MYA", ylab="magnitude of shift", bty="n", col=add.transparency(1,0.75), ylim=c(-mx,mx), lwd=2.5)
		lines(es[,1],es[,2],type="l", xlim=c(root,0), bty="n", col=add.transparency("gray",0.75), lwd=2.5)
		legend("topleft", legend=c("net diversification","relative extinction"), pt.bg=c(add.transparency(1,0.75),add.transparency("gray",0.75)), pch=22, cex=2, bty="n")
	
	
		# PLOT MAGNITUDE of shifts SEPARATELY for shifts occuring in at least cutoff percent of tree replicates
		# contrast is shifted parameter value - direct ancestral value (positive values are increases in diversification parameter) 
		plot(rs[,1],rs[,2],type="l", xlim=c(root,0), xlab="MYA", ylab="magnitude of shift: net diversification", bty="n", col=add.transparency(1,0.75), lwd=2.5)
		plot(es[,1],es[,2],type="l", xlim=c(root,0), xlab="MYA", ylab="magnitude of shift: relative extinction", bty="n", col=add.transparency(1,0.75), lwd=2.5)
	
		# PLOT SHIFTS according to temporal uncertainty (if 1.0, same shift in same 1 MY bin across all tree replicates)
		tt=c(t.res[t.res>0])
		zz=numeric(length(tt))
		zz[]=1
		tmp=data.frame(time=tt, shift=zz)
		res=shifts_timeline(tmp, root, dt=1, abs=TRUE)
		res[,2]=res[,2]/ncol(t.res)
		plot(res[,1],res[,2],type="l", xlim=c(root,0), xlab="MYA", ylab="shifts", bty="n", lwd=2.5)
	
		## ADD PLOTS FOR PARAMETERS THROUGH TIME
	dev.off()
	
}


