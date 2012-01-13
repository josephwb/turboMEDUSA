medusaMLFitBase(z=z, sp=sp, model="yule", criterion=criterion);
medusaMLFitBase(z=z, sp=sp, model="bd", criterion=criterion);
medusaMLFitBase(z=z, sp=sp, model="mixed", criterion=criterion);
medusaMLFitBase(z=z, sp=sp, model="fixedEpsilon", criterion=criterion);
medusaMLFitBase(z=z, sp=sp, model="fixedR", criterion=criterion);

medusaMLFitBase(z=z, sp=sp, model="fixedB", criterion=criterion);
medusaMLFitBase(z=z, sp=sp, model="fixedD", criterion=criterion);


## Only used for base model
medusaMLFitBase <- function (z, sp, model, criterion)
{
	fit <- getOptimalModelFlavour(partition.id=1, z=z, sp=sp, model=model, criterion=criterion);
	model.fit <- calculateModelFit(fit=fit, z=z);
	
	BD <- as.numeric(getBD(fit$par));
	
	return(list(REpsilon=matrix(fit$par, nrow=1, dimnames=list(NULL,c("r", "epsilon"))), 
		
		BD=matrix(BD, nrow=1, dimnames=list(NULL,c("birth", "death"))),
		
		lnLik=fit$lnLik, aicc=round(model.fit[2], digits=7), num.par=model.fit[3], model=fit$model));
}



mrca <- c("Musaceae", "Poaceae", "Smilacaceae")