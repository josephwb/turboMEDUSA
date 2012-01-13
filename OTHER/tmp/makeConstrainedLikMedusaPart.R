## makeLikMedusaPart: generate a likelihood function for a single partition.
makeConstrainedLikMedusaPart <- function (partition, model, sp)
{
	int  <- partition[is.na(partition[,"n.t"]),,drop=FALSE];
	pend <- partition[!is.na(partition[,"n.t"]),,drop=FALSE];
	n.int <- length(int[,1]);
	n.pend <- length(pend[,1]);
	sum.int.t.len <- sum(int[,"t.len"]);
	int.t.0 <- int[,"t.0"];
	pend.t.len <- pend[,"t.len"];
	
	simple <- FALSE;
	if (all(partition[,"n.t"] == 1, na.rm = TRUE)) {simple <- TRUE;}
		
	if (model == "fixedEpsilon" && simple) {
		epsilon <- sp[2];
		f <- function(pars) {
			r <- pars[1];
			if (r <= 0 | epsilon <= 0 | epsilon >= 1) {return(-Inf);}
			l.int <- n.int * log(r) - r * sum.int.t.len - sum(log(1 - (epsilon * exp(-r * int.t.0))));
			l.pend <- sum(log(1 - epsilon) - log(exp(r * pend.t.len) - epsilon));
			return(l.int + l.pend);
		}
		return(f)
	} else if (model == "fixedEpsilon" && !simple) {
		pend.n.t.minus.1 <- pend[,"n.t"] - 1;
		epsilon <- sp[2];
		f <- function(pars) {
			r <- pars[1];
			if (r <= 0 | epsilon <= 0 | epsilon >= 1) {return(-Inf);}
			l.int <- n.int * log(r) - r * sum.int.t.len - sum(log(1 - (epsilon * exp(-r * int.t.0))));
			ert <- exp(r * pend.t.len);
			B <- (ert - 1) / (ert - epsilon);
			l.pend <- sum(log(1 - B) + pend.n.t.minus.1*log(B));
			return(l.int + l.pend);
		}
		return(f)
	} else if (model == "fixedR" && simple) {
		r <- sp[1];
		f <- function(pars) {
			epsilon <- pars[1];
			if (r <= 0 | epsilon <= 0 | epsilon >= 1) {return(-Inf);}
			l.int <- n.int * log(r) - r * sum.int.t.len - sum(log(1 - (epsilon * exp(-r * int.t.0))));
			l.pend <- sum(log(1 - epsilon) - log(exp(r * pend.t.len) - epsilon));
			return(l.int + l.pend);
		}
		return(f)
	} else if (model == "fixedR" && !simple) {
		pend.n.t.minus.1 <- pend[,"n.t"] - 1;
		r <- sp[1];
		f <- function(pars) {
			epsilon <- pars[1];
			if (r <= 0 | epsilon <= 0 | epsilon >= 1) {return(-Inf);}
			l.int <- n.int * log(r) - r * sum.int.t.len - sum(log(1 - (epsilon * exp(-r * int.t.0))));
			ert <- exp(r * pend.t.len);
			B <- (ert - 1) / (ert - epsilon);
			l.pend <- sum(log(1 - B) + pend.n.t.minus.1*log(B));
			return(l.int + l.pend);
		}
		return(f)
	} else if (model == "fixedD" && simple) {
		fixD <- sp[2];
		f <- function(pars) {
			r <- pars[1] + fixD;
			epsilon <- fixD/pars[2];
			if (r <= 0 | epsilon <= 0 | epsilon >= 1) {return(-Inf);}
			l.int <- n.int * log(r) - r * sum.int.t.len - sum(log(1 - (epsilon * exp(-r * int.t.0))));
			l.pend <- sum(log(1 - epsilon) - log(exp(r * pend.t.len) - epsilon));
			return(l.int + l.pend);
		}
		return(f)
	}  else if (model == "fixedD" && !simple) {
		pend.n.t.minus.1 <- pend[,"n.t"] - 1;
		fixD <- sp[2];
		f <- function(pars) {
			r <- pars[1] + fixD;
			epsilon <- fixD/pars[2];
			if (r <= 0 | epsilon <= 0 | epsilon >= 1) {return(-Inf);}
			l.int <- n.int * log(r) - r * sum.int.t.len - sum(log(1 - (epsilon * exp(-r * int.t.0))));
			ert <- exp(r * pend.t.len);
			B <- (ert - 1) / (ert - epsilon);
			l.pend <- sum(log(1 - B) + pend.n.t.minus.1*log(B));
			return(l.int + l.pend);
		}
		return(f)
	} else if (model == "fixedB" && simple) {
		fixB <- sp[1];
		f <- function(pars) {
			r <- fixB - pars[1];
			epsilon <- fixB * pars[2];
			if (r <= 0 | epsilon <= 0 | epsilon >= 1) {return(-Inf);}
			l.int <- n.int * log(r) - r * sum.int.t.len - sum(log(1 - (epsilon * exp(-r * int.t.0))));
			l.pend <- sum(log(1 - epsilon) - log(exp(r * pend.t.len) - epsilon));
			return(l.int + l.pend);
		}
		return(f)
	} else if (model == "fixedB" && !simple) {
		pend.n.t.minus.1 <- pend[,"n.t"] - 1;
		fixB <- sp[1];
		f <- function(pars) {
			r <- fixB - pars[1];
			epsilon <- fixB * pars[2];
			if (r <= 0 | epsilon <= 0 | epsilon >= 1) {return(-Inf);}
			l.int <- n.int * log(r) - r * sum.int.t.len - sum(log(1 - (epsilon * exp(-r * int.t.0))));
			ert <- exp(r * pend.t.len);
			B <- (ert - 1) / (ert - epsilon);
			l.pend <- sum(log(1 - B) + pend.n.t.minus.1*log(B));
			return(l.int + l.pend);
		}
		return(f)
	}
}
