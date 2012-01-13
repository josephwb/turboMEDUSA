makeLikMedusaPart <- function (partition, model, fixPar)
{
	i.int <- is.na(partition[,"n.t"])
	i.pend <- !(i.int)
	n.int <- sum(i.int);
	
	simple <- FALSE;
	if (all(partition[,"n.t"] == 1, na.rm = TRUE)) {simple <- TRUE;}
	
	if (simple) {
		if (model == "yule") {
			sum.t <- sum(partition[,"t.len"]);
			f <- function(pars) {
				if (pars <= 0) return(-Inf);
				r <- pars[1];
				return(n.int * log(r) - r * sum.t);
			}
			return(f)
		} 
		
		sum.int.t.len <- sum(partition[i.int,"t.len"]);
		int.t.0 <- partition[i.int,"t.0"];
		pend.t.len <- partition[i.pend,"t.len"];
		
		if (model == "bd") {
			f <- function(pars) {
				r <- pars[1];
				epsilon <- pars[2];
				if (r <= 0 | epsilon <= 0 | epsilon >= 1) {return(-Inf);}
				l.int <- n.int * log(r) - r * sum.int.t.len - sum(log(1 - (epsilon * exp(-r * int.t.0))));
				l.pend <- sum(log(1 - epsilon) - log(exp(r * pend.t.len) - epsilon));
				return(l.int + l.pend);
			}
			return(f)
		} else if (model == "fixedD") {
			fixD <- fixPar;
			f <- function(pars) {
				r <- pars[1] + fixD;
				epsilon <- fixD/pars[2];
				if (r <= 0 | epsilon <= 0 | epsilon >= 1) {return(-Inf);}
				l.int <- n.int * log(r) - r * sum.int.t.len - sum(log(1 - (epsilon * exp(-r * int.t.0))));
				l.pend <- sum(log(1 - epsilon) - log(exp(r * pend.t.len) - epsilon));
				return(l.int + l.pend);
			}
			return(f)
		} else if (model == "fixedEpsilon") {
			epsilon <- fixPar;
			f <- function(pars) {
				r <- pars[1];
				if (r <= 0 | epsilon <= 0 | epsilon >= 1) {return(-Inf);}
				l.int <- n.int * log(r) - r * sum.int.t.len - sum(log(1 - (epsilon * exp(-r * int.t.0))));
				l.pend <- sum(log(1 - epsilon) - log(exp(r * pend.t.len) - epsilon));
				return(l.int + l.pend);
			}
			return(f)
		} else if (model == "fixedR") { # can make this simpler; precalculate most things
			r <- fixPar;
			f <- function(pars) {
				epsilon <- pars[1];
				if (r <= 0 | epsilon <= 0 | epsilon >= 1) {return(-Inf);}
				l.int <- n.int * log(r) - r * sum.int.t.len - sum(log(1 - (epsilon * exp(-r * int.t.0))));
				l.pend <- sum(log(1 - epsilon) - log(exp(r * pend.t.len) - epsilon));
				return(l.int + l.pend);
			}
			return(f)
		} else if (model == "fixedB") {
			fixB <- fixPar;
			f <- function(pars) {
				r <- fixB - pars[1];
				epsilon <- fixB * pars[2];
				if (r <= 0 | epsilon <= 0 | epsilon >= 1) {return(-Inf);}
				l.int <- n.int * log(r) - r * sum.int.t.len - sum(log(1 - (epsilon * exp(-r * int.t.0))));
				l.pend <- sum(log(1 - epsilon) - log(exp(r * pend.t.len) - epsilon));
				return(l.int + l.pend);
			}
			return(f)
		}
	} else {
		sum.int.t.len <- sum(partition[i.int,"t.len"]);
		pend.t.len <- partition[i.pend,"t.len"];
		pend.n.t.minus.1 <- partition[i.pend,"n.t"] - 1;
		
		if (model == "yule") {
			sum.pend.t.len <- sum(pend.t.len);
			f <- function(pars) {
				r <- pars[1];
				if (r <= 0) return(-Inf);
				return(n.int * log(r) - r * sum.int.t.len + sum(-pend.t.len * r + pend.n.t.minus.1*log(1 - exp(-pend.t.len * r))));
			}
			return(f)
		}
		
		int.t.0 <- partition[i.int,"t.0"];
		
		 if (model == "bd") {
			f <- function(pars) {
				r <- pars[1];
				epsilon <- pars[2];
				if (r <= 0 | epsilon <= 0 | epsilon >= 1) {return(-Inf);}
				l.int <- n.int * log(r) - r * sum.int.t.len - sum(log(1 - (epsilon * exp(-r * int.t.0))));
				ert <- exp(r * pend.t.len);
				B <- (ert - 1) / (ert - epsilon);
				l.pend <- sum(log(1 - B) + pend.n.t.minus.1*log(B));
				return(l.int + l.pend);
			}
			return(f)
		} else if (model == "fixedD") {
			fixD <- fixPar;
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
		} else if (model == "fixedEpsilon") {
			epsilon <- fixPar;
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
		} else if (model == "fixedR") { # can make this simpler; precalculate most things
			r <- fixPar;
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
		} else if (model == "fixedB") {
			fixB <- fixPar;
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
}
