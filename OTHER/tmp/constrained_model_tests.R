fixB <- fixPar;
f.fb.1par <- function(pars) {
	d <- pars[1];
	r <- fixB - d;
	epsilon <- d/fixB;
	if (r <= 0 | epsilon <= 0 | epsilon >= 1) {return(-Inf);}
	l.int <- n.int * log(r) - r * sum.int.t.len - sum(log(1 - (epsilon * exp(-r * int.t.0))));
	l.pend <- sum(log(1 - epsilon) - log(exp(r * pend.t.len) - epsilon));
	return(l.int + l.pend);
}

lik <- f.fb.1par;
foo <- function (x) {-lik(pars=exp(x));}
fit <- optimize(f=foo, interval=c(-100, log(fixPar)));
d <- exp(fit$minimum);
par <- c((fixPar - d), (d / fixPar));
par;

########################################

fixB <- fixPar;
f.fb.2par.simple <- function(pars) {
#	cat("pars1=", pars[1], "; pars2 =", pars[2], "\n")
	r <- fixB - pars[1];
	epsilon <- fixB * pars[2];
#	epsilon <- pars[2]/fixB;
	cat("r =", r, "; epsilon=", epsilon, "\n")
	if (r <= 0 | epsilon <= 0 | epsilon >= 1 | is.na(epsilon)) {return(-Inf);}
	l.int <- n.int * log(r) - r * sum.int.t.len - sum(log(1 - (epsilon * exp(-r * int.t.0))));
	l.pend <- sum(log(1 - epsilon) - log(exp(r * pend.t.len) - epsilon));
	return(l.int + l.pend);
}

lik <- f.fb.2par;
foo <- function (x) {-lik(pars=exp(x));}
#sp = c((fixPar/2),fixPar/4); # make sure initial value of r evaluated is < b
sp = c((fixPar/2),0.001); # make sure initial value of r evaluated is < b
fit <- optim(fn=foo, par=log(sp), method="N", control=list(maxit=5000));
#par=c((fixPar - exp(fit$par[1])), (exp(fit$par[2])/fixPar));
par=c((fixPar - exp(fit$par[1])), (exp(fit$par[2])*fixPar));
par;

fixB <- fixPar;
f.fb.2par <- function(pars) {
	r <- fixB - pars[1];
	epsilon <- pars[2]/fixB;

	if (r <= 0 | epsilon <= 0 | epsilon >= 1) {return(-Inf);}
	l.int <- n.int * log(r) - r * sum.int.t.len - sum(log(1 - (epsilon * exp(-r * int.t.0))));
	ert <- exp(r * pend.t.len);
	B <- (ert - 1) / (ert - epsilon);
	l.pend <- sum(log(1 - B) + pend.n.t.minus.1*log(B));
	return(l.int + l.pend);
}

###########################################
		
f.bd <- function(pars) {
	r <- pars[1];
	epsilon <- pars[2];
	if (r <= 0 | epsilon <= 0 | epsilon >= 1) {return(-Inf);}
	l.int <- n.int * log(r) - r * sum.int.t.len - sum(log(1 - (epsilon * exp(-r * int.t.0))));
	ert <- exp(r * pend.t.len);
	B <- (ert - 1) / (ert - epsilon);
	l.pend <- sum(log(1 - B) + pend.n.t.minus.1*log(B));
	return(l.int + l.pend);
}

lik <- f.bd;
foo <- function (x) {-lik(pars=exp(x));}
sp=c(0.05, 0.5); # make sure initial value of r evaluated is < b
fit <- optim(fn=foo, par=log(sp), method="N", control=list(maxit=5000));
par=exp(fit$par);
par;


########################################

fixD <- fixPar;
f.fd.1par <- function(pars) {
	b <- pars[1];
	r <- b - fixD;
	epsilon <- fixD/b;
	if (r <= 0 | epsilon <= 0 | epsilon >= 1) {return(-Inf);}
	l.int <- n.int * log(r) - r * sum.int.t.len - sum(log(1 - (epsilon * exp(-r * int.t.0))));
	ert <- exp(r * pend.t.len);
	B <- (ert - 1) / (ert - epsilon);
	l.pend <- sum(log(1 - B) + pend.n.t.minus.1*log(B));
	return(l.int + l.pend);
}

lik <- f.fd.1par;
foo <- function (x) {-lik(pars=exp(x));}
fit <- optimize(f=foo, interval=c(log(fixPar),1));
b <- exp(fit$minimum);
par <- c((b - fixPar), (fixPar / b));
par;


########################################


fixD <- fixPar;
f.fd.2par <- function(pars) {
	r <- pars[1] - fixD;
	epsilon <- fixD/pars[2];
	if (r <= 0 | epsilon <= 0 | epsilon >= 1) {return(-Inf);}
	l.int <- n.int * log(r) - r * sum.int.t.len - sum(log(1 - (epsilon * exp(-r * int.t.0))));
	ert <- exp(r * pend.t.len);
	B <- (ert - 1) / (ert - epsilon);
	l.pend <- sum(log(1 - B) + pend.n.t.minus.1*log(B));
	return(l.int + l.pend);
}

lik <- f.fd.2par;
foo <- function (x) {-lik(pars=exp(x));}
sp=c(0.05, 0.5); # make sure initial value of r evaluated is < b
fit <- optim(fn=foo, par=log(sp), method="N", control=list(maxit=5000));
par=c((exp(fit$par[1]) - fixPar),(fixPar/exp(fit$par[2])));
par;

