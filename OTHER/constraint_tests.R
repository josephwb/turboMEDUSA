f.fixD <- function(pars) {
				b <- pars[1];
				r <- b - fixD;
				epsilon <- fixD/b;
				if (r <= 0 | epsilon < 0 | epsilon >= 1) {return(-Inf);}
				l.int <- n.int * log(r) - r * sum.int.t.len - sum(log(1 - (epsilon * exp(-r * int.t.0))));
				ert <- exp(r * pend.t.len);
				B <- (ert - 1) / (ert - epsilon);
				l.pend <- sum(log(1 - B) + pend.n.t.minus.1*log(B));
				return(l.int + l.pend);
}

foo.fd <- function (x) {-f.fixD(pars=exp(x));}
x <- 1;
maxVal <- fixD + x;
fit.fd <- optimize(f=foo.fd, interval=c(log(fixD), log(maxVal)));
while (fit.fd$objective == Inf)
{
	x <- x/2;
	maxVal <- fixD + x;
	fit.fd <- optimize(f=foo.fd, interval=c(log(fixD), log(maxVal)));
}
fd.b <- exp(fit$minimum);
par.fd <- c((fd.b - fixPar), (fixPar / fd.b));

###########################################################

f.fixB <- function(pars) {
				d <- pars[1];
				r <- fixB - d;
				epsilon <- d/fixB;
#				cat("r =", r, "; epsilon =", epsilon, "\n")
				if (r <= 0 | epsilon < 0 | epsilon >= 1) {return(-Inf);}
				l.int <- n.int * log(r) - r * sum.int.t.len - sum(log(1 - (epsilon * exp(-r * int.t.0))));
				ert <- exp(r * pend.t.len);
				B <- (ert - 1) / (ert - epsilon);
				l.pend <- sum(log(1 - B) + pend.n.t.minus.1*log(B));
				return(l.int + l.pend);
			}

foo.fb <- function (x) {-f.fixB(pars=exp(x));}
fit.fb <- optimize(f=foo.fb, interval=c(-50, log(fixB)));
fb.d <- exp(fit.fb$minimum);
par.fb <- c((fixB - fb.d), (fb.d / fixB));

###########################################################

f.bd <- function(pars) {
				r <- pars[1];
				epsilon <- pars[2];
				if (r <= 0 | epsilon < 0 | epsilon >= 1) {return(-Inf);}
				l.int <- n.int * log(r) - r * sum.int.t.len - sum(log(1 - (epsilon * exp(-r * int.t.0))));
				ert <- exp(r * pend.t.len);
				B <- (ert - 1) / (ert - epsilon);
				l.pend <- sum(log(1 - B) + pend.n.t.minus.1*log(B));
				return(l.int + l.pend);
			}

foo.bd <- function (x) {-f.bd(pars=exp(x));}
sp=c(0.05, 0.5)
fit.bd <- optim(fn=foo.bd, par=log(sp), method="N", control=list(maxit=5000));
par.bd <- exp(c(fit.bd$par));
BD <- getBD(par.bd);