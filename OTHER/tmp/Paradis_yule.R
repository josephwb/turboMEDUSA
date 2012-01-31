yule <- function (phy, use.root.edge = FALSE) 
{
    if (!is.binary.tree(phy)) 
        stop("tree must be dichotomous to fit the Yule model.")
    bt <- rev(sort(branching.times(phy)))
    ni <- cumsum(rev(table(bt))) + 1
    X <- sum(phy$edge.length)
    nb.node <- phy$Nnode
    if (!is.null(phy$root.edge) && use.root.edge) {
        X <- X + phy$root.edge
        ni <- c(1, ni)
    } else nb.node <- nb.node - 1
    lambda <- nb.node/X
    se <- lambda/sqrt(nb.node)
    loglik <- -lambda * X + lfactorial(phy$Nnode) + nb.node * log(lambda)
    obj <- list(lambda = lambda, se = se, loglik = loglik)
    class(obj) <- "yule"
    obj
}

nb.node = num.nodes - 1
X = sum of edge lengths
phy$Nnode = num.nodes
ni = birth time of lineage i

Rabosky (2007): L = N.int * log(r) - r * sum.t
Paradis (2003): L = sum.t * (-r) + (N.int - 1) * log(r) + lfactorial(N.int)

Rabosky
pureBirth <- function (x) 
{
    if (!is.numeric(x)) 
        stop("object x not of class 'numeric'")
    res <- list()
    x <- rev(sort(x))
    temp <- yuleint2(x, x[1], 0)
    res$LH <- temp$LH
    res$aic <- (-2 * res$LH) + 2
    res$r1 <- temp$smax
    return(res)
}
yuleint2  <- function (x, st1, st2) 
{
    nvec <- 2:(length(x) + 1)
    nv <- x[(x < st1) & (x >= st2)]
    lo <- max(nvec[x >= st1])
    up <- max(nvec[x >= st2])
    res <- list()
    if (st1 <= x[1]) { # isn't this always true?
        nv <- c(st1, nv) - st2
    } else {
        nv <- nv - st2
    }
    smax <- (up - lo)/(lo * nv[1] + sum(nv[2:(up - lo + 1)]))
    s1 <- sum(log(lo:(up - 1)))
    s2 <- (up - lo) * log(smax)
    s3 <- -(lo * nv[1] + sum(nv[2:(up - lo + 1)])) * smax
    res$smax <- smax
    res$LH <- s1 + s2 + s3
    return(res)
}