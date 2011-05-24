\name{turboMEDUSA}
\alias{turboMEDUSA}
\docType{package}
\title{
MEDUSA: Modeling Evolutionary Diversification Using Stepwise AIC
}
\description{
Fits piecewise birth-death models to an ultrametric phylogenetic tree according to phylogenetic (edge-length) and taxonomic (richness) likelihoods. A reimplementation of the original \pkg{\link{geiger}} MEDUSA method, meant especially for the analysis of large trees. Optimal model size is determined via a stepwise AIC approach. Run with \code{\link{runTurboMEDUSA}} and summarize results with \code{\link{summarize.turboMEDUSA}}.
}
\details{
\tabular{ll}{
Package: \tab turboMEDUSA\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2011-05-19\cr
License: \tab GPL version 2 or greater\cr
LazyLoad: \tab yes\cr
 }
}

\author{
Luke J. Harmon, Daniel L. Rabosky, Richard G. FitzJohn, and Joseph W. Brown

Maintainer: Joseph W. Brown <josephwb@uidaho.edu>
}

\references{
Alfaro, ME, F Santini, C Brock, H Alamillo, A Dornburg, DL Rabosky, G Carnevale, and LJ Harmon. 2009. Nine exceptional radiations plus high turnover explain species diversity in jawed vertebrates. \emph{Proceedings of the National Academy of Sciences} \bold{106}: 13410-13414.
}