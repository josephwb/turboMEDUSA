\name{MEDUSA-package}
\alias{MEDUSA-package}
\docType{package}
\title{
MEDUSA: Modeling Evolutionary Diversification Using Stepwise AIC
}
\description{
Fits piecewise birth-death models to an ultrametric phylogenetic tree according to phylogenetic (edge-length) and taxonomic (richness) likelihoods. A reimplementation of the original \pkg{\link{geiger}} MEDUSA method, meant especially for the analysis of large trees. Optimal model size is determined via a stepwise AIC approach. Run with \code{\link{MEDUSA}}. For a single tree,  summarize results with \code{\link{medusaSummary}}. Publication-quality trees with coloured branches (reflecting piecewise model membership) can be plotted using \code{\link{plotPrettyTree}} For a distrbution of trees, summarize results with \code{\link{multiMedusaSummary}}, and plot fancy figures with \code{\link{plotMultiMedusa}}.
}
\details{
\tabular{ll}{
Package: \tab MEDUSA\cr
Type: \tab Package\cr
Version: \tab 0.93-4-12\cr
Date: \tab 2012-09-19\cr
License: \tab GPL version 2 or greater\cr
LazyLoad: \tab yes\cr
 }
}

\references{
Alfaro, ME, F Santini, C Brock, H Alamillo, A Dornburg, DL Rabosky, G Carnevale, and LJ Harmon. 2009. Nine exceptional radiations plus high turnover explain species diversity in jawed vertebrates. \emph{Proceedings of the National Academy of Sciences} \bold{106}: 13410-13414.
}

\author{
Joseph W. Brown, Richard G. FitzJohn, Michael E. Alfaro, and Luke J. Harmon.

Maintainer: Joseph W. Brown <josephwb@uidaho.edu>}