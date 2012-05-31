\name{multiMedusaSummary}
\alias{multiMedusaSummary}
\title{
Summarize the fit of piecewise diversification models across a distribution of trees
}
\description{
Summarize results from \code{\link{MEDUSA}} analysis across a distribution of trees.
}
\usage{
multiMedusaSummary(res, conTree, cutOff=0.05, plotModelSizes=FALSE, cex=0.5, ...)
}

\arguments{
  \item{res}{
An object of class \code{"multiMedusa"} produced from a \code{\link{MEDUSA}} analysis across a distribution of trees.
}
  \item{conTree}{
Tree of class \code{"phylo"} on which to summarize results. Although in principle any tree can be used, ideally this tree should share as many edges as possible with those across the distribution of analyzed trees e.g. some sort of consensus or maximum clade credibility tree. Rate shift positons and average edge-specific rates are mapped onto this tree. May be pruned using richness information (a component of \code{results}) to correspond with tip sampling of trees that were analyzed. 
}
  \item{cutOff}{
Determines the threshold (expressed as a proportion) for which splits are recorded on \code{conTree}. Default is to keep shifts that occur in > 5% of trees.
}
  \item{plotModelSizes}{
Boolean value that determines whether a histogram of the number of piecewise diversification models per tree is plotted. Default is no plot.
}
  \item{cex}{
Specifies the font size for printing tip labels. Default is \code{cex=0.5}; may need to decrease for larger trees. Only has meaning if \code{plotTree=TRUE}.
}
  \item{\dots}{
Potential additional arguments for plotting etc.
}
}

\value{
An object of class \code{"multiMedusaSummary"}. This is a list with elements:
\item{$model.sizes}{A matrix listing branch times and richnesses, coded by model membership. The model used is either the optimal model (the default) or a user-selected mdoel.}
\item{$num.trees}{The number of trees analyzed.}
\item{$shift.summary}{The position in the tree where shifts in models occur, together with their frequecny across trees.}
\item{$summary.tree}{The tree (possibly pruned from the original tree using richness information) on which the results were summarized. Of class \code{"phylo"}. In addition to standard \code{"phylo"} elements, includes a matrix \code{$rates}. This matrix is oredered according to \code{summary.tree$edge}, and includes columns: \code{r.mean}, \code{r.median}, \code{r.sd}, \code{eps.mean}, \code{eps.median}, \code{eps.sd}, and \code{freq}. \code{freq} corresponds to the frequency with which an edge is present accross all trees.}
This object can be passed to \code{plotPrettyTree} to produce publication-quality trees.
}

\author{
Joseph W. Brown, Richard G. FitzJohn, Michael E. Alfaro, and Luke J. Harmon.

Maintainer: Joseph W. Brown <josephwb@uidaho.edu>}

\seealso{
\code{\link{MEDUSA}}
\code{\link{plotPrettyTree}}
}