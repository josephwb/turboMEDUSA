\name{plotPrettyTree}
\alias{plotPrettyTree}
\title{
MEDUSA: Modeling Evolutionary Diversification Using Stepwise AIC
}
\description{
Plot a pretty colour-coded tree, where edge colours indicate piecewise model membership from an\code{\link{runTurboMEDUSA}} analysis.
}
\usage{
plotPrettyTree(treeParameters, time=TRUE, node.labels=FALSE, cex=0.5, ...)
}

\arguments{
  \item{treeParameters}{
The list produced from a \code{\link{summarizeTurboMEDUSA}} summary.
}
  \item{time}{
Whether a timescale should be plotted. Possibly not wanted for large trees. Default = true;
}
  \item{node.labels}{
Whether node labels should be plotted on the tree to indicate the position of shifts. 
}
  \item{cex}{
Specifies the font size for printing tip labels. Default is \code{cex=0.5}; may need to decrease for larger trees. Only has meaning if \code{plotTree=TRUE}.
}
  \item{\dots}{
Potential additional arguments for plotting etc.
}
}

\value{
None.
}

\references{
Alfaro, ME, F Santini, C Brock, H Alamillo, A Dornburg, DL Rabosky, G Carnevale, and LJ Harmon. 2009. Nine exceptional radiations plus high turnover explain species diversity in jawed vertebrates. \emph{Proceedings of the National Academy of Sciences} \bold{106}: 13410-13414.
}

\author{
Luke J. Harmon, Daniel L. Rabosky, Richard G. FitzJohn, and Joseph W. Brown

Maintainer: Joseph W. Brown <josephwb@uidaho.edu>}

\seealso{
\code{\link{runTurboMEDUSA}}
\code{\link{summarizeTurboMEDUSA}}
}