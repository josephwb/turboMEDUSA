\name{medusaSummary}
\alias{medusaSummary}
\title{
Summarize the fit of a piecewise diversification model
}
\description{
Summarize results from a \code{\link{MEDUSA}} analysis.
}
\usage{
medusaSummary(results, plotTree=TRUE, time=TRUE, node.labels=TRUE, colourOffset=0, 
	cex=0.5, plotSurface=FALSE, printTitle=TRUE, n.points=100, ...)
}

\arguments{
  \item{results}{
An object of class \code{"medusa"} produced from a \code{\link{MEDUSA}} analysis.
}
  \item{plotTree}{
Boolean value that determines whether the fitted tree is plotted. If \code{plotTree=TRUE} (the default) the tree is plotted with edges coloured according to which piecewise model they belong to, and positions (nodes) of break points are potentially indicated (see \code{node.labels} below). 
}
  \item{time}{
Boolean value that indicates whether a timescale is printed below the tree. Only has meaning if \code{plotTree=TRUE}. For large trees it might be desirable to set \code{time=FALSE} to make tip labels legible.
}
  \item{node.labels}{
Boolean value that indicates whether node breaks are plotted on the tree. Only has meaning if \code{plotTree=TRUE}.
}
  \item{colourOffset}{
A plotting option for the colouring of branches by piecewise model membership. An integer value > 1; sets the colour cycling from a different starting point. Only has meaning if \code{plotTree=TRUE}.
}
  \item{cex}{
Specifies the font size for printing tip labels. Default is \code{cex=0.5}; may need to decrease for larger trees. Only has meaning if \code{plotTree=TRUE}.
}
  \item{plotSurface}{
Boolean value that indicates whether likelihood surface(s) for the selected piecewise model is plotted. If \code{plotSurface=TRUE} the optimal model is evaluated at \code{n.points} (see below), and contours are plotted. Default is no plot. 
}
  \item{printTitle}{
Whether to print titles on surface plots. 
}
  \item{n.points}{
The number of values at which r and epsilon are evaluated for the selected piecewise model (hence, a total of \code{n.points} x \code{n.points} parameter combinations considered for model \code{"bd"}). A greater value of \code{n.points} gives better resolution of the surface, but comes with a computational burden. Default is \code{n.points=100}. Only has meaning if \code{plotSurface=TRUE}.
}
  \item{\dots}{
Potential additional arguments for plotting etc.
}
}

\value{
An object of class \code{"medusa.sum"}. This is a list with elements:
\item{$z}{A matrix listing branch times and richnesses, coded by model membership. The model used is either the optimal model (the default) or a user-selected mdoel.}
\item{$mm}{Colour-coding for edges.}
\item{$break.point}{The position in the tree where shifts in models occur.}
\item{$phy}{The tree (possibly pruned from the original tree using richness information).}
This object can be passed to \code{plotPrettyTree} to produce publication-quality trees.
}

\author{
Joseph W. Brown, Richard G. FitzJohn, Michael E. Alfaro, and Luke J. Harmon.

Maintainer: Joseph W. Brown <josephwb@uidaho.edu>}

\seealso{
\code{\link{MEDUSA}}
\code{\link{plotPrettyTree}}
}