\name{plotPrettyTree}
\alias{plotPrettyTree}
\title{
Produce a publication-quality tree figure indicating the position of rate shifts
}
\description{
Plot a publication-quality colour-coded tree, where edge colours indicate piecewise diversification model membership from a \code{\link{MEDUSA}} analysis. A wrapper for \code{\link{plot.phylo}} in the APE package.
}
\usage{
plotPrettyTree(treeParameters, time=TRUE, node.labels=FALSE, margin=FALSE, cex=0.5,
	label.offset=0, font=3, colourTipLabels=FALSE, ...)
}

\arguments{
  \item{treeParameters}{
An object of class \code{"medusa.sum"} produced from a \code{\link{medusaSummary}} af a \code{\link{MEDUSA}} analysis. Contains the tree and piecewise diversification model membership.
}
  \item{time}{
Whether a timescale should be plotted. Possibly not wanted for large trees. Default=true;
}
  \item{node.labels}{
Whether node labels should be plotted on the tree to indicate the position and order of diversification rate shifts. 
}
  \item{margin}{
Whether plot should include margins. \code{time} (above) overrules this option, as a margin is required for plotting a timescale.
}
  \item{cex}{
Specifies the font size for printing tip labels. Default is \code{cex=0.5}; may need to decrease for larger trees.
}
  \item{label.offset}{
Specifies how far tip labels are offset rightward from the default \code{plot.phylo()} settings.
}
  \item{font}{
Specifies the font to display for tip labels. The default \code{font=3} plots labels in italic. Options are \code{font=1} (plain text), \code{font=2} (bold), \code{font=3} (italic, the default), and \code{font=4} (bold italic).
}
  \item{colourTipLabels}{
Determines whether tip labels should be coloured according to their piecewise diversification model membership. Default=FALSE.
}
  \item{\dots}{
Potential additional arguments to be passed to \code{\link{plot.phylo}}.
}
}

\value{
None.
}

\author{
Joseph W. Brown, Richard G. FitzJohn, Michael E. Alfaro, and Luke J. Harmon.

Maintainer: Joseph W. Brown <josephwb@uidaho.edu>}

\seealso{
\code{\link{MEDUSA}}
\code{\link{medusaSummary}}
}