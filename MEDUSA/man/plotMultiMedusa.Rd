\name{plotMultiMedusa}
\alias{plotMultiMedusa}
\title{
Plot the summary tree from a multiMEDUSA analysis.
}
\description{
Plot a consensus tree annotated with results from a \code{\link{MEDUSA}} analysis across a distribution of trees.
}
\usage{
plotMultiMedusa(summary, treeRearrange="down", annotateShift=TRUE, annotateRate="r.median",
	plotRichnesses=TRUE, richPlot="log", time=TRUE, tip.cex=0.3, shiftScale=1, label.offset=0.5, font=3,
	shift.leg.pos="left", power=1.5, ...)
}

\arguments{
  \item{summary}{
An object of class \code{"multiMedusaSummary"} produced from a \code{\link{multiMedusaSummary}} summary across a distribution of trees.
}
  \item{treeRearrange}{
Ladderize the tree for nicer viewing. Default is \code{treeRearrange="up"}. Possible values are \code{"up"}, \code{"down"}, or \code{NULL} (for no rearrangement of edges).
}
  \item{plotRichnesses}{
Boolean value that determines whether bar plot of extant tip richnesses is plotted. Default is \code{true}.
}
  \item{richPlot}{
Flavour of richness to plot. Only has meaning if \code{plotRichnesses=TRUE}. By default, richnesses are plotted as \code{"log"}. Alternative is \code{"raw"} i.e. straight richnesses are plotted.
}
  \item{annotateShift}{
Boolean value that determines whether shift positions are plotted on the tree. Default is \code{true}.
}
  \item{annotateRate}{
Specifies which rate to annotate the tree with. Possibilities are: \code{"r.mean"}, \code{"r.median"}, \code{"r.sd"}, \code{"eps.mean"}, \code{"eps.median"}, and \code{"eps.sd"}.
}
  \item{time}{
Boolean value that indicates whether a timescale is printed below the tree. For large trees it might be desirable to set \code{time=FALSE} to make tip labels legible.
}
  \item{tip.cex}{
Specifies the font size for printing tip labels. Default is \code{cex=0.5}; may need to decrease for larger trees. Only has meaning if \code{plotTree=TRUE}.
}
  \item{shiftScale}{
Alter the size of plotted shift proportions. Takes numeric value.
}
  \item{label.offset}{
Specifies how far tip labels are offset rightward from the default \code{plot.phylo()} settings.
}
  \item{font}{
Specifies the font to display for tip labels. The default \code{font=3} plots labels in italic. Options are \code{font=1} (plain text), \code{font=2} (bold), \code{font=3} (italic, the default), and \code{font=4} (bold italic).
}
  \item{shift.leg.pos}{
Specifies the position of the shift legend. Default is \code{shift.leg.pos="bottomleft"}.
}
  \item{power}{
Alters the colour scheme for annotating rates.
}
  \item{\dots}{
Potential additional arguments for plotting etc.
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
\code{\link{multiMedusaSummary}}
}