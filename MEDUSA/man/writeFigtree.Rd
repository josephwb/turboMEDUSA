\name{writeFigtree}
\alias{writeFigtree}
\title{
Write a Nexus-flavoured tree file with rate and shift annotations.
}
\description{
Write a summary tree annotated with results from a \code{\link{MEDUSA}} analysis across a distribution of trees.
}
\usage{
writeFigtree(summary, file="", digits=10)
}

\arguments{
  \item{summary}{
An object of class \code{"multiMedusaSummary"} produced from a \code{\link{MEDUSA}} analysis.
}
  \item{file}{
The file to write the tree to. Default=screen;
}
  \item{digits}{
How many digits to write. Default=10. 
}
}

\value{
Tree string.
}

\author{
Joseph W. Brown, Richard G. FitzJohn, Michael E. Alfaro, and Luke J. Harmon.

Maintainer: Joseph W. Brown <josephwb@uidaho.edu>}

\seealso{
\code{\link{MEDUSA}}
\code{\link{medusaSummary}}
}