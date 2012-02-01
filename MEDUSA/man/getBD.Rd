\name{getBD}
\alias{getBD}
\title{Get per-lineage birth and death rates}
\description{Get per-lineage birth (speciation) and death (extinction) rates from inferred net diversification and relative extinction rates.}
\usage{
getBD(r, epsilon)
}
\arguments{
  \item{r}{Net diversification rate (b-d)}
  \item{epsilon}{Relative extinction rate (d/b)}
}
\details{Useful if one is in interested in b (lambda) & d (mu) rather than r and epsilon.}

\value{
Returns a list containing:
\item{$b }{Per-lineage birth (speciation) rate}
\item{$d }{Per-lineage death (extinction) rate}
}