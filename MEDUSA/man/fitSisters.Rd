\name{fitSisters}
\alias{fitSisters}
\alias{print.sisterFit}
\title{
Fit piecewise diversification models to sister clades in a time-calibrated phylogeny
}
\description{
Fits individual diversification models to the two sister clades branching from a specified node.
}

\usage{
fitSisters(phy, richness=NULL, node=NULL, tips=NULL, model="mixed",
	criterion="aicc", epsilon=NULL, r=NULL, b=NULL, d=NULL, plotSurface=FALSE,
	initialR=0.05, initialE=0.5, verbose=TRUE, mc=FALSE, numCores=NULL, ...)
}

\arguments{
  \item{phy}{
Ultrametric phylogenetic tree. Assumed to be of class \code{phylo} or (if multiple trees) \code{multiPhylo}.
}
  \item{richness}{
Optional; only required if tree is NOT completely sampled. A dataframe with minimally two columns \code{taxon} and \code{n.taxa}. The identifier \code{taxon} must match exactly with a tip.label in the phylogeny \code{phy}; this is checked. May also include an \code{exemplar} column, used for renaming incompletely-sampled clades that require collapsing; in this case, \code{exemplar} must match exactly with a tip label and \code{taxon} will be the name for that clade (say, genus, family, etc.). If no richness information is provided then it is assumed tips represent single species with complete sampling.
}
  \item{node}{
The node where the two clades of interest share a common ancestor. Value may be \code{node="root"} if the root node, or a numeric value based on APE's node numbering system. If the node number is unknown, \code{tips} can be provided instead (see below).
}
  \item{tips}{
If the node number (using APE's numbering system) of a node of interest is not known, tip labels for taxa sharing a MRCA at that node can instead be passed in. Specified as: \code{tips=c("tax1", "tax2")}, where tip labels must be given in quotes.
}
  \item{model}{
The flavour of piecewise models to be fit to the phylogeny. Can be birth-death (\code{model="bd"}), pure-birth (\code{model="yule"}), or mixed (\code{model="mixed"}) (the default). For \code{mixed} model analyses, both BD and Yule models are considered for all potential breaks.
}
  \item{criterion}{
The statistical criterion used to determine whether larger models fit the data significantly better. Options are \code{criterion="aicc"} (default) or \code{criterion="aic"}. As AIC and AICc will select the same breaks for any given model size, this option is only useful if \code{stop="threshold"}. In this case, the default is strongly advised, as for some trees AIC can continue to improve well after AICc has deteriorated.
}
  \item{epsilon}{
A fixed value for the extinction fraction. Will be held constant while r (net diversification) is estimated. Only for use in the birth-death model.
}
  \item{r}{
A fixed value for the net diversification rate. Will be held constant while epsilon (extinction fraction) is estimated. Only for use in the birth-death model.
}
  \item{b}{
A fixed value for the per-lineage speciation rate. Only for use in the birth-death model.
}
  \item{d}{
A fixed value for the per-lineage extinction rate. Will be held constant while epsilon (extinction fraction) is estimated.
}
  \item{plotSurface}{
Whether parameter likelihood surfaces should be plotted. Default is FALSE.
}
  \item{initialR}{
Starting guess for net diversification rate. May need to be decreased for the program to work on large trees.
}
  \item{initialE}{
Starting guess for relative extinction rate.
}
  \item{verbose}{
Print out additional information to screen, such as which taxa are pruned from the original tree given richness information. Mostly useful for debugging.
}
  \item{mc}{
Indicate whether multiple cores should be used (\code{mc=TRUE}). Works only on unix-like machines (i.e. not Windows), and only in non-GUI instantiations of R.
}
  \item{numCores}{
The number of cores to use if \code{mc=TRUE}. By default all cores detected will be used. This option is only useful if one wishes to limit the number of cares devoted to \code{MEDUSA} analysis.
}
  \item{\dots}{
Potential additional arguments for plotting etc.
}
}
\details{
A diversification model is first fit to the entire subtree rooted at node \code{node}. Separate diversification models are then fit to the two clades emanating from \code{node}. MLEs and associated precision intervals are reported for parameters in both the \code{Base} and \code{Sisters} models. If a single tree in passed in, confidence intervals are calculated using profile likelihoods. If multiple trees are analyzed, the following are reported for each parameter: mean, min, max, and standard deviation.

In addition, the relative fit of the \code{Base} and \code{Sisters} models is also reported using AIC weights. In the case of multiple trees, the mean and precision of model weights is presented, as well as the proportion of times each model is the better fit model.
}

\value{
An object of class \code{"sisterFit"}. Returned invisibly. This is a list with elements:
\item{$parameterSummary}{A dataframe listing parameter estimates and precision statistics.}
\item{$modelSummary}{A dataframe summarizing the relative fit of a single model vs. the implemented 2-piece sister-lineages model.}
\item{$fixPar}{The fixed parameter value used in a constrained model. If a constrained model is not implemented, this will be NULL.}
\item{$models}{A list containing the base and sister-clade diversificaton models estimated for each tree. Includes tree(s) (which may have been pruned using richness information).}
}

\author{
Joseph W. Brown.

Maintainer: Joseph W. Brown <josephwb@uidaho.edu>}