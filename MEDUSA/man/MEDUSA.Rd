\name{MEDUSA}
\alias{MEDUSA}
\title{
Fit piecewise diversification models to a time-calibrated phylogeny
}
\description{
Fits piecewise birth-death models to an ultrametric phylogenetic tree according to phylogenetic (edge-length) and taxonomic (richness) likelihoods. Optimal model size is determined via a stepwise AIC approach.
}

\usage{
MEDUSA(phy, richness=NULL, model="mixed", modelLimit=20, stop="threshold",
	shiftCut="both", criterion="aicc", stepBack=TRUE, preserveModelFlavour=FALSE, epsilon=NULL,
	r=NULL, b=NULL, d=NULL, fixThreshold=NULL, initialR=0.05, initialE=0.5,
	verbose=TRUE, mc=FALSE, numCores=NULL, ...)
}

\arguments{
  \item{phy}{
Ultrametric phylogenetic tree. Assumed to be of class \code{phylo} or (if multiple trees) \code{multiPhylo}.
}
  \item{richness}{
Optional; only required if tree is NOT completely sampled. A dataframe with minimally two columns \code{taxon} and \code{n.taxa}. The identifier \code{taxon} must match exactly with a tip.label in the phylogeny \code{phy}; this is checked. May also include an \code{exemplar} column, used for renaming incompletely-sampled clades that require collapsing; in this case, \code{exemplar} must match exactly with a tip label and \code{taxon} will be the name for that clade (say, genus, family, etc.). If no richness information is provided then it is assumed tips represent single species with complete sampling.
}
  \item{model}{
The flavour of piecewise models to be fit to the phylogeny. Can be birth-death (\code{model="bd"}), pure-birth (\code{model="yule"}), or mixed (\code{model="mixed"}) (the default). For \code{mixed} model analyses, both BD and Yule models are considered for all potential breaks.
}
  \item{modelLimit}{
The number of piecewise models to be fit to the phylogeny. The program may overrule this value and fit fewer models is the tree is small (such that the denominator in the AIC correction factor becomes undefined).
}
  \item{stop}{
The criterion for when the analysis should be terminated. \code{stop="modelLimit"} will terminate after \code{modelLimit} piecewise models have been fit to the the phylogeny (see above). The default criterion is \code{stop="threshold"}. Here, \code{threshold} refers to the improvement in AIC score that should be considered significant given the number of tips in the phylogeny; analysis will terminate when the addition of a piecewise BD model does not improve AIC score by \code{threshold}.
}
  \item{shiftCut}{
Determines where rate shifts are placed on the tree. Three options are available. For the first option \code{shiftCut="stem"}, when a subtree is selected for submodel consideration, it is assumed that the shift in diversification occurred in the branch leading to the clade. The second option \code{shiftCut="node"} assumes the shift occurred at the node representing the most recent common ancestor of the clade (that is, omitting the leading branch). The final and default option \code{shiftCut="both"} considers both possibilities for shift placement, and carries the shift with the higher AICc score forward to piecewise model selection.
}
  \item{criterion}{
The statistical criterion used to determine whether larger models fit the data significantly better. Options are \code{criterion="aicc"} (default) or \code{criterion="aic"}. As AIC and AICc will select the same breaks for any given model size, this option is only useful if \code{stop="threshold"}. In this case, the default is strongly advised, as for some trees AIC can continue to improve well after AICc has deteriorated.
}
  \item{stepBack}{
Determines whether parameter removal should be considered. Default is TRUE.
}
  \item{preserveModelFlavour}{
Only for use when using \code{model="mixed"}. When a clade has been fit with model i (either \code{"yule"} or \code{"bd"}), and a rate shift is introduced within the clade, this option determines whether the basal rate is reestimated under model i (\code{preserveModelFlavour=TRUE}) or alternative models flavours are considered (\code{preserveModelFlavour=FALSE}). The default is to not preserve model flavour.
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
  \item{fixThreshold}{
User-provided value for what improvement in AIC score is deemed significant. A value provided for \code{fixThreshold} will overrule the default value determined through simulations, whether the provided value is more liberal or conservative than the default. This value is used to determine when to terminate the analysis when \code{stop="threshold"} (see above).
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
A series of diversification models is fit to a combination of phylogenetic and taxonomic data. The input is a phylogenetic tree with branch lengths proportional to time showing the relationship amongst clades, and the diversity for living species in all of those clades. All taxa missing from the tree must be assigned to one of the tip clades in the richness matrix.

The algorithm first fits a single diversification model to the entire dataset. Then, a series of breaks are added, so that different parts of the tree evolve with different parameter values (per-lineage net diversification and relative extinction rates). The program first compares all single-breakpoint models with the overall model, and selects the best one. Then all possible two-breakpoint models are compared with the best single-breakpoint model, and so on.

If \code{phy} is of class \code{multiPhylo} then the analysis will be carried out on each tree in succession.

Summarize results with \code{\link{medusaSummary}}.
}

\value{
An object of class \code{"medusa"}. Returned invisibly. This is a list with elements:
\item{$desc}{A list of all descendants for each node in the tree; for summarizing results only.}
\item{$phy}{The phylogeny that was analyzed. May have been pruned (using richness information) from originally tree passed in. Mainly for summarizing results.}
\item{$modelSummary}{A dataframe containing break and fit values for optimal models at each size considered.}
\item{$fixPar}{The fixed parameter value used in a constrained model. If a constrained model is not implemented, this will be NULL.}
\item{$threshold}{The correct AIC threshold given the size of the tree analyzed. Used in summarizing results and (possibly) terminating the analysis.}
\item{$models}{A list containing the optimal models at each successive size. For each model contains:
  \itemize{
  \item{\code{$z}}{ A matrix listing branch times and richnesses; for summarizing results only.}
  \item{\code{$par}}{ A matrix containing estimates of r and epsilon for all piecewise models.}
  \item{\code{$lnLik.part}}{ The likelihoods of each piecewise model.}
  \item{\code{$lnLik}}{ The total likelihood for the piecewise model.}
  \item{\code{$split.at}}{ The node(s) where individual BD models originate.}
  \item{\code{$aic}}{ AIC for the piecewise model.}
  \item{\code{$aicc}}{ AICc for the piecewise model.}
  \item{\code{$num.par}}{ The number of parameters (total r + epsilon + breakpoints) in the piecewise model.}}
  }
If \code{phy} is of class \code{multiPhylo} then the returned object is of class \code{"multiMedusa"}, and is a list of class \code{"medusa"} objects.
}

\references{
Alfaro, ME, F Santini, C Brock, H Alamillo, A Dornburg, DL Rabosky, G Carnevale, and LJ Harmon. 2009. Nine exceptional radiations plus high turnover explain species diversity in jawed vertebrates. \emph{Proceedings of the National Academy of Sciences} \bold{106}: 13410-13414.
}

\author{
Joseph W. Brown, Richard G. FitzJohn, Michael E. Alfaro, and Luke J. Harmon.

Maintainer: Joseph W. Brown <josephwb@uidaho.edu>}

\seealso{
\code{\link{medusaSummary}}
\code{\link{plotPrettyTree}}
}