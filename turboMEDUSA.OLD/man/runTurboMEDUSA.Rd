\name{runTurboMEDUSA}
\alias{runTurboMEDUSA}
\title{
MEDUSA: Modeling Evolutionary Diversification Using Stepwise AIC
}
\description{
Fits piecewise birth-death models to an ultrametric phylogenetic tree according to phylogenetic (edge-length) and taxonomic (richness) likelihoods. Optimal model size is determined via a stepwise AIC approach.
}

\usage{
runTurboMEDUSA(phy, richness = NULL, model.limit = 20, stop = "model.limit", model="bd",
	criterion = "aicc", cutAtStem = TRUE, initial.r = 0.05, initial.e = 0.5, plotFig = FALSE,
	nexus = FALSE, verbose = TRUE, mc = FALSE, num.cores = NULL, ...)
}

\arguments{
  \item{phy}{
Ultrametric phylogenetic tree. Assumed to be of class phylo (i.e. ape format). Nexus trees can be passed in using \code{nexus=TRUE} (see below).
}
  \item{richness}{
Optional; only required if tree is NOT completely sampled. A dataframe with minimally two columns \code{taxon} and \code{n.taxa}. The identifier \code{taxon} must match exactly with a tip.label in the phylogeny \code{phy}; this is checked. May also include an \code{exemplar} column, used for renaming incompletely-sampled clades that require collapsing; in this case, \code{exemplar} must match exactly with a tip label and \code{taxon} will be the name for that clade (say, genus, family, etc.). If no richness information is provided then it is assumed tips represent single species with complete sampling.
}
  \item{model.limit}{
The number of piecewise models to be fit to the phylogeny. The program may overrule this value and fit fewer models is the tree is small (such that the denominator in the AIC correction factor becomes undefined).
}
  \item{model}{
The flavour of piecewise models to be fit to the phylogeny. Can be birth-death (\code{model="bd"}), the default, pure-birth (\code{model="yule"}), or mixed (\code{model="mixed"}). For \code{mixed} model analyses, both BD and Yule models are considered for all potential breaks.
}
  \item{stop}{
The criterion for when the analysis should be terminated. \code{stop="model.limit"} will terminate after \code{model.limit} piecewise models have been fit to the the phylogeny (see above). The alternative criterion is \code{stop="threshold"}. Here \code{threshold} refers to the increase in AIC score that should be considered significant given the number of tips in the phylogeny; analysis will terminate when the addition of a piecewise BD model does not improve AIC score by \code{threshold}.
}
  \item{criterion}{
The statistical criterion used to determine whether larger models fit the data significantly better. Options are \code{criterion="aicc"} (default) or \code{criterion="aic"}. As AIC and AICc will select the same breaks for any given model size, this option is only useful if \code{stop="threshold"}. In this case, the default is strongly advised, as for some trees AIC can continue to improve well after AICc has deteriorated.
}
  \item{cutAtStem}{
Determines where rate shifts are placed on the tree. Three options are available. For the default option \code{cutAtStem=TRUE}, when a subtree is selected for submodel consideration, it is assumed that the shift in diversification occurred in the branch leading to the clade. The second option \code{cutAtStem=FALSE} assumes the shift occurred at the node representing the most recent common ancestor of the clade (that is, omitting the leading branch). The final option \code{cutAtStem="both"} considers both possibilities for shift placement, and carries the shift with the higher AICc score forward to piecewise model selection.
}
  \item{initial.r}{
Starting guess for net diversification rate. May need to be decreased for the program to work on large trees.
}
  \item{initial.e}{
Starting guess for relative extinction rate.
}
  \item{plotFig}{
Boolean value that indicates whether a plot of AIC and AICc against model size should be generated. Default is for no plot.
}
  \item{nexus}{
Boolean value that indicates whether the tree being passed in is in Nexus format. Default \code{nexus=FALSE} (i.e. the tree being passed in is in phylo format)
}
  \item{verbose}{
Print out additional information to screen, such as which taxa are pruned from the original tree given richness information. Mostly useful for debugging.
}
  \item{mc}{
Indicate whether multiple cores should be used (\code{mc=TRUE}). Works only on unix-like machines (i.e. not Windows), and only in non-GUI instantiations of R.
}
  \item{num.cores}{
The number of cores to use if \code{mc=TRUE}. By default all cores detected will be used. This option is only useful if one wants to limits the number of cares devoted to a particular R job.
}
  \item{\dots}{
Potential additional arguments for plotting etc.
}
}
\details{
A series of diversification models is fit to a combination of phylogenetic and taxonomic data. The input is a phylogenetic tree with branch lengths proportional to time showing the relationship amongst clades, and the diversity for living species in all of those clades. All taxa missing from the tree have to be assigned to one of the tip clades in the richness matrix.

The algorithm first fits a single diversification model to the entire dataset. Then, a series of breaks are added, so that different parts of the tree evolve with different parameter values (per-lineage net diversification and relative extinction rates). The program first compares all single-breakpoint models with the overall model, and selects the best one. Then all possible two-breakpoint models are compared with the best single-breakpoint model, and so on.

Summarize results with \code{\link{summarizeTurboMEDUSA}}.
}

\value{
Returned list has elements:
\item{$z}{A matrix listing branch times and richnesses; for summarizing results only.}
\item{$anc}{A list of all ancestors; for summarizing results only.}
\item{$phy}{The phylogeny that was analyzed. May have been pruned (from richness information) from originally tree passed in. Mainly for summarizing results.}
\item{$modelSummary}{A dataframe containing break and fit values for optimal models at each size considered.}
\item{$threshold}{The correct AIC threshold given the size of the tree analyzed. Used in summarizing results and (possibly) terminating the analysis.}
\item{$models}{A list containing the optimal models at each successive size. For each model contains:
  \itemize{
  \item{\code{$par}}{ A matrix containing estimates of r and epsilon for all piecewise models.}
  \item{\code{$lnLik.part}}{ The likelihoods of each piecewise model.}
  \item{\code{$lnLik}}{ The total likelihood for the piecewise model.}
  \item{\code{$split.at}}{ The node(s) where individual BD models originate.}
  \item{\code{$aic}}{ AIC for the piecewise model.}
  \item{\code{$aicc}}{ AICc for the piecewise model.}
  \item{\code{$num.par}}{ The number of parameters (total r + epsilon + breakpoints) in the piecewise model.}}
  }
}

\references{
Alfaro, ME, F Santini, C Brock, H Alamillo, A Dornburg, DL Rabosky, G Carnevale, and LJ Harmon. 2009. Nine exceptional radiations plus high turnover explain species diversity in jawed vertebrates. \emph{Proceedings of the National Academy of Sciences} \bold{106}: 13410-13414.
}

\author{
Luke J. Harmon, Daniel L. Rabosky, Richard G. FitzJohn, and Joseph W. Brown

Maintainer: Joseph W. Brown <josephwb@uidaho.edu>}

\seealso{
\code{\link{summarizeTurboMEDUSA}}
\code{\link{plotPrettyTree}}
}