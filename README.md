MEDUSA: Modeling Evolutionary Diversification Using Stepwise AIC
===============
### Recent Changes To Note
**19 March 2013:** Replaced apparently corrupted code from last push.

**7 January 2013:** Fixed a bug in multiMedusaSummary when any individual analyzed tree exhibited no significant rate shifts. Those affected may have seen an error along the lines of:

	Error in if (i.cuts[d] == "node") { :  missing value where TRUE/FALSE needed

after which MEDUSA would stubbornly refuse to continue. The function now works, although it is not exceedingly fast.

**14 November 2012:** Running multi-tree analyses require that tip.labels ordering is consistent across all trees. [This is, I believe, ensured when reading in Nexus files possessing a translation table, but not otherwise.]. I have been using ape's ape:::.compressTipLabel() function to do this. However, for large and/or many trees, this can take an enormous amount of time (an analysis might seem to hang at "Managing tip label ordering across trees..."). I've reimplemented the function and (especially when multiple cores are used) it speeds things up considerably. If you are analyzing multiple trees, make sure to use version 0.93-4.17 or higher.
Overview
---------------
Fits piecewise diversification models from time-calibrated phylogenetic tree(s) and complete extant species richness. Will eventually be made available on CRAN (http://cran.r-project.org/) in the Geiger package version 2.0.

### Confusing Name Genealogy
MEDUSA was formerly named "turboMEDUSA" (previously distributed as both scripts and R packages); these are the same method, although the current version is much more general. The previous version of MEDUSA available in the current version Geiger (published in Alfaro et al. 2009. PNAS) is an earlier implementation. Yes, I know this is confusing!

### The Current Method
Besides being faster and more efficient, the current version of MEDUSA does a number things the older versions cannot:

1. MEDUSA can fit both birth-death (BD) and Yule (previously not available) diversification models.
2. "Mixed" models: both Yule and BD models are considered for each possible rate shift position. This is not a true "mixed" model (in which parameter estimates are weighted means across models); rather, it just considers both model flavours and proceeds with the one that has the best AIC score. Using this, you'll likely find your tree painted with Yule and BD models.
3. MEDUSA has always had the strict assumption that shifts occur somewhere within the branch leading to the mrca node of a clade, so rate shifts were always placed at stem nodes. I've generalized this to allow shifts to occur at crown nodes as well. Using this more general shift model often leads to novel inferences of rate shifts, both in their placement and magnitude. For example, if you have a long internal branch with many short terminal branches, the older MEDUSAs (MEDUSAE?!?), by placing the shift at the base of the long branch, would always infer high turnover (the idea being the internal branch is long because extinction erased linages coming off of it, and there exist many short terminal branches because they simply have not gone extinct yet). An alternative interpretation (and one that most people intuit when seeing this pattern) is that of an adaptive radiation near the tips (i.e. at the crown node). The new MEDUSA considers both scenarios.
4. MEDUSA has also always been a strictly parameter-addition algorithm. Oftentimes, a basal clade may be fit with bd; when subsequent shifts are introduced within that clade it is possible that that basal part of the clade is better characterized by yule rather than bd. The new MEDUSA considers such parameter removal, allowing better AIC scores and potentially more/different shifts.
5. More general to #4, MEDUSA now considers removing previously fit rate shifts. As above, when a basal rate shift is fit and subsequent shifts are introduced within that clade, the basal shift may become unnecessary. The parameters regained from removing this shift can again lead to better AIC scores and potentially more/different shifts.
6. The ability to fix parameter values (r, epsilon, b, or d).
7. Calculation of confidence intervals on parameter estimates.
8. The ability to analyze a distribution of trees (e.g. from a BEAST analysis).

Installation
---------------
At the moment, MEDUSA is not available on CRAN, although it will be, as a component of GEIGER 2.0. To install, put the non-decompressed *.tar.gz file (from the 'Stable' directory) in your R working directory. You will need to install some R dependencies first.

In R, type:

	install.packages("colorspace");
	install.packages("shape");

You may need to specify a download location if you have not previously set up a default location.

To install MEDUSA proper, type in R:

<<<<<<< HEAD
	install.packages("MEDUSA_0.93-4-21.tar.gz", repos=NULL, type="source", INSTALL_opts="--byte-compile");

If you instead prefer the command line (Windows prompt, or Mac/Unix Terminal), navigate to the directory containing the MEDUSA *.tar.gz file and type in the prompt:

	R CMD install --byte-compile MEDUSA_0.93-4-21.tar.gz

The "--byte-compile" option is optional and regards performance. In order to use this option, you must have a fortran compiler. Compilers are available on the CRAN homepage (e.g. for Mac: http://cran.r-project.org/bin/macosx/tools/). If you cannot be bothered with this, just omit the option:

	install.packages("MEDUSA_0.93-4-21.tar.gz", repos=NULL, type="source");
=======
	install.packages("MEDUSA_0.93-4-21.tar.gz", repos=NULL, type="source", INSTALL_opts="--byte-compile");

If you instead prefer the command line (Windows prompt, or Mac/Unix Terminal), navigate to the directory containing the MEDUSA *.tar.gz file and type in the prompt:

	R CMD install --byte-compile MEDUSA_0.93-4-21.tar.gz

The "--byte-compile" option is optional and regards performance. In order to use this option, you must have a fortran compiler. Compilers are available on the CRAN homepage (e.g. for Mac: http://cran.r-project.org/bin/macosx/tools/). If you cannot be bothered with this, just omit the option:

	install.packages("MEDUSA_0.93-4-21.tar.gz", repos=NULL, type="source");

If you are having problems, please look at [Issues](https://github.com/josephwb/turboMEDUSA#issues) below. Also, please provide me with detailed information (operating system, R version, APE version, etc.) so that the process can be made painless for future users.

Usage
--------------
To run MEDUSA, you will need:

An ultrametric tree (ideally time-calibrated). Can be a single tree or a distribution. Read in as (depending on tree format):

	phy <- read.tree("treeFileName"); # for newick-formatted tree(s)

or

	phy <- read.nexus("treeFileName"); # for nexus-formatted tree(s)

Richness information. Optional. Only required if your tree(s) has incomplete sampling. If so, tips will represent sampled species + unsampled species. I gather you know this, as this was a focal impetus for the development of MEDUSA in the first place, but let me know if you require further explanation (and see [these slides](https://sites.google.com/site/macroevolutioninr/course-materials/MEDUSA_intro.pdf)). The richness file (I prefer comma-delimited myself) should contain two columns, with the column names of "taxon" and "n.taxa". Read in as:

	richness <- read.csv("richnessFileName"); # for comma-delimited data

or

	richness <- read.table("richnessFileName", header=TRUE, strip.white=TRUE); # for tab-delimited data

To access the MEDUSA man pages, type in R:

	?MEDUSA

At the simplest (i.e. using all default options), with your tree(s) 'phy' and richness 'richness', just type:

	res <- MEDUSA(phy, richness);

If you have a particularly difficult problem (i.e. large tree, or large number of trees), you will benefit from using multiple processing cores, using the R package "multicore". Unfortunately, however, multicore has some strong restrictions. In addition to requiring a POSIX-compliant OS (i.e NOT Windows), it also won't work if you are using graphical devices (e.g. some form of R console). So multi-core processing can only be done in a strictly command-line environment (e.g. Terminal in Mac). There is negligible overhead involved using multiple cores, so scaling is (as far as I can see it) linear with the number of processors. If compatible to your system, install multicore as:

	install.packages("multicore");

and run MEDUSA as:

	res <- MEDUSA(phy, richness, mc=TRUE);

This will work whether you have one or several trees. 'res' contains all of the results (some components of which are invisible, as they are not of direct interest to a user). You will want to save this.

To summarize results for a single tree analysis, type:

	summ <- medusaSummary(res);

To summarize results across a distribution of trees, you will need a single optimal/consensus tree (say, "conTree") to summarize results upon. Read in as above. Now, summarize:

	summ <- multiMedusaSummary(res, conTree);

To save stuff, type:

	save(phy, richness, res, summ, file="My_MEDUSA_results.RData");


To learn more, there are some slides here:

https://sites.google.com/site/macroevolutioninr/course-materials/MEDUSA_intro.pdf

Issues
--------------
MEDUSA installation seems to have some problems with some flavour(s?) of Windows. While I have not yet found the source of the problem, the following seems to help. First, download the most recent version of MEDUSA from the 'Stable' directory. Double-click on the *.tar.gz file you downloaded. This should expand to a directory called simply 'MEDUSA'; this contains all of the code in uncompressed form. Using your Windows command prompt, navigate to the directory containing 'MEDUSA', and type the following:

	R CMD build MEDUSA

This should remake the *.tar.gz file "MEDUSA_0.93-4-21.tar.gz" (or whatever version you have). Now, type:

	R CMD check MEDUSA_0.93-4-21.tar.gz

This will give semi-detailed error reports (if present). The hope is that something on your system doesn't like how my system compresses things; from previous cases, doing it yourself may ensure that everything will work. Now, assuming no errors (or, even if there are, we'll plod ahead anyway), type:

	R CMD install --byte-compile MEDUSA_0.93-4-21.tar.gz

Does it work? If not, please send me the output of the R CMD check MEDUSA_0.93-4-21.tar.gz command.
=======
This should remake the *.tar.gz file "MEDUSA_0.93-4-21.tar.gz" (or whatever version you have). Now, type:

	R CMD check MEDUSA_0.93-4-21.tar.gz

This will give semi-detailed error reports (if present). The hope is that something on your system doesn't like how my system compresses things; from previous cases, doing it yourself may ensure that everything will work. Now, assuming no errors (or, even if there are, we'll plod ahead anyway), type:

	R CMD install --byte-compile MEDUSA_0.93-4-21.tar.gz

Does it work? If not, please send me the output of the R CMD check MEDUSA_0.93-4-21.tar.gz command.
