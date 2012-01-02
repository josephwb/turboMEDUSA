## Takes in tree summary data from summarizeTurboMEDUSA
## treeParameters <- list(z, edge.colour, break.pts, phy, labels)
## make it possible to change default colors
plotPrettyTree <- function (treeParameters, time=TRUE, node.labels=FALSE, margin=FALSE, cex=0.5, label.offset=0,
	font=3, colourTipLabels=FALSE, ...)
{
	edge.colour <- treeParameters$edge.colour;
	labels <- treeParameters$labels;
	break.pts <- treeParameters$break.pts;
	phy <- treeParameters$phy;
	z <- treeParameters$z;
	colour <- NULL;
	
	dev.new();
	
	if (colourTipLabels)
	{
		for (i in 1:length(phy$tip.label))
		{
			colour[i] <- as.integer(z[which(z[,"dec"] == i),"partition"]);
		}
	}
	if (time) {margin=TRUE;}
	plot.phylo(phy, edge.color=z[edge.colour,"partition"], no.margin=!margin, cex=cex, label.offset=label.offset, tip.color=colour, font=font, ...);
	if (time)
	{
		axisPhylo(cex.axis=0.75);
		mtext("Divergence Time (MYA)", at=(max(get("last_plot.phylo", envir = .PlotPhyloEnv)$xx)*0.5), side = 1, line = 2, cex=0.75);
	}
	if (node.labels)
	{
		for (i in  1:length(break.pts))
		{
			nodelabels(i, node= break.pts[i], frame = "c", font = 1, cex=0.5);
		}
	}
}
