require(turboMEDUSA)


require(geiger)

model.limit=20; stop="model.limit"; model="mixed";
	criterion="aicc"; shiftCut="both"; initialR=0.05; initialE=0.5; plotFig=FALSE; nexus=FALSE;
	verbose=TRUE; mc=FALSE; num.cores=NULL

load("/Users/josephwb/Projects/R_working/turboMEDUSA/turboMEDUSA_BAD/data/whales.RData")
phy <- whales$phy
richness <- whales$richness

load("/Users/josephwb/Projects/R_working/turboMEDUSA/Mammal_results.RData")

runTurboMEDUSA(phy, richness) -> res.normal
 
runTurboMEDUSA(phy, richness, shiftCut="node") -> res.node

runTurboMEDUSA(tmp$phy, tmp$richness)

phy <- mam.res$phy