pkgname <- "turboMEDUSA"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('turboMEDUSA')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("carnivora")
### * carnivora

flush(stderr()); flush(stdout())

### Name: carnivora
### Title: Carnivora Phylogeny
### Aliases: carnivora
### Keywords: datasets

### ** Examples

data(carnivora)
str(carnivora)



cleanEx()
nameEx("whales")
### * whales

flush(stderr()); flush(stdout())

### Name: whales
### Title: Whale Phylogeny
### Aliases: whales
### Keywords: datasets

### ** Examples

data(whales)
str(whales)



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
