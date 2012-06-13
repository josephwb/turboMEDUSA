pkgname <- "MEDUSA"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('MEDUSA')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("whales")
### * whales

flush(stderr()); flush(stdout())

### Name: whales
### Title: Whale phylogeny
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
