## ----setup, echo=FALSE, results="hide"-----------------------------------
# knitr::opts_chunk$set(tidy=FALSE, cache=TRUE, dev="png", message=FALSE,
# error=FALSE, warning=TRUE)

## ------------------------------------------------------------------------
library(NormalyzerDE)
#getwd()
normalyzer
normalyzer("data.tsv", "vignette_run", designMatrix="design.tsv", outputDir="testout")

