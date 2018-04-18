## ----setup, echo=FALSE, results="hide"-----------------------------------
# knitr::opts_chunk$set(tidy=FALSE, cache=TRUE, dev="png", message=FALSE,
# error=FALSE, warning=TRUE)

## ------------------------------------------------------------------------
head(readLines("design.tsv"))

## ------------------------------------------------------------------------
head(readLines("data.tsv"))

## ------------------------------------------------------------------------
library(NormalyzerDE)
normalyzer("data.tsv", "vignette_run", designMatrix="design.tsv", outputDir="testout")

## ------------------------------------------------------------------------
full_df <- read.csv("data.tsv", sep="\t")
design_df <- read.csv("design.tsv", sep="\t")
head(full_df)
head(design_df)

## ------------------------------------------------------------------------
sample_names <- as.character(design_df$sample)
data_m <- as.matrix(full_df[, sample_names])
retention_times <- full_df$Average.RT

head(data_m)

## ------------------------------------------------------------------------
typeof(data_m)

print("Rows and columns of data")
dim(data_m)

print("Number of retention times")
length(retention_times)

## ------------------------------------------------------------------------
performCyclicLoessNormalization <- function(rawMatrix) {
    log2Matrix <- log2(rawMatrix)
    normMatrix <- limma::normalizeCyclicLoess(log2Matrix, method="fast")
    colnames(normMatrix) <- colnames(rawMatrix)
    normMatrix
}

## ------------------------------------------------------------------------
rt_norm_m <- getRTNormalizedMatrix(data_m, retention_times, performCyclicLoessNormalization, stepSizeMinutes=1, windowMinCount=100)

## ------------------------------------------------------------------------
global_norm_m <- performCyclicLoessNormalization(data_m)
dim(rt_norm_m)
dim(global_norm_m)
head(rt_norm_m, 1)
head(global_norm_m, 1)

## ------------------------------------------------------------------------
layered_rt_norm_m <- getSmoothedRTNormalizedMatrix(data_m, retention_times, performCyclicLoessNormalization, stepSizeMinutes=1, windowMinCount=100, frameShifts=3, merge_method="mean")

dim(layered_rt_norm_m)
head(layered_rt_norm_m, 1)

## ------------------------------------------------------------------------
normObj <- getVerifiedNormalyzerObject(...)
jobDir <- setupJobDir(...)

## ------------------------------------------------------------------------
normObj <- getVerifiedNormalyzerObject(...)

## ------------------------------------------------------------------------
normalyzerResultsObject <- normMethods(...)

## ------------------------------------------------------------------------

## ------------------------------------------------------------------------
writeNormalizedDatasets(...)

## ------------------------------------------------------------------------
generatePlots(...)

