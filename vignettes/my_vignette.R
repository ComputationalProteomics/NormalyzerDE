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
fullDf <- read.csv("data.tsv", sep="\t")
designDf <- read.csv("design.tsv", sep="\t")
head(fullDf)
head(designDf)

## ------------------------------------------------------------------------
sampleNames <- as.character(designDf$sample)
dataMat <- as.matrix(fullDf[, sampleNames])
retentionTimes <- fullDf$Average.RT

head(dataMat)

## ------------------------------------------------------------------------
typeof(dataMat)

print("Rows and columns of data")
dim(dataMat)

print("Number of retention times")
length(retentionTimes)

## ------------------------------------------------------------------------
performCyclicLoessNormalization <- function(rawMatrix) {
    log2Matrix <- log2(rawMatrix)
    normMatrix <- limma::normalizeCyclicLoess(log2Matrix, method="fast")
    colnames(normMatrix) <- colnames(rawMatrix)
    normMatrix
}

## ------------------------------------------------------------------------
rtNormMat <- getRTNormalizedMatrix(dataMat, retentionTimes, performCyclicLoessNormalization, stepSizeMinutes=1, windowMinCount=100)

## ------------------------------------------------------------------------
globalNormMat <- performCyclicLoessNormalization(dataMat)
dim(rtNormMat)
dim(globalNormMat)
head(rtNormMat, 1)
head(globalNormMat, 1)

## ------------------------------------------------------------------------
layeredRtNormMat <- getSmoothedRTNormalizedMatrix(dataMat, retentionTimes, performCyclicLoessNormalization, stepSizeMinutes=1, windowMinCount=100, frameShifts=3, mergeMethod="mean")

dim(layeredRtNormMat)
head(layeredRtNormMat, 1)

## ------------------------------------------------------------------------
jobName <- "vignette_run"
normObj <- getVerifiedNormalyzerObject("data.tsv", jobName, "design.tsv")

## ------------------------------------------------------------------------
normResults <- normMethods(normObj)

## ------------------------------------------------------------------------
normResultsWithEval <- analyzeNormalizations(normResults)

## ------------------------------------------------------------------------

## ------------------------------------------------------------------------
jobDir <- setupJobDir("vignette_run", ".")
writeNormalizedDatasets(normResultsWithEval, jobDir)

## ------------------------------------------------------------------------
generatePlots(normResultsWithEval, jobDir)

