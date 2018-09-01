## ----setup, echo=FALSE, results="hide"-----------------------------------
# knitr::opts_chunk$set(tidy=FALSE, cache=TRUE, dev="png", message=FALSE,
# error=FALSE, warning=TRUE)

## ------------------------------------------------------------------------
library(NormalyzerDE)
outDir <- tempdir()
design_fp <- system.file(package="NormalyzerDE", "extdata", "design.tsv")
data_fp <- system.file(package="NormalyzerDE", "extdata", "data.tsv")
print(design_fp)
print(data_fp)
normalyzer(jobName="vignette_run", designPath=design_fp, dataPath=data_fp, 
           outputDir=outDir)

## ------------------------------------------------------------------------
outDir <- tempdir()
normMatrixPath <- paste(outDir, "vignette_run/Loess-G-normalized.txt", sep="/")
normalyzerDE("vignette_run", 
             design_fp, 
             normMatrixPath,
             outputDir=outDir, 
             comparisons=c("1-2", "1-3"))

## ------------------------------------------------------------------------
fullDf <- read.csv(data_fp, sep="\t")
designDf <- read.csv(design_fp, sep="\t")
head(fullDf, 1)
head(designDf, 1)

## ------------------------------------------------------------------------
sampleNames <- as.character(designDf$sample)
typeof(sampleNames)

## ------------------------------------------------------------------------
dataMat <- as.matrix(fullDf[, sampleNames])
retentionTimes <- fullDf$Average.RT

head(dataMat, 1)

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
rtNormMat <- getRTNormalizedMatrix(dataMat, 
                                   retentionTimes, 
                                   performCyclicLoessNormalization, 
                                   stepSizeMinutes=1, 
                                   windowMinCount=100)

## ------------------------------------------------------------------------
globalNormMat <- performCyclicLoessNormalization(dataMat)
dim(rtNormMat)
dim(globalNormMat)
head(rtNormMat, 1)
head(globalNormMat, 1)

## ------------------------------------------------------------------------
layeredRtNormMat <- getSmoothedRTNormalizedMatrix(
    dataMat, 
    retentionTimes, 
    performCyclicLoessNormalization, 
    stepSizeMinutes=1, 
    windowMinCount=100, 
    windowShifts=3, 
    mergeMethod="mean")

dim(layeredRtNormMat)
head(layeredRtNormMat, 1)

## ------------------------------------------------------------------------
jobName <- "vignette_run"
designDf <- loadDesign(design_fp)
dataDf <- loadData(data_fp)
normObj <- getVerifiedNormalyzerObject(jobName, designDf, dataDf)

## ------------------------------------------------------------------------
normResults <- normMethods(normObj)

## ------------------------------------------------------------------------
normResultsWithEval <- analyzeNormalizations(normResults)

## ------------------------------------------------------------------------
jobDir <- setupJobDir("vignette_run", tempdir())
writeNormalizedDatasets(normResultsWithEval, jobDir)

## ------------------------------------------------------------------------
generatePlots(normResultsWithEval, jobDir)

## ------------------------------------------------------------------------
bestNormMatPath <- paste(jobDir, "RT-Loess-normalized.txt", sep="/")
fullDf <- utils::read.csv(bestNormMatPath, sep="\t")
designDf <- utils::read.csv(design_fp, sep="\t")
nst <- setupStatisticsObject(designDf, fullDf)

## ------------------------------------------------------------------------
comparisons <- c("1-2", "2-3")
nst <- calculateContrasts(nst, comparisons, condCol="group")

## ------------------------------------------------------------------------
annotDf <- generateAnnotatedMatrix(nst)
utils::write.table(annotDf, file=paste(jobDir, "stat_table.tsv", sep="/"))
generateStatsReport(nst, "Vignette stats", jobDir)

