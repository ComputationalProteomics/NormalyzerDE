#' Calculate statistics measures
#' 
#' @param nr Normalyzer results object with calculated results.
#' @param comparisons Target sample contrasts to run.
#' @param varFilterFrac Perform variance filtering before tests.
#'
#' @return Normalyzer results with attached statistics object.
#' @export
#' @examples
#' normObj <- getVerifiedNormalyzerObject("data.tsv", "job_name", "design.tsv")
#' normResults <- normMethods(normObj)
#' normStats <- calculateStatistics(normResults)
calculateStatistics <- function(dataFp, designFp, comparisons, limmaTest=TRUE, varFilterFrac=1, logTrans=FALSE, robustLimma=FALSE) {
    
    nst <- setupStatisticsObject(dataFp, designFp, comparisons, logTrans=logTrans)
    
    if (limmaTest) {
        nst <- calculatePairwiseComparisonsLimma(nst, comparisons, "group", robustLimma=robustLimma)
    }
    else {
        nst <- calculatePairwiseComparisons(nst, comparisons, "group")
    }
    
    nst
}

setupStatisticsObject <- function(dataFp, designFp, comparisons, sampleCol="sample", 
                                  conditionCol="group", batchCol=NULL, logTrans=logTrans) {

    fullDf <- read.csv(dataFp, sep="\t")
    designDf <- read.csv(designFp, sep="\t")
    designDf[, sampleCol] <- as.character(designDf[, sampleCol])
    
    dataCols <- designDf[, sampleCol]
    annotMat <- as.matrix(fullDf[, which(!colnames(fullDf) %in% dataCols)])
    dataMat <- as.matrix(fullDf[, which(colnames(fullDf) %in% dataCols)], )
    
    if (logTrans) {
        dataMat <- log2(dataMat)
    }
    
    # nst <- NormalyzerStatistics(dummyStr="test")
    nst <- NormalyzerStatistics(annotMat=annotMat, dataMat=dataMat, designDf=designDf)
    nst
}
