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
calculateStatistics <- function(dataFp, designFp, comparisons, limmaTest=TRUE, varFilterFrac=1, logTrans=FALSE, robustLimma=FALSE, type="limma") {
    
    nst <- setupStatisticsObject(dataFp, designFp, comparisons, logTrans=logTrans)
    
    if (limmaTest) {
        nst <- calculatePairwiseComparisonsLimma(nst, comparisons, "group", robustLimma=robustLimma)
    }
    else {
        nst <- calculateContrasts(nst, comparisons, "group", type=type)
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

outputAnnotatedMatrix <- function(nst, outPath) {
    
    pairwiseHead <- paste(names(nst@pairwiseCompsP), "p", sep="_")
    pMat <- data.frame(nst@pairwiseCompsP)
    colnames(pMat) <- pairwiseHead
    
    pairwiseHeadFdr <- paste(names(nst@pairwiseCompsFdr), "p.adj", sep="_")
    fdrMat <- data.frame(nst@pairwiseCompsFdr)
    colnames(fdrMat) <- pairwiseHeadFdr
    
    pairwiseHeadAve <- paste(names(nst@pairwiseCompsAve), "ave.expr", sep="_")
    aveMat <- data.frame(nst@pairwiseCompsAve)
    colnames(aveMat) <- pairwiseHeadAve
    
    pairwiseHeadFold <- paste(names(nst@pairwiseCompsFold), "fold", sep="_")
    foldMat <- data.frame(nst@pairwiseCompsFold)
    colnames(foldMat) <- pairwiseHeadFold
    
    outDf <- cbind(nst@annotMat, pMat, fdrMat, aveMat, foldMat, nst@dataMat)
    outDf
}

outputStatsReport <- function(nst, outPath) {
    
}