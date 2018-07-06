# #' Calculate statistics measures
# #' 
# #' @param nr Normalyzer results object with calculated results.
# #' @param comparisons Target sample contrasts to run.
# #' @param varFilterFrac Perform variance filtering before tests.
# #'
# #' @return Normalyzer results with attached statistics object.
# #' @export
# #' @examples
# #' normObj <- getVerifiedNormalyzerObject("job_name", "design.tsv", "data.tsv")
# #' normResults <- normMethods(normObj)
# #' normStats <- calculateStatistics(normResults)
# calculateStatistics <- function(dataFp, designFp, comparisons, limmaTest=TRUE, varFilterFrac=1, 
#                                 logTrans=FALSE, robustLimma=FALSE, type="limma", batchCol=NULL) {
#     
#     nst <- setupStatisticsObject(dataFp, designFp, comparisons, logTrans=logTrans)
#     nst <- calculateContrasts(nst, comparisons, condCol="group", type=type, batchCol=batchCol)
# 
#     nst
# }

filterLowRep <- function(df, groups, leastRep = 2) {
    rowMeetThresContrast <- apply(df, 1, allReplicatesHaveValuesContrast, 
                                     groups = groups, minCount = leastRep)
    filteredDf <- df[rowMeetThresContrast, ]
    filteredDf
}

allReplicatesHaveValuesContrast <- function(row, groups, minCount) {
    names(row) <- groups
    repCounts <- table(names(stats::na.omit(row)))
    length(repCounts) == length(unique(groups)) && min(repCounts) >= minCount
}

setupStatisticsObject <- function(designFp, dataFp, comparisons, sampleCol="sample", 
                                  conditionCol="group", batchCol=NULL, logTrans=logTrans,
                                  leastRepCount=2) {

    fullDf <- utils::read.csv(dataFp, sep="\t")
    designDf <- utils::read.csv(designFp, sep="\t")
    designDf[, sampleCol] <- as.character(designDf[, sampleCol])

    dataCols <- designDf[, sampleCol]
    annotMat <- as.matrix(fullDf[, which(!colnames(fullDf) %in% dataCols)])
    dataMat <- as.matrix(fullDf[, dataCols], )
    
    rownames(dataMat) <- 1:nrow(dataMat)
    dataMatNAFiltered <- filterLowRep(dataMat, designDf[, conditionCol], leastRep=leastRepCount)
    naFilterContrast <- rownames(dataMat) %in% rownames(dataMatNAFiltered)
      
    if (logTrans) {
        dataMat <- log2(dataMat)
    }
    
    nst <- NormalyzerStatistics(
        annotMat=annotMat, 
        dataMat=dataMat, 
        designDf=designDf, 
        contrasts=comparisons,
        filteredDataMat=dataMatNAFiltered,
        filteringContrast=naFilterContrast
    )
    
    nst
}

generateAnnotatedMatrix <- function(nst) {
    
    pairwiseHead <- paste(names(nst@pairwiseCompsP), "PValue", sep="_")
    pMat <- data.frame(nst@pairwiseCompsP)
    colnames(pMat) <- pairwiseHead
    
    pairwiseHeadFdr <- paste(names(nst@pairwiseCompsFdr), "AdjPVal", sep="_")
    fdrMat <- data.frame(nst@pairwiseCompsFdr)
    colnames(fdrMat) <- pairwiseHeadFdr
    
    pairwiseHeadFold <- paste(names(nst@pairwiseCompsFold), "log2FoldChange", sep="_")
    foldMat <- data.frame(nst@pairwiseCompsFold)
    colnames(foldMat) <- pairwiseHeadFold
    
    aveMat <- data.frame(nst@pairwiseCompsAve)
    
    outDf <- cbind(nst@annotMat, pMat, fdrMat, foldMat, featureAvg=aveMat[, 1], nst@dataMat)
    outDf
}

#' Generate full output report plot document
#' 
#' @param nst NormalyzerDE statistics object.
#' @param jobName Name of processing run.
#' @param jobDir Path to output directory.
#' @param plotRows Number of plot rows.
#' @param plotCols Number of plot columns.
#' @return None
#' @export
#' @examples
#' normObj <- getVerifiedNormalyzerObject("job_name", "design.tsv", "data.tsv")
#' normResults <- normMethods(normObj)
#' normResultsWithEval <- analyzeNormalizations(normResults)
#' generatePlots(normResultsWithEval, "path/to/output")
generateStatsReport <- function(nst, jobName, jobDir, plotRows=3, plotCols=4) {
    
    nrows <- plotRows + 2
    ncols <- plotCols + 2
    
    currentLayout <- grid::grid.layout(nrow=nrows, ncol=ncols,
                                       heights=c(0.1, rep(3 / (nrows-2), (nrows - 2)), 0.1), 
                                       widths=c(0.1, rep(4 / (ncols-2), (ncols - 2)), 0.1), 
                                       default.units=c("null", "null"))
    
    currentFont <- "Helvetica"
    setupPlotting(jobName, jobDir, "Norm-stats-report")
    plotFrontPage(jobName, currentFont)
    pageNo <- 2
    plotContrastPHists(nst, jobName, currentLayout, pageNo)
    
    grDevices::dev.off()
}

#' Generate P-histograms for each contrast
#' 
#' @param nst NormalyzerDE statistics object.
#' @param jobName Name of processing run.
#' @param currentLayout Layout used for document.
#' @param pageno Current page number.
#' @return None
plotContrastPHists <- function(nst, jobName, currentLayout, pageno) {

    contrastPLists <- nst@pairwiseCompsP
    histPlots <- list()
    
    for (i in 1:length(contrastPLists)) {
        contrast <- names(contrastPLists)[i]
        pVals <- contrastPLists[[contrast]]
        df <- data.frame(pVals=pVals)
        histPlots[[i]] <- ggplot2::ggplot(df) + 
            ggplot2::geom_histogram(ggplot2::aes(pVals), na.rm=TRUE, binwidth=0.01) +
            ggplot2::ggtitle(paste("Contrast:", contrast))
    }
    
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout=currentLayout))
    printPlots(histPlots, "HistPlots", pageno, jobName, currentLayout)  
}


reduceTechnicalReplicates <- function(dataMat, techRepGroups) {
    
    uniqueGroups <- unique(techRepGroups)
    indices <- lapply(uniqueGroups, function(i) { which(techRepGroups %in% i) })
    
    collDataMat <- as.matrix(data.frame(lapply(indices, function(inds) {rowMeans(dataMat[, inds], na.rm = TRUE)})))
    colnames(collDataMat) <- uniqueGroups
    collDataMat
}


reduceDesign <- function(designMat, techRepGroups) {
    
    uniqueGroups <- unique(techRepGroups)
    indices <- lapply(uniqueGroups, function(i) { which(techRepGroups %in% i) })
    collDesignMatList <- lapply(indices, function(inds) { designMat[inds[1],] })
    collDesignMat <- do.call(rbind.data.frame, collDesignMatList)
    collDesignMat
}

