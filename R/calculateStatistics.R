
#' Filter rows with lower than given number of replicates for any condition
#' 
#' @param df Dataframe with expression data to filter
#' @param groups Condition groups header
#' @param leastRep Minimum number of replicates in each group
#'   to retain
#' @return collDesignDf Reduced design matrix
#' @keywords internal
filterLowRep <- function(df, groups, leastRep = 2) {
    
    allReplicatesHaveValuesContrast <- function(row, groups, minCount) {
        names(row) <- groups
        repCounts <- table(names(stats::na.omit(row)))
        length(repCounts) == length(unique(groups)) && min(repCounts) >= minCount
    }
    
    rowMeetThresContrast <- apply(
        df, 
        1,
        allReplicatesHaveValuesContrast, 
        groups = groups, 
        minCount = leastRep)
    
    filteredDf <- df[rowMeetThresContrast, ]
    filteredDf
}


#' Remove technical replicates from data matrix while collapsing
#' sample values into average values
#' 
#' @param dataMat NormalyzerDE data matrix
#' @param techRepGroups Technical replicates vector
#' @return collDataMat Reduced data matrix
#' @export
#' @examples
#' techRep <- c("a", "a", "b", "b", "c", "c", "d", "d")
#' testData <- data.frame(
#'     c(1,1,1), 
#'     c(1,2,1), 
#'     c(3,3,3), 
#'     c(5,3,3), 
#'     c(5,5,4), 
#'     c(5,5,5), 
#'     c(7,7,7), 
#'     c(7,9,7))
#' colnames(testData) <- c("a1", "a2", "b1", "b2", "c1", "c2", "d1", "d2")
#' statObj <- reduceTechnicalReplicates(testData, techRep)
reduceTechnicalReplicates <- function(dataMat, techRepGroups) {
    
    uniqueGroups <- unique(techRepGroups)
    indices <- lapply(uniqueGroups, function(i) { which(techRepGroups %in% i) })
    
    collDataMat <- as.matrix(data.frame(lapply(
        indices, 
        function(inds) {rowMeans(dataMat[, inds], na.rm = TRUE)})))
    colnames(collDataMat) <- uniqueGroups
    collDataMat
}

#' Remove technical replicates from design matrix.
#' Technical replicates are specified as duplicate strings.
#' The first sample name corresponding for each technical replicate group
#' is retained.
#' 
#' @param designMat NormalyzerDE design matrix
#' @param techRepGroups Technical replicates vector
#' @return collDesignDf Reduced design matrix
#' @export
#' @examples
#' designDf <- data.frame(
#'   sample=c("a1", "a2", "b1", "b2", "c1", "c2", "d1", "d2"),
#'   group=c(rep("A", 4), rep("B", 4)),
#'   techrep=c("a", "a", "b", "b", "c", "c", "d", "d")
#' )
#' statObj <- reduceDesignTechRep(designDf, designDf$techrep)
reduceDesignTechRep <- function(designMat, techRepGroups) {
    
    uniqueGroups <- unique(techRepGroups)
    indices <- lapply(uniqueGroups, function(i) { which(techRepGroups %in% i) })
    collDesignMatList <- lapply(indices, function(inds) { designMat[inds[1],] })
    collDesignDf <- do.call(rbind.data.frame, collDesignMatList)
    collDesignDf <- droplevels(collDesignDf)
    rownames(collDesignDf) <- seq_len(nrow(collDesignDf))
    collDesignDf
}

#' Setup NormalyzerDE statistics object
#' 
#' @param designDf Design matrix.
#' @param fullDf Full data matrix containing expression data and annotation data.
#' @param sampleCol Name of column containing sample IDs.
#' @param conditionCol Name of column containing condition levels.
#' @param batchCol Optional name of column containing batch effect levels.
#' @param logTrans Option for log transforming input data.
#' @param leastRepCount Lowest number of replicate required.
#' @return Prepared statistics object.
#' @export
#' @examples 
#' data(example_design)
#' data(example_stat_data)
#' setupStatisticsObject(example_design, example_stat_data)
setupStatisticsObject <- function(designDf, fullDf, sampleCol="sample", 
                                  conditionCol="group", batchCol=NULL, logTrans=FALSE,
                                  leastRepCount=2) {

    dataCols <- as.character(designDf[, sampleCol])
    annotMat <- as.matrix(fullDf[, which(!colnames(fullDf) %in% dataCols)])
    dataMat <- as.matrix(fullDf[, dataCols], )
    
    rownames(dataMat) <- seq_len(nrow(dataMat))
    dataMatNAFiltered <- filterLowRep(dataMat, designDf[, conditionCol], leastRep=leastRepCount)
    naFilterContrast <- rownames(dataMat) %in% rownames(dataMatNAFiltered)
      
    if (logTrans) {
        dataMat <- log2(dataMat)
    }
    
    nst <- NormalyzerStatistics(
        annotMat=annotMat, 
        dataMat=dataMat, 
        designDf=designDf, 
        filteredDataMat=dataMatNAFiltered,
        filteringContrast=naFilterContrast
    )
    
    nst
}

#' Generate annotated dataframe from statistics object
#' 
#' @param nst NormalyzerDE statistics object.
#' @return outDf Annotated statistics matrix
#' @export
#' @examples
#' data(example_stat_data)
#' data(example_design)
#' statObj <- setupStatisticsObject(example_design, example_stat_data)
#' statObj <- calculateContrasts(statObj, comparisons=c("1-2", "2-3"), condCol="group", type="limma")
#' annotDf <- generateAnnotatedMatrix(statObj)
generateAnnotatedMatrix <- function(nst) {
    
    pairwiseHead <- paste(
        names(nst@pairwiseCompsP), 
        "PValue", 
        sep="_")
    pMat <- data.frame(nst@pairwiseCompsP)
    colnames(pMat) <- pairwiseHead
    
    pairwiseHeadFdr <- paste(
        names(nst@pairwiseCompsFdr), 
        "AdjPVal", 
        sep="_")
    fdrMat <- data.frame(nst@pairwiseCompsFdr)
    colnames(fdrMat) <- pairwiseHeadFdr
    
    pairwiseHeadFold <- paste(
        names(nst@pairwiseCompsFold), 
        "log2FoldChange", 
        sep="_")
    foldMat <- data.frame(nst@pairwiseCompsFold)
    colnames(foldMat) <- pairwiseHeadFold
    
    aveMat <- data.frame(nst@pairwiseCompsAve)
    outDf <- cbind(nst@annotMat, pMat, fdrMat, foldMat, 
                   featureAvg=aveMat[, 1], nst@dataMat)
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
#' data(example_stat_data)
#' data(example_design)
#' statObj <- setupStatisticsObject(example_design, example_stat_data)
#' statObj <- calculateContrasts(statObj, comparisons=c("1-2", "2-3"), 
#'   condCol="group", type="limma")
#' outputDir <- tempdir()
#' generateStatsReport(statObj, "jobName", outputDir)
generateStatsReport <- function(nst, jobName, jobDir, plotRows=3, plotCols=4) {
    
    nrows <- plotRows + 2
    ncols <- plotCols + 2
    
    currentLayout <- grid::grid.layout(
        nrow=nrows, ncol=ncols,
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
#' @keywords internal
plotContrastPHists <- function(nst, jobName, currentLayout, pageno) {

    contrastPLists <- nst@pairwiseCompsP
    histPlots <- list()
    
    for (i in seq_len(length(contrastPLists))) {
        contrast <- names(contrastPLists)[i]
        pVals <- contrastPLists[[contrast]]
        df <- data.frame(pVals=pVals)
        histPlots[[i]] <- ggplot2::ggplot(df) + 
            ggplot2::geom_histogram(
                ggplot2::aes(pVals), na.rm=TRUE, binwidth=0.01) +
            ggplot2::ggtitle(paste("Contrast:", contrast))
    }
    
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout=currentLayout))
    printPlots(histPlots, "HistPlots", pageno, jobName, currentLayout)  
}


