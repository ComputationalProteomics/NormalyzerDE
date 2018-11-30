#' Remove technical replicates from data matrix 
#' 
#' Collapses sample values into their average. If only one value is present
#' due to NA-values in other technical replicates, then that value is used.
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
#' 
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

#' Generate an annotated data frame from statistics object
#' 
#' Extracts key values (p-value, adjusted p-value, log2-fold change and
#' average expression values) from an NormalyzerStatistics instance
#' and appends these to the annotation- and data-matrices
#' 
#' @param nst NormalyzerDE statistics object.
#' @return outDf Annotated statistics matrix
#' @export
#' @examples
#' data(example_stat_summarized_experiment)
#' statObj <- NormalyzerStatistics(example_stat_summarized_experiment)
#' statObj <- calculateContrasts(statObj, comparisons=c("1-2", "2-3"), condCol="group", type="limma")
#' annotDf <- generateAnnotatedMatrix(statObj)
generateAnnotatedMatrix <- function(nst) {
    
    pairwiseHead <- paste(
        names(pairwiseCompsP(nst)), "PValue", sep="_"
    )
    pMat <- data.frame(pairwiseCompsP(nst))
    colnames(pMat) <- pairwiseHead
    
    pairwiseHeadFdr <- paste(
        names(pairwiseCompsFdr(nst)), "AdjPVal", sep="_"
    )
    fdrMat <- data.frame(pairwiseCompsFdr(nst))
    colnames(fdrMat) <- pairwiseHeadFdr
    
    pairwiseHeadFold <- paste(
        names(pairwiseCompsFold(nst)), "log2FoldChange", sep="_"
    )
    foldMat <- data.frame(pairwiseCompsFold(nst))
    colnames(foldMat) <- pairwiseHeadFold
    
    aveMat <- data.frame(pairwiseCompsAve(nst))
    outDf <- cbind(annotMat(nst), pMat, fdrMat, foldMat, 
                   featureAvg=aveMat[, 1], dataMat(nst))
    outDf
}

#' Generate full output report plot document. Plots p-value histograms
#' for each contrast in the NormalyzerStatistics instance and writes these
#' to a PDF report.
#' 
#' @param nst NormalyzerDE statistics object.
#' @param jobName Name of processing run.
#' @param jobDir Path to output directory.
#' @param plotRows Number of plot rows.
#' @param plotCols Number of plot columns.
#' @return None
#' @export
#' @examples
#' data(example_stat_summarized_experiment)
#' statObj <- NormalyzerStatistics(example_stat_summarized_experiment)
#' statObj <- calculateContrasts(statObj, comparisons=c("1-2", "2-3"), 
#'   condCol="group", type="limma")
#' outputDir <- tempdir()
#' generateStatsReport(statObj, "jobName", outputDir)
generateStatsReport <- function(nst, jobName, jobDir, plotRows=3, plotCols=4) {
    
    nrows <- plotRows + 1
    ncols <- plotCols + 1
    
    currentLayout <- grid::grid.layout(
        nrow=nrows, 
        ncol=ncols,
        heights=c(0.1, rep(3 / (nrows - 2), (nrows - 2)), 0.1), 
        widths=c(0.1, rep(4 / (ncols - 2), (ncols - 2)), 0.1), 
        default.units=c("null", "null"))
    
    currentFont <- "Helvetica"
    setupPlotting(jobName, jobDir, "Norm-stats-report")
    plotFrontPage(jobName, currentFont)
    
    pageNo <- 2
    plotContrastPHists(nst, jobName, currentLayout, pageNo)
    
    pageNo <- 3
    plotSigScatter(nst, jobName, currentLayout, pageNo, type="MA")

    pageNo <- 4
    plotSigScatter(nst, jobName, currentLayout, pageNo, type="Vulcano")

    grDevices::dev.off()
}

#' Takes an NormalyzerStatistics instance and generates and prints a p-value
#' histogram for each onto the viewport
#' 
#' @param nst NormalyzerDE statistics object.
#' @param jobName Name of processing run.
#' @param currentLayout Layout used for document.
#' @param pageno Current page number.
#' @return None
#' @keywords internal
plotContrastPHists <- function(nst, jobName, currentLayout, pageno) {

    contrastPLists <- pairwiseCompsP(nst)
    histPlots <- list()
    
    for (i in seq_along(contrastPLists)) {
        contrast <- names(contrastPLists)[i]
        pVals <- contrastPLists[[contrast]]
        df <- data.frame(pVals=pVals)
        histPlots[[i]] <- ggplot2::ggplot(df, ggplot2::aes(pVals)) + 
            ggplot2::theme_classic() + 
            ggplot2::geom_histogram(na.rm=TRUE, binwidth=0.01, color="#268BD2", fill="#268BD2") +
            ggplot2::ggtitle(paste("Contrast:", contrast)) + 
            ggplot2::xlab("P-values") + 
            ggplot2::ylab("Count") +
            ggplot2::theme(
                axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust=0.5), 
                legend.position="none"
            )
    }
    
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout=currentLayout))
    printPlots(histPlots, "HistPlots", pageno, jobName, currentLayout)  
}

#' Takes an NormalyzerStatistics instance and generates and prints a 
#' vulcano plot
#' 
#' @param nst NormalyzerDE statistics object.
#' @param jobName Name of processing run.
#' @param currentLayout Layout used for document.
#' @param pageno Current page number.
#' @param type Specify whether to plot 'Vulcano' or 'MA'.
#' @param sigThres FDR threshold for DE coloring.
#' @return None
#' @keywords internal
plotSigScatter <- function(nst, jobName, currentLayout, pageno, type="Vulcano", sigThres=0.1) {
    
    contrastPLists <- pairwiseCompsP(nst)
    contrastFDRLists <- pairwiseCompsFdr(nst)
    contrastFoldLists <- pairwiseCompsFold(nst)
    contrastAveLists <- pairwiseCompsAve(nst)
    plots <- list()
    
    for (i in seq_along(contrastPLists)) {
        
        contrast <- names(contrastPLists)[i]
        pVals <- contrastPLists[[contrast]]
        fold <- contrastFoldLists[[contrast]]
        expression <- contrastAveLists[[contrast]]
        sig <- contrastFDRLists[[contrast]] < sigThres
        
        df <- data.frame(pVals=pVals, fold=fold, sig=sig, expression=expression)
        df <- df[complete.cases(df),]
        
        if (type == "Vulcano") {
            plots[[i]] <- ggplot2::ggplot(df, ggplot2::aes(fold, -log10(pVals), color=sig)) + 
                ggplot2::geom_point(size=1, alpha=0.5, na.rm=TRUE) +
                ggplot2::xlab("Fold") + 
                ggplot2::ylab("log10(Pvalue)")
        }
        else if (type == "MA") {
            plots[[i]] <- ggplot2::ggplot(df, ggplot2::aes(expression, fold, color=sig)) + 
                ggplot2::geom_point(size=1, alpha=0.5, na.rm=TRUE) +
                ggplot2::xlab("Expression") + 
                ggplot2::ylab("Fold")
        }
        else {
            stop("Unknown plot type: ", type)
        }
        
        plots[[i]] <- plots[[i]] + 
            ggplot2::ggtitle(paste("Contrast:", contrast)) + 
            ggplot2::theme_classic() +
            ggplot2::labs(color=paste("FDR <", sigThres))
    }
    
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout=currentLayout))
    printPlots(plots, type, pageno, jobName, currentLayout)  
}


