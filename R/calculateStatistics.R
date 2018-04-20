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
calculateStatistics <- function(dataFp, designFp, comparisons, limmaTest=TRUE, varFilterFrac=1, 
                                logTrans=FALSE, robustLimma=FALSE, type="limma", batchCol=NULL) {
    
    nst <- setupStatisticsObject(dataFp, designFp, comparisons, logTrans=logTrans)
    
    nst <- calculateContrasts(nst, comparisons, condCol="group", type=type, batchCol=batchCol)

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
    
    nst <- NormalyzerStatistics(annotMat=annotMat, dataMat=dataMat, designDf=designDf, contrasts=comparisons)
    nst
}

generateAnnotatedMatrix <- function(nst) {
    
    pairwiseHead <- paste(names(nst@pairwiseCompsP), "p", sep="_")
    pMat <- data.frame(nst@pairwiseCompsP)
    colnames(pMat) <- pairwiseHead
    
    pairwiseHeadFdr <- paste(names(nst@pairwiseCompsFdr), "p.adj", sep="_")
    fdrMat <- data.frame(nst@pairwiseCompsFdr)
    colnames(fdrMat) <- pairwiseHeadFdr
    
    pairwiseHeadFold <- paste(names(nst@pairwiseCompsFold), "fold", sep="_")
    foldMat <- data.frame(nst@pairwiseCompsFold)
    colnames(foldMat) <- pairwiseHeadFold
    
    aveMat <- data.frame(nst@pairwiseCompsAve)
    
    outDf <- cbind(nst@annotMat, pMat, fdrMat, foldMat, featureAvg=aveMat[, 1], nst@dataMat)
    outDf
}

#' Generate full output report plot document
#' 
#' @param nr Normalyzer results object.
#' @param jobdir Path to output directory for run.
#' @param plot_rows Number of plot rows.
#' @param plot_cols Number of plot columns.
#' @return None
#' @export
#' @examples
#' normObj <- getVerifiedNormalyzerObject("data.tsv", "job_name", "design.tsv")
#' normResults <- normMethods(normObj)
#' normResultsWithEval <- analyzeNormalizations(normResults)
#' generatePlots(normResultsWithEval, "path/to/output")
generateStatsReport <- function(nst, jobName, jobDir, plot_rows=3, plot_cols=4) {
    
    nrows <- plot_rows + 2
    ncols <- plot_cols + 2
    
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
#' @param nr Normalyzer results object.
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



reduce_technical_replicates <- function(dataMat, techRepGroups) {
    
    uniqueGroups <- unique(techRepGroups)
    indices <- lapply(uniqueGroups, function(i) { which(techRepGroups %in% i) })
    
    collDataMat <- as.matrix(data.frame(lapply(indices, function(inds) {rowMeans(dataMat[, inds])})))
    colnames(collDataMat) <- uniqueGroups
    collDataMat
}

reduce_design <- function(designMat, techRepGroups) {
    
    uniqueGroups <- unique(techRepGroups)
    indices <- lapply(uniqueGroups, function(i) { which(techRepGroups %in% i) })
    collDesignMatList <- lapply(indices, function(inds) { designMat[inds[1],] })
    collDesignMat <- do.call(rbind.data.frame, collDesignMatList)
    collDesignMat
}

# reduce_technical_replicates = function(data_m, sample_names, sample_numbers) {
#     
#     merge_df <- data.frame(matrix(NA, nrow=nrow(data_m), ncol=0))
#     unique_samples <- unique(sample_numbers)
#     
#     for (sample_nbr in unique_samples[!is.na(unique_samples)]) {
#         sample_indices <- which(!is.na(match(sample_numbers, sample_nbr)))
#         target_names <- sample_names[sample_indices]
#         data_cols <- data_m[, target_names, drop=F]
#         if (ncol(data_cols) > 1) {
#             mean_col <- data.frame(apply(data_cols, 1, function(row) { 
#                 row_mean <- mean(na.omit(row)) 
#                 if (!is.nan(row_mean)) {
#                     row_mean
#                 }
#                 else {
#                     NA
#                 }
#             }))
#             colnames(mean_col) <- target_names[1]
#             merge_df <- cbind(merge_df, mean_col)
#         }
#         else {
#             merge_df <- cbind(merge_df, data_cols)
#         }
#     }
#     
#     # browser()
#     
#     # Include samples if any started as NA
#     # no_sample_nbr_df <- data_m[which(is.na(sample_numbers))]
#     # merge_m <- cbind(no_sample_nbr_df, merge_df)
#     
#     # Reorder to maintain original order
#     # orig <- colnames(data_m)
#     # updated <- colnames(merge_m)
#     # new_names_orig_order <- orig[sort(match(updated, orig))]
#     # merge_m <- merge_m[, new_names_orig_order]
#     
#     merge_m
# }

