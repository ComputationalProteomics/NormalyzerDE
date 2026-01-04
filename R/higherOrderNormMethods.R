#' Perform RT-segmented normalization by performing the supplied normalization
#' over retention-time sliced data
#' 
#' The function orders the retention times and steps through them using the
#' supplied step size (in minutes). If smaller than a fixed lower boundary 
#' the window is expanded to ensure a minimum amount of data in each 
#' normalization step. An offset can be specified which can be used to perform
#' multiple RT-segmentations with partial overlapping windows.
#' 
#' @param rawMatrix Target matrix to be normalized
#' @param retentionTimes Vector of retention times corresponding to rawMatrix
#' @param normMethod The normalization method to apply to the time windows
#' @param stepSizeMinutes Size of windows to be normalized
#' @param windowMinCount Minimum number of values for window to not be expanded.
#' @param offset Whether time window should shifted half step size
#' @param noLogTransform Don't log-transform the data
#' @return Normalized matrix
#' @export
#' @examples
#' data(example_data_small)
#' data(example_design_small)
#' data(example_data_only_values)
#' dataMat <- example_data_only_values
#' retentionTimes <- as.numeric(example_data[, "Average.RT"])
#' performCyclicLoessNormalization <- function(rawMatrix) {
#'     log2Matrix <- log2(rawMatrix)
#'     normMatrix <- limma::normalizeCyclicLoess(log2Matrix, method="fast")
#'     colnames(normMatrix) <- colnames(rawMatrix)
#'     normMatrix
#' }
#' rtNormMat <- getRTNormalizedMatrix(dataMat, retentionTimes, 
#' performCyclicLoessNormalization, stepSizeMinutes=1, windowMinCount=100)
getRTNormalizedMatrix <- function(rawMatrix, retentionTimes, normMethod, 
                                  stepSizeMinutes=1, windowMinCount=100, 
                                  offset=0, noLogTransform=FALSE) {
    
    # Key variables:
    #    targetSliceIndices 
    #        Indices in raw matrix for rows falling within retention time 
    #        interval
    #    normalizationSliceIndices 
    #        Rows used for normalization, can include wider interval than 
    #        target slice if target slice isn't containing enough data
    #    indicesOfInterest
    #        Target indices within the normalization slice window
    
    if (!is(rawMatrix, "matrix")) {
        stop("Type of rawMatrix is expected to be matrix, received: ", 
             class(rawMatrix)
        )
    }
    
    sortedRetentionTimes <- sort(retentionTimes)

    startVal <- min(retentionTimes, na.rm=TRUE)
    endVal <- max(retentionTimes, na.rm=TRUE)

    if (offset) {
        startVal <- startVal - stepSizeMinutes * offset
    }
    
    processedRowsList <- list()
    rowNumbersList <- list()
    sliceIndex <- 0

    for (windowStart in seq(startVal, endVal, stepSizeMinutes)) {
        
        windowEnd <- windowStart + stepSizeMinutes
        normalizationStartRT <- windowStart
        normalizationEndRT <- windowEnd
        targetSliceIndices <- which(
            retentionTimes >= windowStart & retentionTimes < windowEnd
        )
        
        if (length(targetSliceIndices) == 0) {
            next
        }
        else if (length(targetSliceIndices) < windowMinCount) {
            
            normalizationRange <- getWidenedRTRange(
                windowStart, 
                windowEnd, 
                windowMinCount, 
                retentionTimes,
                sortedRetentionTimes=sortedRetentionTimes
            )
            normalizationStartRT <- normalizationRange[1]
            normalizationEndRT <- normalizationRange[2]
        }
        
        normalizationSliceIndices <- which(
            retentionTimes >= normalizationStartRT & 
                retentionTimes <= normalizationEndRT
        )
        normalizationRows <- rawMatrix[normalizationSliceIndices,, drop=FALSE]
        
        if (noLogTransform) {
            processedNormalizationRows <- normMethod(normalizationRows, noLogTransform=noLogTransform)
        }
        else {
            processedNormalizationRows <- normMethod(normalizationRows)
        }
        
        rownames(processedNormalizationRows) <- rownames(normalizationRows)
        
        if (length(targetSliceIndices) < length(normalizationSliceIndices)) {
            
            indicesOfInterest <- which(
                normalizationSliceIndices %in% targetSliceIndices
            )
            normalizedTargetRows <- processedNormalizationRows[indicesOfInterest,, drop=FALSE]
        }
        else {
            normalizedTargetRows <- processedNormalizationRows
        }
        
        sliceIndex <- sliceIndex + 1
        rowNumbersList[[sliceIndex]] <- targetSliceIndices
        processedRowsList[[sliceIndex]] <- normalizedTargetRows
    }
    
    if (sliceIndex == 0) {
        return(matrix(, ncol=ncol(rawMatrix), nrow=0))
    }

    rowNumbers <- unlist(rowNumbersList, use.names=FALSE)
    processedRows <- do.call(rbind, processedRowsList)

    orderedProcessedRows <- processedRows[order(rowNumbers), , drop=FALSE]
    orderedProcessedRows
}

#' Pick datapoints before and after window until a minimum number is reached
#' Expects the start and end retention times to match actual retention times
#' present in the data
#' 
#' @param rtStart Original retention time start point
#' @param rtEnd Original retention time end point
#' @param minimumDatapoints Required number of datapoints to fulfill
#' @param retentionTimes Vector with all retention times
#' @return Vector with start and end of new RT range
#' @keywords internal
getWidenedRTRange <- function(rtStart, rtEnd, minimumDatapoints, retentionTimes,
                              sortedRetentionTimes=NULL, allowTooWideData=FALSE) {
    
    sortedRts <- if (is.null(sortedRetentionTimes)) {
        sort(retentionTimes)
    }
    else {
        sortedRetentionTimes
    }
    currentRTSlice <- sortedRts[sortedRts >= rtStart & sortedRts < rtEnd] 
    
    if (length(currentRTSlice) == 0) {
        stop("Selected retention time slice doesn't contain any data")
    }
    
    if (length(currentRTSlice) > minimumDatapoints) {
        if (allowTooWideData) {
            return(c(rtStart, rtEnd))
        }
        else {
            stop("Number of datapoints exceed minimum, add option ",
                 "'allowTooWideData' to process anyway")
        }
    }
    
    # Get single element if multiple with exactly same RT
    startIndex <- utils::head(which(sortedRts == min(currentRTSlice)), 1)
    endIndex <- utils::tail(which(sortedRts == max(currentRTSlice)), 1)

    currentCount <- length(currentRTSlice)
    remainingCount <- minimumDatapoints - currentCount
    
    pickBefore <- floor(remainingCount / 2)
    pickAfter <- ceiling(remainingCount / 2)
    
    totalBefore <- length(sortedRts[sortedRts < rtStart])
    totalAfter <- length(sortedRts[sortedRts >= rtEnd])
    
    stopifnot(remainingCount == pickBefore + pickAfter)
    stopifnot(totalBefore + totalAfter + length(currentRTSlice) == length(retentionTimes))
    
    if (pickBefore > totalBefore && pickAfter > totalAfter) {
        stop("Not enough values in dataset to do RT normalization with current 
             minimum datapoints setting - Please adjust settings")
    }
    else if (pickBefore > totalBefore) {
        diff <- pickBefore - totalBefore
        pickAfter <- pickAfter + diff
        pickBefore <- pickBefore - diff
    }
    else if (pickAfter > totalAfter) {
        diff <- pickAfter - totalAfter
        pickBefore <- pickBefore + diff
        pickAfter <- pickAfter - diff
    }
    
    newStartRtIndex <- startIndex - pickBefore
    newEndRtIndex <- endIndex + pickAfter
    
    if (newEndRtIndex - newStartRtIndex + 1 > length(sortedRts)) {
        stop("Requested minimum window size (", 
             newEndRtIndex - newStartRtIndex + 1, 
             ") exceeds total number of datapoints (", 
             length(sortedRts),
             ")")
    }
    
    widenedSlice <- sortedRts[newStartRtIndex:newEndRtIndex]

    if (length(widenedSlice) != minimumDatapoints) {
        stop("Widened slice should equal to minimum number of data points")
    }
    
    stopifnot(length(widenedSlice) == minimumDatapoints)

    widenedStartRt <- min(widenedSlice)
    widenedEndRt <- max(widenedSlice)
            
    c(widenedStartRt, widenedEndRt)
}


#' Generate multiple RT time-window normalized matrices where one is shifted.
#' Merge them using a specified method (mean or median) and return the result.
#' 
#' Uses the function getRTNormalizedMatrix to generate multiple normalized
#' matrices which are shifted respective to each other and finally merged into
#' a single matrix. This could potentially reduce effect of fluctuations
#' within individual windows.
#' 
#' @param rawMatrix Target matrix to be normalized
#' @param retentionTimes Vector of retention times corresponding to rawMatrix
#' @param normMethod The normalization method to apply to the time windows
#' @param stepSizeMinutes Size of windows to be normalized
#' @param windowShifts Number of frame shifts.
#' @param windowMinCount Minimum number of features within window.
#' @param mergeMethod Layer merging approach. Mean or median.
#' @param noLogTransform Don't log transform the input
#' @return Normalized matrix
#' @export
#' @examples
#' 
#' data(example_data_small)
#' data(example_data_only_values)
#' data(example_design_small)
#' retentionTimes <- as.numeric(example_data[, "Average.RT"])
#' dataMat <- example_data_only_values
#' performCyclicLoessNormalization <- function(rawMatrix) {
#'     log2Matrix <- log2(rawMatrix)
#'     normMatrix <- limma::normalizeCyclicLoess(log2Matrix, method="fast")
#'     colnames(normMatrix) <- colnames(rawMatrix)
#'     normMatrix
#' }
#' rtNormMat <- getSmoothedRTNormalizedMatrix(dataMat, retentionTimes, 
#'     performCyclicLoessNormalization, stepSizeMinutes=1, windowMinCount=100, 
#'     windowShifts=2, mergeMethod="median")
getSmoothedRTNormalizedMatrix <- function(
    rawMatrix, retentionTimes, normMethod, stepSizeMinutes, 
    windowShifts=2, windowMinCount=100, mergeMethod="mean", noLogTransform=FALSE) {
    
    matrices <- list()
    
    for (i in seq_len(windowShifts)) {
        
        fracShift <- (i - 1) * 1 / windowShifts
        matrices[[i]] <- getRTNormalizedMatrix(
            rawMatrix, 
            retentionTimes, 
            normMethod,
            stepSizeMinutes=stepSizeMinutes, 
            windowMinCount=windowMinCount, 
            offset=fracShift,
            noLogTransform=noLogTransform
        )
    }

    if (mergeMethod == "mean") {
        combinedMatrices <- getCombinedMatrix(matrices, mean)
    }
    else if (mergeMethod == "median") {
        combinedMatrices <- getCombinedMatrix(matrices, stats::median)
    }
    else {
        stop("Unknown merge method: ", mergeMethod)
    }
    
    colnames(combinedMatrices) <- colnames(rawMatrix)
    combinedMatrices
}

#' Merge multiple dataframes using provided function
#' 
#' @param mList List containing dataframes of same shape
#' @param combFunc Function performing elementwise merge of matrices
#' @return combinedMatrix A single dataframe with combined data
#' @keywords internal
getCombinedMatrix <- function(mList, combFunc) {
    
    matrixCount <- length(mList)

    if (matrixCount == 0) {
        stop("Expected at least one matrix to merge")
    }
    if (matrixCount == 1) {
        return(mList[[1]])
    }

    rows <- nrow(mList[[1]])
    cols <- ncol(mList[[1]])

    if (any(vapply(mList, function(mat) { !all(dim(mat) == c(rows, cols)) }, logical(1)))) {
        stop("All matrices must have the same dimensions to merge")
    }

    if (identical(combFunc, mean)) {
        combinedMatrix <- Reduce("+", mList) / matrixCount
        colnames(combinedMatrix) <- colnames(mList[[1]])
        return(combinedMatrix)
    }

    isMedian <- identical(combFunc, stats::median) || identical(combFunc, median)
    if (!isMedian) {
        stop("Unknown merge function. Only mean and median are supported.")
    }

    if (matrixCount == 2) {
        combinedMatrix <- (mList[[1]] + mList[[2]]) / 2
        colnames(combinedMatrix) <- colnames(mList[[1]])
        return(combinedMatrix)
    }

    if (matrixCount == 3) {
        m1 <- mList[[1]]
        m2 <- mList[[2]]
        m3 <- mList[[3]]
        combinedMatrix <- m1 + m2 + m3 - pmin(m1, m2, m3) - pmax(m1, m2, m3)
        colnames(combinedMatrix) <- colnames(mList[[1]])
        return(combinedMatrix)
    }

    stackedValues <- do.call(cbind, lapply(mList, as.vector))
    medians <- matrixStats::rowMedians(stackedValues, na.rm=FALSE)

    combinedMatrix <- matrix(medians, nrow=rows, ncol=cols)
    colnames(combinedMatrix) <- colnames(mList[[1]])
    combinedMatrix
}

