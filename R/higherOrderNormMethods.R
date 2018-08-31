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
#' @return Normalized matrix
#' @export
#' @examples
#' data(example_stat_data)
#' data(example_design)
#' sampleNames <- as.character(example_design$sample)
#' dataMat <- as.matrix(example_stat_data[, sampleNames])
#' retentionTimes <- example_stat_data$Average.RT
#' performCyclicLoessNormalization <- function(rawMatrix) {
#'     log2Matrix <- log2(rawMatrix)
#'     normMatrix <- limma::normalizeCyclicLoess(log2Matrix, method="fast")
#'     colnames(normMatrix) <- colnames(rawMatrix)
#'     normMatrix
#' }
#' rtNormMat <- getRTNormalizedMatrix(dataMat, retentionTimes, 
#'   performCyclicLoessNormalization, stepSizeMinutes=1, windowMinCount=100)
getRTNormalizedMatrix <- function(rawMatrix, retentionTimes, normMethod, 
                                  stepSizeMinutes=1, windowMinCount=100, 
                                  offset=0) {
    
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
    
    sortedRT <- sort(retentionTimes)
    
    startVal <- min(retentionTimes, na.rm=TRUE)
    endVal <- max(retentionTimes, na.rm=TRUE)
    rowNumbers <- c()

    if (offset) {
        startVal <- startVal - stepSizeMinutes * offset
    }
    
    processedRows <- matrix(, ncol=ncol(rawMatrix), nrow=0)

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
                retentionTimes
            )
            normalizationStartRT <- normalizationRange[1]
            normalizationEndRT <- normalizationRange[2]
        }
        
        normalizationSliceIndices <- which(
            retentionTimes >= normalizationStartRT & 
                retentionTimes <= normalizationEndRT
        )
        normalizationRows <- rawMatrix[normalizationSliceIndices,, drop=FALSE]
        processedNormalizationRows <- normMethod(normalizationRows)
        
        rownames(processedNormalizationRows) <- rownames(normalizationRows)
        
        
        if (length(targetSliceIndices) < length(normalizationSliceIndices)) {
            
            indicesOfInterest <- which(
                normalizationSliceIndices %in% targetSliceIndices
            )
            normalizedTargetRows <- processedNormalizationRows[indicesOfInterest,]
            # normalizedTargetRows <- processedNormalizationRows[indicesOfInterest,, drop=FALSE]
        }
        else {
            normalizedTargetRows <- processedNormalizationRows
        }
        
        rowNumbers <- c(rowNumbers, targetSliceIndices)
        processedRows <- rbind(processedRows, normalizedTargetRows)
    }
    
    orderedProcessedRows <- processedRows[order(rowNumbers), ]
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
                              allowTooWideData=FALSE) {
    
    sortedRts <- sort(retentionTimes)
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
    startIndex <- utils::tail(which(sortedRts == min(currentRTSlice)), 1)
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
#' @return Normalized matrix
#' @export
#' @examples
#' 
#' data(example_stat_data)
#' data(example_design)
#' sampleNames <- as.character(example_design$sample)
#' dataMat <- as.matrix(example_stat_data[, sampleNames])
#' retentionTimes <- example_stat_data$Average.RT
#' performCyclicLoessNormalization <- function(rawMatrix) {
#'     log2Matrix <- log2(rawMatrix)
#'     normMatrix <- limma::normalizeCyclicLoess(log2Matrix, method="fast")
#'     colnames(normMatrix) <- colnames(rawMatrix)
#'     normMatrix
#' }
#' rtNormMat <- getSmoothedRTNormalizedMatrix(dataMat, retentionTimes, 
#'     performCyclicLoessNormalization, stepSizeMinutes=1, windowMinCount=100, 
#'     windowShifts=3, mergeMethod="median")
getSmoothedRTNormalizedMatrix <- function(
    rawMatrix, retentionTimes, normMethod, stepSizeMinutes, 
    windowShifts=2, windowMinCount=100, mergeMethod="mean") {
    
    matrices <- list()
    
    for (i in seq_len(windowShifts)) {
        
        fracShift <- (i - 1) * 1 / windowShifts
        matrices[[i]] <- getRTNormalizedMatrix(
            rawMatrix, 
            retentionTimes, 
            normMethod,
            stepSizeMinutes=stepSizeMinutes, 
            windowMinCount=windowMinCount, 
            offset=fracShift)
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
    rows <- nrow(mList[[1]])
    cols <- ncol(mList[[1]])
    mLength <- rows * cols
    combinedMatrix <- matrix(0, nrow=rows, ncol=cols)
    
    # Iterate over each element position
    for (i in seq_len(mLength)) {
        elemVals <- vapply(mList, function(mat) {mat[[i]]}, 0)
        targetVal <- combFunc(elemVals)
        combinedMatrix[i] <- targetVal
    }
    
    colnames(combinedMatrix) <- colnames(mList[[1]])
    combinedMatrix
}



