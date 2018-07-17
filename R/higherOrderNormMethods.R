#' Normalize over slices in retention time, using either datapoints within the
#' time interval, or if too few, use neighboring datapoints for normalization
#' 
#' Key variables:
#'     targetSliceIndices - Indices in raw matrix for rows falling within retention time interval
#'     normalizationSliceIndices - Rows used for normalization, can include wider interval than target slice
#'                                 if target slice isn't containing enough data
#'     indicesOfInterest - Target indices within the normalization slice window
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
#' rtNormMat <- getRTNormalizedMatrix(dataMat, retentionTimes, performCyclicLoessNormalization, stepSizeMinutes=1, windowMinCount=100)
getRTNormalizedMatrix <- function(rawMatrix, retentionTimes, normMethod, stepSizeMinutes=1, windowMinCount=100, offset=0) {
    
    sortedRT <- sort(retentionTimes)
    
    startVal <- min(stats::na.omit(retentionTimes))
    endVal <- max(stats::na.omit(retentionTimes))
    rowNumbers <- c()

    if (offset) {
        startVal <- startVal - stepSizeMinutes * offset
    }
    
    processedRows <- matrix(, ncol=ncol(rawMatrix), nrow=0)

    for (windowStart in seq(startVal, endVal, stepSizeMinutes)) {
        
        windowEnd <- windowStart + stepSizeMinutes
        normalizationStartRT <- windowStart
        normalizationEndRT <- windowEnd
        targetSliceIndices <- which(retentionTimes >= windowStart & retentionTimes < windowEnd)
        
        if (length(targetSliceIndices) == 0) {
            # print(paste("No values found, skipping window from", windowStart, "to", windowEnd))
            next
        }
        else if (length(targetSliceIndices) < windowMinCount) {
            
            normalizationRange <- getWidenedRTRange(windowStart, windowEnd, windowMinCount, retentionTimes)
            normalizationStartRT <- normalizationRange[1]
            normalizationEndRT <- normalizationRange[2]
        }
        
        # Likely issue: Retrieving wrong values (although the number seems OK)
        normalizationSliceIndices <- which(retentionTimes >= normalizationStartRT & retentionTimes <= normalizationEndRT)
        normalizationRows <- rawMatrix[normalizationSliceIndices,, drop=FALSE]
        processedNormalizationRows <- normMethod(normalizationRows)
        
        rownames(processedNormalizationRows) <- rownames(normalizationRows)
        
        
        if (length(targetSliceIndices) < length(normalizationSliceIndices)) {
            
            # THIS WILL ALWAYS RETURN LIST OF NUMBERS CORRESPONDING TO LENGTH OF targetSliceIndices !!!
            # indicesOfInterest <- which(targetSliceIndices %in% normalizationSliceIndices)
            indicesOfInterest <- which(normalizationSliceIndices %in% targetSliceIndices)
            normalizedTargetRows <- processedNormalizationRows[indicesOfInterest,]
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
#' 
#' @param rtStart Original retention time start point
#' @param rtEnd Original retention time end point
#' @param minimumDatapoints Required number of datapoints to fullfil
#' @param retentionTimes Vector with all retention times
#' @return Vector with start and end of new RT range
getWidenedRTRange <- function(rtStart, rtEnd, minimumDatapoints, retentionTimes) {
    
    sortedRts <- sort(retentionTimes)
    currentRTSlice <- sortedRts[sortedRts >= rtStart & sortedRts < rtEnd] 
    
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
    
    widenedSlice <- sortedRts[newStartRtIndex:newEndRtIndex]

    if (length(widenedSlice) != minimumDatapoints) {
        stop(paste("Widened slice should equal to minimum number of data points"))
    }
    
    stopifnot(length(widenedSlice) == minimumDatapoints)

    widenedStartRt <- min(widenedSlice)
    widenedEndRt <- max(widenedSlice)
            
    c(widenedStartRt, widenedEndRt)
}


#' Generate two RT time-window normalized matrices where one is shifted.
#' Then return the mean of these matrices.
#' 
#' @param rawMatrix Target matrix to be normalized
#' @param retentionTimes Vector of retention times corresponding to rawMatrix
#' @param normMethod The normalization method to apply to the time windows
#' @param stepSizeMinutes Size of windows to be normalized
#' @param frameShifts Number of frame shifts.
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
#'     frameShifts=3, mergeMethod="median")
getSmoothedRTNormalizedMatrix <- function(rawMatrix, retentionTimes, normMethod, stepSizeMinutes, 
                                          frameShifts=2, windowMinCount=100, mergeMethod="mean") {
    
    matrices <- list()

    for (i in seq_len(frameShifts)) {
        
        fracShift <- (i - 1) * 1 / frameShifts
        matrices[[i]] <- getRTNormalizedMatrix(rawMatrix, retentionTimes, normMethod, 
                                               stepSizeMinutes, windowMinCount=windowMinCount, offset=fracShift)
    }

    if (mergeMethod == "mean") {
        combinedMatrices <- getCombinedMatrix(matrices, mean)
    }
    else if (mergeMethod == "median") {
        combinedMatrices <- getCombinedMatrix(matrices, stats::median)
    }
    else {
        stop(paste("Unknown merge method:", mergeMethod))
    }
    
    colnames(combinedMatrices) <- colnames(rawMatrix)
    combinedMatrices
}

getCombinedMatrix <- function(mList, combFunc) {
    
    matrixCount <- length(mList)
    rows <- nrow(mList[[1]])
    cols <- ncol(mList[[1]])
    mLength <- rows * cols
    combinedMatrix <- matrix(0, nrow=rows, ncol=cols)
    
    for (i in seq_len(mLength)) {
        elemVals <- sapply(mList, function(mat) {mat[[i]]})
        targetVal <- combFunc(elemVals)
        combinedMatrix[i] <- targetVal
    }
    
    combinedMatrix
}







