#' Normalize over slices in retention time, using either datapoints within the
#' time interval, or if too few, use neighboring datapoints for normalization
#' 
#' @param rawMatrix Target matrix to be normalized
#' @param retentionTimes Vector of retention times corresponding to rawMatrix
#' @param normMethod The normalization method to apply to the time windows
#' @param stepSizeMinutes Size of windows to be normalized
#' @param offset Whether time window should shifted half step size
#' @return Normalized matrix
getRTNormalizedMatrix <- function(rawMatrix, retentionTimes, normMethod, stepSizeMinutes, windowMinCount=50, offset=0) {
    
    # targetSliceIndices - Indices in raw matrix for rows falling within retention time interval
    # normalizationSliceIndices - Rows used for normalization, can include wider interval than target slice
    #                              if target slice isn't containing enough data
    # indicesOfInterest - Target indices within the normalization slice window
    
    sortedRT <- sort(retentionTimes)
    
    startVal <- min(na.omit(retentionTimes))
    endVal <- max(na.omit(retentionTimes))
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
        
        normalizationSliceIndices <- which(retentionTimes >= normalizationStartRT & retentionTimes <= normalizationEndRT)
        normalizationRows <- rawMatrix[normalizationSliceIndices,, drop=FALSE]
        processedNormalizationRows <- normMethod(normalizationRows)
        
        if (length(targetSliceIndices) < length(normalizationSliceIndices)) {
            indicesOfInterest <- which(targetSliceIndices %in% normalizationSliceIndices)
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
    
    startIndex <- which(sortedRts == min(currentRTSlice))
    endIndex <- which(sortedRts == max(currentRTSlice))

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

    stopifnot(length(widenedSlice) == minimumDatapoints)

    widenedStartRt <- min(widenedSlice)
    widenedEndRt <- max(widenedSlice)
            
    c(widenedStartRt, widenedEndRt)
}


# getRTAndFixNormalizedMatrix <- function(rawMatrix, retentionTimes, normMethod, stepSizeMinutes, windowMinCount=1, offset=0) {
#     
#     sortedRT <- sort(retentionTimes)
#     
#     startVal <- min(na.omit(retentionTimes))
#     endVal <- max(na.omit(retentionTimes))
#     rowNumbers <- c()
#     
#     if (offset) {
#         startVal <- startVal - stepSizeMinutes * offset
#     }
#     
#     processedRows <- matrix(, ncol=ncol(rawMatrix), nrow=0)
#     windowStart <- startVal
#     windowEnd <- startVal + stepSizeMinutes
#     
#     while (windowEnd < endVal) {
#         
#         windowEnd <- windowStart + stepSizeMinutes
#         firstTimeSlicedIndices <- which(retentionTimes >= windowStart & retentionTimes < windowEnd)
#         
#         if (length(firstTimeSlicedIndices) < windowMinCount) {
#             
#             startPos <- which(sortedRT >= windowStart)[1]
#             if (startPos + windowMinCount < endVal) {
#                 endPos <- startPos + windowMinCount
#             }
#             else {
#                 endPos <- length(retentionTimes)
#             }
#             
#             windowStart <- sortedRT[startPos]
#             windowEnd <- sortedRT[endPos]
#             nextWindowStart <- sortedRT[endPos + 1]
#             
#             fixSlicedIndices <- c(startPos:endPos)
#             posSliced <- sortedRT[fixSlicedIndices]
#             timeSlicedIndices <- match(posSliced, retentionTimes)
#         }
#         else {
#             timeSlicedIndices <- firstTimeSlicedIndices
#             nextWindowStart <- windowEnd
#         }
#         
#         rowNumbers <- c(rowNumbers, timeSlicedIndices)
#         currentRows <- rawMatrix[timeSlicedIndices,, drop=FALSE]
#         
#         processedSlice <- normMethod(currentRows)
#         processedRows <- rbind(processedRows, processedSlice)
#         
#         windowStart <- nextWindowStart
#     }
#     
#     orderedProcessedRows <- processedRows[order(rowNumbers), ]
#     orderedProcessedRows
# }


#' Generate two RT time-window normalized matrices where one is shifted.
#' Then return the mean of these matrices.
#' 
#' @param rawMatrix Target matrix to be normalized
#' @param retentionTimes Vector of retention times corresponding to rawMatrix
#' @param normMethod The normalization method to apply to the time windows
#' @param stepSizeMinutes Size of windows to be normalized
#' @param offset Whether time window should shifted half step size
#' @return Normalized matrix
getSmoothedRTNormalizedMatrix <- function(rawMatrix, retentionTimes, normMethod, stepSizeMinutes, 
                                          frame_shifts=2, win_size_min=50, merge_method="mean", verbose=FALSE) {
    
    matrices <- list()

    if (verbose) {
        print(paste("win_size_min", win_size_min, "frame_shifts", frame_shifts, "merge_method", merge_method))
    }
    
    for (i in 1:frame_shifts) {
        frac_shift <- (i - 1) * 1 / frame_shifts
        matrices[[i]] <- getRTNormalizedMatrix(rawMatrix, retentionTimes, normMethod, 
                                               stepSizeMinutes, windowMinCount=win_size_min, offset=frac_shift)
    }
    
    # mean, median, anonymous...
    if (merge_method == "mean") {
        combinedMatrices <- getCombinedMatrix(matrices, mean)
    }
    else if (merge_method == "median") {
        combinedMatrices <- getCombinedMatrix(matrices, median)
    }
    else {
        stop(paste("Unknown merge method:", merge_method))
    }
    
    combinedMatrices
}

getCombinedMatrix <- function(m_list, comb_func) {
    
    matrix_count <- length(m_list)
    rows <- nrow(m_list[[1]])
    cols <- ncol(m_list[[1]])
    m_length <- rows * cols
    combined_matrix <- matrix(0, nrow=rows, ncol=cols)
    
    for (i in 1:m_length) {
        elem_vals <- sapply(m_list, function(mat) {mat[[i]]})
        target_val <- comb_func(elem_vals)
        combined_matrix[i] <- target_val
    }
    
    combined_matrix
}







