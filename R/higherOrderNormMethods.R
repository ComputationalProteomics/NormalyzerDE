#' Perform desired normalization on time-windows, before summing together
#' The idea is to reduce noise introduced from LC fluctuations
#' 
#' @param rawMatrix Target matrix to be normalized
#' @param retentionTimes Vector of retention times corresponding to rawMatrix
#' @param normMethod The normalization method to apply to the time windows
#' @param stepSizeMinutes Size of windows to be normalized
#' @param offset Whether time window should shifted half step size
#' @return Normalized matrix
getRTNormalizedMatrix <- function(rawMatrix, retentionTimes, normMethod, stepSizeMinutes, windowMinCount=1, offset=0) {
    
    sortedRT <- sort(retentionTimes)
    
    startVal <- min(na.omit(retentionTimes))
    endVal <- max(na.omit(retentionTimes))
    rowNumbers <- c()
    
    if (offset) {
        startVal <- startVal - stepSizeMinutes * offset
    }
    
    processedRows <- matrix(, ncol=ncol(rawMatrix), nrow=0)
    windowStart <- startVal
    windowEnd <- startVal + stepSizeMinutes
    
    while (windowEnd < endVal) {
        
        windowEnd <- windowStart + stepSizeMinutes
        firstTimeSlicedIndices <- which(retentionTimes >= windowStart & retentionTimes < windowEnd)
        
        if (length(firstTimeSlicedIndices) < windowMinCount) {
            
            startPos <- which(sortedRT >= windowStart)[1]
            if (startPos + windowMinCount < endVal) {
                endPos <- startPos + windowMinCount
            }
            else {
                endPos <- length(retentionTimes)
            }
            
            windowStart <- sortedRT[startPos]
            windowEnd <- sortedRT[endPos]
            nextWindowStart <- sortedRT[endPos + 1]
            
            fixSlicedIndices <- c(startPos:endPos)
            posSliced <- sortedRT[fixSlicedIndices]
            timeSlicedIndices <- match(posSliced, retentionTimes)
        }
        else {
            timeSlicedIndices <- firstTimeSlicedIndices
            nextWindowStart <- windowEnd
        }

        rowNumbers <- c(rowNumbers, timeSlicedIndices)
        currentRows <- rawMatrix[timeSlicedIndices,, drop=FALSE]
        
        processedSlice <- normMethod(currentRows)
        processedRows <- rbind(processedRows, processedSlice)
        
        windowStart <- nextWindowStart
    }
    
    orderedProcessedRows <- processedRows[order(rowNumbers), ]
    orderedProcessedRows
}

#' Generate two RT time-window normalized matrices where one is shifted.
#' Then return the mean of these matrices.
#' 
#' @param rawMatrix Target matrix to be normalized
#' @param retentionTimes Vector of retention times corresponding to rawMatrix
#' @param normMethod The normalization method to apply to the time windows
#' @param stepSizeMinutes Size of windows to be normalized
#' @param offset Whether time window should shifted half step size
#' @return Normalized matrix
getSmoothedRTNormalizedMatrix <- function(rawMatrix, retentionTimes, normMethod, stepSizeMinutes, frame_shifts=2, win_size_min=1) {
    
    matrices <- list()

    for (i in 1:frame_shifts) {
        frac_shift <- (i - 1) * 1 / frame_shifts
        matrices[[i]] <- getRTNormalizedMatrix(rawMatrix, retentionTimes, normMethod, stepSizeMinutes, windowMinCount=win_size_min, offset=frac_shift)
    }
    
    # mean, median, anonymous...
    combinedMatrices <- getCombinedMatrix(matrices, median)
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







