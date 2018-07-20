#' Perform normalizations on Normalyzer dataset
#' 
#' @param nds Normalyzer dataset object.
#' @param forceAll Force all methods to run despite not qualifying for thresholds.
#' @param normalizeRetentionTime Perform retention time based normalization methods.
#' @param retentionTimeWindow Default window size for retention times.
#' @param quiet Prevent diagnostic output
#' 
#' @param rtStepSizeMinutes Retention time normalization window size.
#' @param rtWindowMinCount Minimum number of datapoints in each retention-time
#'   segment.
#' @param rtWindowShifts Number of layered retention time normalized windows.
#' @param rtWindowMergeMethod Merge approach for layered retention time windows.
#' 
#' @return Returns Normalyzer results object with performed analyzes assigned
#'  as attributes
#' @export
#' @examples
#' data(example_data)
#' data(example_design)
#' normObj <- getVerifiedNormalyzerObject("job_name", example_design, example_data)
#' normResults <- normMethods(normObj)
normMethods <- function(nds, forceAll=FALSE, normalizeRetentionTime=TRUE, 
                        quiet=FALSE, rtStepSizeMinutes=1, rtWindowMinCount=100, 
                        rtWindowShifts=1, rtWindowMergeMethod="mean") {
    
    nr <- NormalyzerResults(nds=nds)
    nr <- initializeResultsObject(nr)
    nr <- performNormalizations(
        nr, 
        forceAll=forceAll, 
        rtNorm=normalizeRetentionTime, 
        rtStepSizeMinutes=rtStepSizeMinutes, 
        rtWindowMinCount=rtWindowMinCount, 
        rtWindowShifts=rtWindowShifts, 
        rtWindowMergeMethod=rtWindowMergeMethod, 
        quiet=quiet)
    
    return(nr)
}

#' Normalization adjusting based on total sample intensity
#' 
#' @param rawMatrix Target matrix to be normalized
#' @return Normalized and log-transformed matrix
#' @export
#' @examples
#' data(example_data_only_values)
#' normMatrix <- globalIntensityNormalization(example_data_only_values)
globalIntensityNormalization <- function(rawMatrix) {
    
    colSums <- colSums(rawMatrix, na.rm=TRUE)
    colSumsMedian <- stats::median(colSums)
    normMatrix <- matrix(nrow=nrow(rawMatrix), ncol=ncol(rawMatrix), byrow=TRUE)
    
    normFunc <- function(colIndex) { (rawMatrix[rowIndex, colIndex] / colSums[colIndex]) * colSumsMedian }
    
    for (rowIndex in seq_len(nrow(rawMatrix))) {
        normMatrix[rowIndex, ] <- vapply(seq_len(ncol(rawMatrix)), normFunc, 0)
    }
    
    normLog2Matrix <- log2(normMatrix)
    colnames(normLog2Matrix) <- colnames(rawMatrix)
    normLog2Matrix
}

#' Normalization adjusting towards sample median intensity
#' 
#' @param rawMatrix Target matrix to be normalized
#' @return Normalized and log-transformed matrix
#' @export
#' @examples
#' data(example_data_only_values)
#' normMatrix <- medianNormalization(example_data_only_values)
medianNormalization <- function(rawMatrix) {
    
    colMedians <- apply(rawMatrix, 2, FUN="median", na.rm=TRUE)
    meanColMedian <- mean(colMedians, na.rm=TRUE)
    normMatrix <- matrix(nrow=nrow(rawMatrix), ncol=ncol(rawMatrix), byrow=TRUE)
    normFunc <- function(colIndex) { (rawMatrix[rowIndex, colIndex] / colMedians[colIndex]) * meanColMedian }
    
    for (rowIndex in seq_len(nrow(rawMatrix))) {
        normMatrix[rowIndex, ] <- vapply(seq_len(ncol(rawMatrix)), normFunc, 0)
    }
    
    normLog2Matrix <- log2(normMatrix)
    colnames(normLog2Matrix) <- colnames(rawMatrix)
    normLog2Matrix
}

#' Normalization adjusting towards sample mean intensity
#' 
#' @param rawMatrix Target matrix to be normalized
#' @return Normalized and log-transformed matrix
#' @export
#' @examples
#' data(example_data_only_values)
#' normMatrix <- meanNormalization(example_data_only_values)
meanNormalization <- function(rawMatrix) {
    
    colMeans <- apply(rawMatrix, 2, FUN="mean", na.rm=TRUE)
    avgColMean <- mean(colMeans, na.rm=TRUE)
    normMatrix <- matrix(nrow=nrow(rawMatrix), ncol=ncol(rawMatrix), byrow=TRUE)
    
    normFunc <- function(colIndex) { (rawMatrix[rowIndex, colIndex] / colMeans[colIndex]) * avgColMean }
    
    for (rowIndex in seq_len(nrow(rawMatrix))) {
        normMatrix[rowIndex, ] <- vapply(seq_len(ncol(rawMatrix)), normFunc, 0)
    }
    
    normLog2Matrix <- log2(normMatrix)
    colnames(normLog2Matrix) <- colnames(rawMatrix)
    normLog2Matrix
}

#' Perform Variance Stabilizing Normalization
#' 
#' @param rawMatrix Target matrix to be normalized
#' @return Normalized matrix
#' @export
#' @examples
#' data(example_data_only_values)
#' normMatrix <- performVSNNormalization(example_data_only_values)
performVSNNormalization <- function(rawMatrix) {
    
    normMatrix <- suppressMessages(vsn::justvsn(rawMatrix))
    colnames(normMatrix) <- colnames(rawMatrix)
    normMatrix
}

#' Perform quantile normalization, where overall distribution is assumed
#' to be the same, and using a 'average' distribution as measure to normalize
#' again
#' 
#' @param rawMatrix Target matrix to be normalized
#' @return Normalized matrix
#' @export
#' @examples
#' data(example_data_only_values)
#' normMatrix <- performQuantileNormalization(example_data_only_values)
performQuantileNormalization <- function(rawMatrix) {
    
    log2Matrix <- log2(rawMatrix)
    normMatrix <- preprocessCore::normalize.quantiles(log2Matrix, copy=TRUE)
    
    colnames(normMatrix) <- colnames(rawMatrix)
    normMatrix
}

#' Median absolute deviation normalization
#' Scales values with MAD and adds it with logged sample median
#' 
#' @param rawMatrix Target matrix to be normalized
#' @return Normalized matrix
#' @export
#' @examples
#' data(example_data_only_values)
#' normMatrix <- performSMADNormalization(example_data_only_values)
performSMADNormalization <- function(rawMatrix) {
    
    log2Matrix <- log2(rawMatrix)
    sampleLog2Median <- apply(log2Matrix, 2, "median", na.rm=TRUE)
    sampleMAD <- apply(log2Matrix, 2, function(col) stats::mad(col, na.rm=TRUE))
    madMatrix <- t(apply(log2Matrix, 1, function(row) ((row - sampleLog2Median) / sampleMAD)))
    
    madPlusMedianMatrix <- madMatrix + mean(sampleLog2Median)
    colnames(madPlusMedianMatrix) <- colnames(rawMatrix)
    
    madPlusMedianMatrix
}

#' Cyclic Loess normalization
#' Local regression based on k-nearest-neighbor model
#' 
#' @param rawMatrix Target matrix to be normalized
#' @return Normalized matrix
#' @export
#' @examples
#' data(example_data_only_values)
#' normMatrix <- performCyclicLoessNormalization(example_data_only_values)
performCyclicLoessNormalization <- function(rawMatrix) {
    
    log2Matrix <- log2(rawMatrix)
    normMatrix <- limma::normalizeCyclicLoess(log2Matrix, method="fast")
    colnames(normMatrix) <- colnames(rawMatrix)
    
    normMatrix
}

#' Global linear regression normalization
#' 
#' @param rawMatrix Target matrix to be normalized
#' @return Normalized matrix
#' @export
#' @examples
#' data(example_data_only_values)
#' normMatrix <- performGlobalRLRNormalization(example_data_only_values)
performGlobalRLRNormalization <- function(rawMatrix) {
    
    log2Matrix <- log2(rawMatrix)
    sampleLog2Median <- apply(log2Matrix, 1, "median", na.rm=TRUE)
    isFirstSample <- TRUE
    
    for (colIndex in seq_len(ncol(log2Matrix))) {
        
        lrFit <- MASS::rlm(as.matrix(log2Matrix[, colIndex])~sampleLog2Median, na.action=stats::na.exclude)
        coeffs <- lrFit$coefficients
        coefIntercept <- coeffs[1]
        coefSlope <- coeffs[2]
        
        if (isFirstSample) {
            globalFittedRLR <- (log2Matrix[, colIndex] - coefIntercept) / coefSlope
            isFirstSample <- FALSE
        }
        else {
            globalFittedRLR <- cbind(globalFittedRLR, (log2Matrix[, colIndex] - coefIntercept) / coefSlope)
        }
    }
    
    colnames(globalFittedRLR) <- colnames(rawMatrix)
    globalFittedRLR
}

#' Do no normalization (For debugging purposes)
#' 
#' @param rawMatrix Target matrix to be normalized
#' @return Normalized matrix
#' @keywords internal
performNoNormalization <- function(rawMatrix) {
    rawMatrix
}
