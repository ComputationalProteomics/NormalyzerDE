#' Perform normalizations on Normalyzer dataset
#' 
#' @param nds Normalyzer dataset object.
#' @param forceAll Force all methods to run despite not qualifying for thresholds.
#' @param normalizeRetentionTime Perform retention time based normalization methods.
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
#' data(example_summarized_experiment)
#' normObj <- getVerifiedNormalyzerObject("job_name", example_summarized_experiment)
#' normResults <- normMethods(normObj)
normMethods <- function(nds, forceAll=FALSE, normalizeRetentionTime=TRUE, 
                        quiet=FALSE, rtStepSizeMinutes=1, rtWindowMinCount=100, 
                        rtWindowShifts=1, rtWindowMergeMethod="mean") {
    
    nr <- NormalyzerResults(nds=nds)
    # nr <- initializeResultsObject(nr)
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

#' The normalization divides the intensity of each variable in a sample
#' with the sum of intensities of all variables in the sample and multiplies
#' with the median of sum of intensities of all variables in all samples.
#' The normalized data is then log2-transformed.
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

#' Intensity of each variable in a given sample is divided by the median of
#' intensities of all variables in the sample and then multiplied by the mean
#' of median of sum of intensities of all variables in all samples.
#' The normalized data is then log2-transformed.
#' 
#' @param rawMatrix Target matrix to be normalized
#' @return Normalized and log-transformed matrix
#' @export
#' @examples
#' data(example_data_only_values)
#' normMatrix <- medianNormalization(example_data_only_values)
medianNormalization <- function(rawMatrix) {

    colMedians <- matrixStats::colMedians(rawMatrix, na.rm=TRUE)
    meanColMedian <- mean(colMedians, na.rm=TRUE)
    normMatrix <- matrix(nrow=nrow(rawMatrix), ncol=ncol(rawMatrix), byrow=TRUE)
    normFunc <- function(colIndex) { 
        (rawMatrix[rowIndex, colIndex] / colMedians[colIndex]) * meanColMedian 
    }
    
    for (rowIndex in seq_len(nrow(rawMatrix))) {
        normMatrix[rowIndex, ] <- vapply(seq_len(ncol(rawMatrix)), normFunc, 0)
    }
    
    normLog2Matrix <- log2(normMatrix)
    colnames(normLog2Matrix) <- colnames(rawMatrix)
    normLog2Matrix
}

#' Intensity of each variable in a given sample is divided by the mean of sum of
#' intensities of all variables in the sample and then multiplied by the mean
#' of sum of intensities of all variables in all samples. The normalized data
#' is then transformed to log2.
#' 
#' @param rawMatrix Target matrix to be normalized
#' @return Normalized and log-transformed matrix
#' @export
#' @examples
#' data(example_data_only_values)
#' normMatrix <- meanNormalization(example_data_only_values)
meanNormalization <- function(rawMatrix) {
    
    colMeans <- colMeans(rawMatrix, na.rm=TRUE)
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

#' Log2 transformed data is normalized using the function "justvsn" from the
#' VSN package. 
#' 
#' The VSN (Variance Stabilizing Normalization) attempts to transform the data 
#' in such a way that the variance remains nearly constant over the intensity 
#' spectrum
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

#' Quantile normalization is performed by the function "normalize.quantiles"
#' from the package preprocessCore.
#' 
#' It makes the assumption that the data in different samples should originate
#' from an identical distribution. It does this by generating a reference
#' distribution and then scaling the other samples accordingly.
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
#' Normalization subtracts the median and divides the data by the
#' median absolute deviation (MAD).
#' 
#' @param rawMatrix Target matrix to be normalized
#' @return Normalized matrix
#' @export
#' @examples
#' data(example_data_only_values)
#' normMatrix <- performSMADNormalization(example_data_only_values)
performSMADNormalization <- function(rawMatrix) {
    
    log2Matrix <- log2(rawMatrix)
    sampleLog2Median <- matrixStats::colMedians(log2Matrix, na.rm=TRUE)
    sampleMAD <- matrixStats::colMads(log2Matrix, na.rm=TRUE)
    madMatrix <- t(apply(log2Matrix, 1, function(row) ((row - sampleLog2Median) / sampleMAD)))
    
    madPlusMedianMatrix <- madMatrix + mean(sampleLog2Median)
    colnames(madPlusMedianMatrix) <- colnames(rawMatrix)
    
    madPlusMedianMatrix
}

#' Cyclic Loess normalization
#' 
#' Log2 transformed data is normalized by Loess method using the function
#' "normalizeCyclicLoess". Further information is available for the function
#' "normalizeCyclicLoess" in the Limma package.
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
#' Log2 transformed data is normalized by robust linear regression using
#' the function "rlm" from the MASS package. 
#' 
#' @param rawMatrix Target matrix to be normalized
#' @return Normalized matrix
#' @export
#' @examples
#' data(example_data_only_values)
#' normMatrix <- performGlobalRLRNormalization(example_data_only_values)
performGlobalRLRNormalization <- function(rawMatrix) {
    
    log2Matrix <- log2(rawMatrix)
    sampleLog2Median <- matrixStats::rowMedians(log2Matrix, na.rm=TRUE)
    
    calculateRLMForCol <- function(colIndex, sampleLog2Median, log2Matrix) {
        
        lrFit <- MASS::rlm(as.matrix(log2Matrix[, colIndex])~sampleLog2Median, na.action=stats::na.exclude)
        coeffs <- lrFit$coefficients
        coefIntercept <- coeffs[1]
        coefSlope <- coeffs[2]
        globalFittedRLRCol <- (log2Matrix[, colIndex] - coefIntercept) / coefSlope
        globalFittedRLRCol
    }
    
    globalFittedRLR <- vapply(
        seq_len(ncol(log2Matrix)),
        calculateRLMForCol,
        rep(0, nrow(log2Matrix)),
        sampleLog2Median=sampleLog2Median,
        log2Matrix=log2Matrix
    )
    
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
