#' Perform normalizations on Normalyzer dataset
#' 
#' @param nds Normalyzer dataset object.
#' @param currentjob Name of the ongoing processing run.
#' @param forceAll Force all methods to run despite not qualifying for thresholds.
#' @param normalizeRetentionTime Perform retention time based normalization methods.
#' @param retentionTimeWindow Default window size for retention times.
#' @param runNormfinder Run the Normfinder normalization method.
#' @return Returns Normalyzer results object with performed analyzes assigned
#'  as attributes
normMethods <- function(nds, currentjob, forceAll=FALSE, normalizeRetentionTime=FALSE, retentionTimeWindow=1, runNormfinder=TRUE) {
    
    nr <- generateNormalyzerResultsObject(nds)
    nr <- performNormalizations(nr, forceAll=forceAll, rtNorm=normalizeRetentionTime, rtWindow=retentionTimeWindow, runNormfinder=runNormfinder)
    
    return(nr)
}

#' Create empty Normalyzer results object from Normalyzer data object
#' 
#' @param nds Normalyzer dataset object.
#' @return Empty Normalyzer results object.
generateNormalyzerResultsObject <- function(nds) {
    nr <- NormalyzerResults(nds=nds)
    nr <- initializeResultsObject(nr)
    nr
}



#' Normalization adjusting based on total sample intensity
#' 
#' @param rawMatrix Target matrix to be normalized
#' @return Normalized and log-transformed matrix
globalIntensityNormalization <- function(rawMatrix) {
    
    colsums <- colSums(rawMatrix, na.rm=TRUE)
    avgcolsum <- stats::median(colsums)
    normMatrix <- matrix(nrow=nrow(rawMatrix), ncol=ncol(rawMatrix), byrow=TRUE)
    
    normFunc <- function(zd) { (rawMatrix[i, zd] / colsums[zd]) * avgcolsum }
    
    for (i in 1:nrow(rawMatrix)) {
        normMatrix[i, ] <- unlist(sapply(1:ncol(rawMatrix), normFunc))
    }
    
    normLog2Matrix <- log2(normMatrix)
    colnames(normLog2Matrix) <- colnames(rawMatrix)
    normLog2Matrix
}

#' Normalization adjusting towards sample median intensity
#' 
#' @param rawMatrix Target matrix to be normalized
#' @return Normalized and log-transformed matrix
medianNormalization <- function(rawMatrix) {
    
    colMedians <- apply(rawMatrix, 2, FUN="median", na.rm=TRUE)
    avgColMedian <- mean(colMedians, na.rm=TRUE)
    
    normMatrix <- matrix(nrow=nrow(rawMatrix), ncol=ncol(rawMatrix), byrow=TRUE)
    
    normFunc <- function(zd) { (rawMatrix[i, zd] / colMedians[zd]) * avgColMedian }
    
    for (i in 1:nrow(rawMatrix)) {
        normMatrix[i, ] <- unlist(sapply(1:ncol(rawMatrix), normFunc))
    }
    
    normLog2Matrix <- log2(normMatrix)
    colnames(normLog2Matrix) <- colnames(rawMatrix)
    normLog2Matrix
}

#' Normalization adjusting towards sample mean intensity
#' 
#' @param rawMatrix Target matrix to be normalized
#' @return Normalized and log-transformed matrix
meanNormalization <- function(rawMatrix) {
    
    colMeans <-apply(rawMatrix, 2, FUN="mean", na.rm=TRUE)
    avgColMean <- mean(colMeans, na.rm=TRUE)
    normMatrix <- matrix(nrow=nrow(rawMatrix), ncol=ncol(rawMatrix), byrow=TRUE)
    
    normFunc <- function(zd) { (rawMatrix[i, zd] / colMeans[zd]) * avgColMean }
    
    for (i in 1:nrow(rawMatrix)) {
        normMatrix[i, ] <- unlist(sapply(1:ncol(rawMatrix), normFunc))
    }
    
    normLog2Matrix <- log2(normMatrix)
    colnames(normLog2Matrix) <- colnames(rawMatrix)
    normLog2Matrix
}

#' Perform Variance Stabilizing Normalization
#' 
#' @param rawMatrix Target matrix to be normalized
#' @return Normalized matrix
performVSNNormalization <- function(rawMatrix) {
    normMatrix <- vsn::justvsn(rawMatrix)
    colnames(normMatrix) <- colnames(rawMatrix)
    normMatrix
}

#' Perform quantile normalization, where overall distribution is assumed
#' to be the same, and using a 'average' distribution as measure to normalize
#' again
#' 
#' @param rawMatrix Target matrix to be normalized
#' @return Normalized matrix
performQuantileNormalization <- function(rawMatrix) {
    
    log2Matrix <- log2(rawMatrix)
    normMatrix <- preprocessCore::normalize.quantiles(log2Matrix, copy=TRUE)
    
    colnames(normMatrix) <- colnames(rawMatrix)
    normMatrix
}

#' Median absolute deviation normalization
#' Scales values with MAD and adds it with logged sample median
#' TODO: Dig into this one a bit more
#' 
#' @param rawMatrix Target matrix to be normalized
#' @return Normalized matrix
performSMADNormalization <- function(rawMatrix) {
    
    log2Matrix <- log2(rawMatrix)
    sampleLog2Median <- apply(log2Matrix, 2, "median", na.rm=TRUE)
    sampleMAD <- apply(log2Matrix, 2, function(x) stats::mad(x, na.rm=TRUE))
    madMatrix <- t(apply(log2Matrix, 1, function(x) ((x - sampleLog2Median) / sampleMAD)))
    
    madPlusMedianMatrix <- madMatrix + mean(sampleLog2Median)
    colnames(madPlusMedianMatrix) <- colnames(rawMatrix)
    
    madPlusMedianMatrix
}

#' Cyclic Loess normalization
#' Local regression based on k-nearest-neighbor model
#' 
#' @param rawMatrix Target matrix to be normalized
#' @return Normalized matrix
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
performGlobalRLRNormalization <- function(rawMatrix) {
    
    log2Matrix <- log2(rawMatrix)
    sampleLog2Median <- apply(log2Matrix, 1, "median", na.rm=TRUE)
    isFirstSample <- TRUE
    
    for (j in 1:ncol(log2Matrix)) {
        
        lrFit <- MASS::rlm(as.matrix(log2Matrix[, j])~sampleLog2Median, na.action=stats::na.exclude)
        coeffs <- lrFit$coefficients
        coeffs2 <- coeffs[2]
        coeffs1 <- coeffs[1]
        
        if (isFirstSample) {
            value_to_assign <- (log2Matrix[, j] - coeffs1) / coeffs2
            globalFittedRLR <- (log2Matrix[, j] - coeffs1) / coeffs2
            isFirstSample <- FALSE
        }
        else {
            globalFittedRLR <- cbind(globalFittedRLR, (log2Matrix[, j] - coeffs1) / coeffs2)
        }
    }
    
    colnames(globalFittedRLR) <- colnames(rawMatrix)
    
    globalFittedRLR
}

#' Do no normalization (For debugging purposes)
#' 
#' @param rawMatrix Target matrix to be normalized
#' @return Normalized matrix
performNoNormalization <- function(rawMatrix) {
    rawMatrix
}


