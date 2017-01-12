#' Perform normalizations on Normalyzer dataset
#' 
#' @param nds Normalyzer dataset object.
#' @param currentjob Name of the ongoing processing run.
#' @param jobdir Path to the output directory
#' @return Returns Normalyzer results object with performed analyzes assigned
#'  as attributes
#' @export
normMethods <- function(nds, currentjob, forceAll=FALSE, normalizeRetentionTime=FALSE, retentionTimeWindow=0.1) {

    nr <- generateNormalyzerResultsObject(nds)
    nr <- performNormalizations(nr, forceAll=forceAll, rtNorm=normalizeRetentionTime, rtWindow=retentionTimeWindow)
    
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

#' Retrieve indices for first or last occurences in vector with replicated 
#' elements
#' 
#' @param targetVector Input vector with replicated elements.
#' @param reverse Look for first or last occurence for each element.
#'  By default looks for the first occurence.
#' @return Vector with indices for each first occurence.
getFirstIndicesInVector <- function(targetVector, reverse=FALSE) {
    
    encounteredNumbers <- c()
    firstIndices <- c()
    
    if (!reverse) {
        startIndex <- 1
        endIndex <- length(targetVector)
    }
    else {
        startIndex <- length(targetVector)
        endIndex <- 1
    }
    
    for (i in startIndex:endIndex) {
        targetValue <- targetVector[i]
        
        if (!is.element(targetValue, encounteredNumbers)) {
            
            encounteredNumbers <- append(encounteredNumbers, targetValue)
            
            if (!reverse) {
                firstIndices <- append(firstIndices, i)
            }
            else {
                firstIndices <- append(i, firstIndices)
            }
        }
    }
    
    firstIndices
}

#' Write normalization matrices to file
#' 
#' @param nr Normalyzer results 
#' @return None
writeNormalizedDatasets <- function(nr, jobdir) {
    
    nds <- nr@nds
    methodnames <- getUsedMethodNames(nr)
    methodlist <- getNormalizationMatrices(nr)
    annotationColumns <- nds@annotationValues
    
    for (sampleIndex in 1:length(methodnames)) {
        
        currentMethod <- methodnames[sampleIndex]
        filePath <- paste(jobdir, "/", currentMethod, "-normalized.txt", sep="")
        outputTable <- cbind(annotationColumns, methodlist[[sampleIndex]])
        utils::write.table(outputTable, file=filePath, sep="\t", row.names=FALSE, col.names=nds@rawData[2,], quote=FALSE)
    }
    
    if (!all(is.na(nr@houseKeepingVars))) {
        hkVarsName <- "housekeeping-variables.tsv"
        hkFilePath <- paste(jobdir, "/", hkVarsName, sep="")
        utils::write.table(file=hkFilePath, nr@houseKeepingVars, sep="\t", 
                           row.names=FALSE, col.names=nds@rawData[2,], quote=FALSE)
    }
    
    rawdata_name <- "submitted_rawdata.tsv"
    rawFilePath <- paste(jobdir, "/", rawdata_name, sep="")
    rawOutputTable <- cbind(annotationColumns, nds@filterrawdata)
    
    utils::write.table(rawOutputTable, file=rawFilePath, sep="\t", row.names=FALSE, col.names=nds@rawData[2,], quote=FALSE)
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

#' Perform desired normalization on time-windows, before summing together
#' The idea is to reduce noise introduced from LC fluctuations
#' 
#' @param rawMatrix Target matrix to be normalized
#' @param retentionTimes Vector of retention times corresponding to rawMatrix
#' @param normMethod The normalization method to apply to the time windows
#' @param stepSizeMinutes Size of windows to be normalized
#' @param offset Whether time window should shifted half step size
#' @return Normalized matrix
getRTNormalizedMatrix <- function(rawMatrix, retentionTimes, normMethod, stepSizeMinutes, offset=FALSE) {
    
    startVal <- min(na.omit(retentionTimes))
    endVal <- max(na.omit(retentionTimes))
    rowNumbers <- c()

    print(paste("Original start time", startVal))
        
    if (offset) {
        startVal <- startVal - stepSizeMinutes / 2
    }
    
    print(paste("New start time", startVal))
    
    processedRows <- matrix(, ncol=ncol(rawMatrix), nrow=0)
    
    for (windowStart in seq(startVal, endVal, stepSizeMinutes)) {
        
        print(paste("Window start: ", windowStart))
        
        windowEnd <- windowStart + stepSizeMinutes
        sliceRows <- which(retentionTimes >= windowStart & retentionTimes < windowEnd)
        rowNumbers <- c(rowNumbers, sliceRows)
        
        currentRows <- rawMatrix[sliceRows,, drop=FALSE]
        
        if (length(currentRows) > 0) {
            
            processedSlice <- normMethod(currentRows)
            processedRows <- rbind(processedRows, processedSlice)
        }
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
getSmoothedRTNormalizedMatrix <- function(rawMatrix, retentionTimes, normMethod, stepSizeMinutes) {
    
    matrixWithoutOffset <- getRTNormalizedMatrix(rawMatrix, retentionTimes, normMethod, stepSizeMinutes, offset=FALSE)
    matrixWithOffset <- getRTNormalizedMatrix(rawMatrix, retentionTimes, normMethod, stepSizeMinutes, offset=TRUE)

    averagedMatrices <- (matrixWithOffset + matrixWithoutOffset) / 2
    averagedMatrices
}


