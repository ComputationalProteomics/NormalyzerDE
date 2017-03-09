#' Perform normalizations on Normalyzer dataset
#' 
#' @param nds Normalyzer dataset object.
#' @param currentjob Name of the ongoing processing run.
#' @param jobdir Path to the output directory
#' @return Returns Normalyzer results object with performed analyzes assigned
#'  as attributes
#' @export
normMethods <- function(nds, currentjob, forceAll=FALSE, normalizeRetentionTime=FALSE, retentionTimeWindow=0.1, runNormfinder=TRUE) {

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
writeNormalizedDatasets <- function(nr, jobdir, include_pvals=FALSE, include_pairwise_comparisons=FALSE) {
    
    nds <- nr@nds
    ner <- nr@ner
    
    methodnames <- getUsedMethodNames(nr)
    methodlist <- getNormalizationMatrices(nr)
    annotationColumns <- nds@annotationValues
    
    for (sampleIndex in 1:length(methodnames)) {
        
        currentMethod <- methodnames[sampleIndex]
        filePath <- paste(jobdir, "/", currentMethod, "-normalized.txt", sep="")

        outputTable <- cbind(annotationColumns, methodlist[[sampleIndex]])
        out_numbering <- nds@rawData[1,]
        out_colnames <- nds@rawData[2,]
        
        if (include_pvals) {
            # outputTable <- cbind(annotationColumns, methodlist[[sampleIndex]])
            # out_colnames <- nds@rawData[2,]
            
            anova_col <- ner@anfdr[,sampleIndex]
            kw_col <- ner@kwfdr[,sampleIndex]
            
            outputTable <- cbind(outputTable, anova_col, kw_col)
            out_numbering <- c(out_numbering, 0, 0)
            out_colnames <- c(out_colnames, 'anova', 'kruskal_wallis')
        }
        
        if (include_pairwise_comparisons) {
            
            for (comp in names(nr@ner@pairwise_comps)) {
                
                out_numbering <- c(out_numbering, 0)
                out_colnames <- c(out_colnames, paste("comp", comp, sep="_"))
                
                comp_col <- ner@pairwise_comps[[comp]][,sampleIndex]
                
                outputTable <- cbind(outputTable, comp_col)
            }
        }

        outputTable <- rbind(out_numbering, out_colnames, outputTable)
        utils::write.table(outputTable, file=filePath, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
        # utils::write.table(outputTable, file=filePath, sep="\t", row.names=FALSE, col.names=out_colnames, quote=FALSE)
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


