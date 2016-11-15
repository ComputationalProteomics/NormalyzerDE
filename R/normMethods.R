DEBUG_PRINTS_ON = T

normMethods <- function(normObj, currentjob, jobdir) {
# normMethods <- function(datafile, currentjob, outputDir=NULL) {
        
    print("DEBUG: normMethods entered !!")
    
    # If one of columns are zero here, are errors propagated?
    colsum <- colSums(normObj@filterrawdata, na.rm=T)

    medofdata <- apply(normObj@filterrawdata, 2, FUN="median", na.rm=T)  
    meanofdata <- apply(normObj@filterrawdata, 2, FUN="mean", na.rm=T) 
    
    normObj <- performBasicNormalizations(normObj, normObj@filterrawdata)
    houseKeepingFlag <- (nrow(normObj@normfinderFilterRawData) < 1000)
    # houseKeepingFlag <- (nrow(normObj@filterRawDataNormfinder) < 1000)
    
    if (houseKeepingFlag) {
        Hkvar <- normfinder(normObj@normfinderFilterRawData, normObj@inputHeaderValues)
        # Hkvar <- normfinder(normObj@filterRawDataNormfinder, normObj@inputHeaderValues)
        normObj <- calculateHKdataForNormObj(normObj, normObj@inputHeaderValues, normObj@sampleReplicateGroups, Hkvar, normObj@filterrawdata)
    }
    normObj <- basicMetricNormalizations(normObj, normObj@filterrawdata, colsum, medofdata, meanofdata)
    methodlist <- getMethodList(houseKeepingFlag, normObj)
    methodnames <- getMethodNames(houseKeepingFlag)
    
    # Perform other norm. if the dataset is not small  
    if (nrow(normObj@filterrawdata) > 50) {
        
        normObj <- performVSNNormalization(normObj, normObj@filterrawdata)
        normObj <- performQuantileNormalization(normObj)
        normObj <- performSMADNormalization(normObj)
        normObj <- performCyclicLoessNormalization(normObj)
        normObj <- performGlobalRLRNormalization(normObj)
        
        ## NORMALIZATION within REPLICATES RLR, VSN and Loess
        normObj <- performReplicateBasedNormalizations(normObj, normObj@sampleReplicateGroups, normObj@filterrawdata)
        
        colnames(normObj@fittedLR) <- colnames(normObj@data2log2)
        colnames(normObj@data2quantile) <- colnames(normObj@data2log2)
        
        if (houseKeepingFlag) {
            methodlist <- list(normObj@data2log2, normObj@data2limloess, normObj@fittedLR, normObj@data2vsnrep, normObj@data2loess, 
                               normObj@globalfittedRLR, normObj@data2vsn, normObj@data2GI, normObj@data2med, normObj@data2mean, 
                               normObj@data2ctrlog, normObj@data2quantile)
            methodnames <- c("Log2", "Loess-R", "RLR-R", "VSN-R", "Loess-G", "RLR-G", "VSN-G", "TI-G", "MedI-G", "AI-G", "NF-G", "Quantile")
        } 
        else {
            methodlist <- list(normObj@data2log2, normObj@data2limloess, normObj@fittedLR, normObj@data2vsnrep, normObj@data2loess, 
                               normObj@globalfittedRLR, normObj@data2vsn, normObj@data2GI, normObj@data2med, normObj@data2mean, normObj@data2quantile)
            methodnames <- c("Log2", "Loess-R", "RLR-R", "VSN-R", "Loess-G", "RLR-G", "VSN-G", "TI-G", "MedI-G", "AI-G", "Quantile")
        }
    }
    
    # Create tmp dir for the current job
    for (sampleIndex in 1:length(methodlist)) {
        write.table(file=paste(jobdir, "/", methodnames[sampleIndex], "-normalized.txt", sep=""), 
                    cbind(normObj@rawData[-(1:2), (1:(length(normObj@inputHeaderValues) - length(normObj@sampleReplicateGroups)))], 
                          methodlist[[sampleIndex]]), sep="\t", row.names=F, col.names=normObj@rawData[2,], quote=F)
    }
    
    if (houseKeepingFlag) {
        write.table(file=paste(jobdir, "/housekeeping-variables.txt", sep=""), Hkvar, sep="\t", row.names=F, col.names=normObj@rawData[2,], quote=F)
    }
    
    write.table(file=paste(jobdir, "/submitted_rawdata.txt", sep=""), 
                cbind(normObj@rawData[-(1:2), (1:(length(normObj@inputHeaderValues) - length(normObj@sampleReplicateGroups)))], normObj@filterrawdata), sep="\t", row.names=F,
                col.names=normObj@rawData[2,], quote=F)
    methodlist <- list(methodlist, methodnames, normObj@rawData, normObj@filterrawdata, normObj@sampleReplicateGroups, houseKeepingFlag)
    
    # sink(NULL)
    
    return(methodlist)
}

performBasicNormalizations <- function(normObj, filterrawdata) {
    if (DEBUG_PRINTS_ON) { print("Function: performBasicNormalizations") }
    
    normObj@data2log2 <- log2(filterrawdata)
    normObj@data2GI <- matrix(nrow=nrow(filterrawdata), ncol=ncol(filterrawdata), byrow=T)
    normObj@data2ctr <- matrix(nrow=nrow(filterrawdata), ncol=ncol(filterrawdata), byrow=T)
    normObj@data2med <- matrix(nrow=nrow(filterrawdata), ncol=ncol(filterrawdata), byrow=T) 
    normObj@data2mean <- matrix(nrow=nrow(filterrawdata), ncol=ncol(filterrawdata), byrow=T)
    
    normObj
}

calculateHKdataForNormObj <- function(normObj, getEDdata, filterED, Hkvar, filterrawdata) {
    if (DEBUG_PRINTS_ON) { print("Function: calculateHKdataForNormObj") }
    
    Hkvartemp <- as.matrix(Hkvar[, -(1:(length(getEDdata) - length(filterED)))])
    class(Hkvartemp) <- "numeric"
    colmedianctr <- apply(Hkvartemp, 2, FUN="mean")
    
    for(i in 1:nrow(filterrawdata)) {
        normObj@data2ctr[i,] <- unlist(sapply(1:ncol(filterrawdata), 
                                              function(zd) { (filterrawdata[i, zd] / colmedianctr[zd]) * (mean(colmedianctr)) }))
    }
    normObj@data2ctrlog <- log2(normObj@data2ctr)
    colnames(normObj@data2ctrlog) <- colnames(normObj@data2log2)
    
    normObj
}

## Normalizing using average sample sum, sample mean and sample median
basicMetricNormalizations <- function(normObj, filterrawdata, colsum, medofdata, meanofdata) {
    if (DEBUG_PRINTS_ON) { print("Function: basicMetricNormalizations") }
    
    avgcolsum <- median(colsum)
    for (i in 1:nrow(filterrawdata)) {
        normObj@data2GI[i, ] <- unlist(sapply(1:ncol(filterrawdata), function(zd) { (filterrawdata[i, zd] / colsum[zd]) * (avgcolsum) }))
        normObj@data2med[i, ] <- unlist(sapply(1:ncol(filterrawdata), function(zd) { (filterrawdata[i, zd] / medofdata[zd]) * mean(medofdata) }))
        normObj@data2mean[i, ] <- unlist(sapply(1:ncol(filterrawdata), function(zd) { (filterrawdata[i, zd] / meanofdata[zd]) * mean(meanofdata) }))
    } 

    # print(head(normObj@data2GI))
    # print(head(normObj@data2med))
    # print(head(normObj@data2mean))
    
    normObj@data2GI <- log2(normObj@data2GI)
    normObj@data2med <- log2(normObj@data2med)
    normObj@data2mean <- log2(normObj@data2mean)
    colnames(normObj@data2GI) <- colnames(normObj@data2log2)
    colnames(normObj@data2med) <- colnames(normObj@data2log2)
    colnames(normObj@data2mean) <- colnames(normObj@data2log2)
    
    normObj
}

getMethodList <- function(HKflag, normObj) {
    if (DEBUG_PRINTS_ON) { print("Function: getMethodList") }
    
    if (HKflag) {
        methodlist <- list(normObj@data2log2, normObj@data2GI, normObj@data2med, normObj@data2mean, normObj@data2ctrlog)
    }
    else {
        methodlist <- list(normObj@data2log2, normObj@data2GI, normObj@data2med, normObj@data2mean)
    }
    methodlist
}

getMethodNames <- function(houseKeepingFlag) {
    if (DEBUG_PRINTS_ON) { print("Function: getMethodNames") }
    
    if (houseKeepingFlag) {
        methodnames <- c("Log2", "TI-G", "MedI-G", "AI-G", "NF-G")
    }
    else {
        methodnames <- c("Log2", "TI-G", "MedI-G", "AI-G")
    }
    methodnames
}

performVSNNormalization <- function(normObj, filterrawdata) {
    if (DEBUG_PRINTS_ON) { print("Function: performVSNNormalization") }
    
    # sink(sinkfile, type="output")
    # cat("\n###Warning messages for VSN-G###\n")
    # sink(sinkfile, type="message")
    normObj@data2vsn <- justvsn(filterrawdata)
    
    normObj
}

performQuantileNormalization <- function(normObj) {
    if (DEBUG_PRINTS_ON) { print("Function: performQuantileNormalization") }
    
    normObj@data2quantile <- normalize.quantiles((normObj@data2log2), copy=T)
    normObj
}

performSMADNormalization <- function(normObj) {
    if (DEBUG_PRINTS_ON) { print("Function: performSMADNormalization") }
    
    mediandata <- apply(normObj@data2log2, 2, "median", na.rm=T)
    maddata <- apply(normObj@data2log2, 2, function(x) mad(x, na.rm=T))
    normObj@data2mad <- t(apply(normObj@data2log2, 1, function(x) ((x - mediandata) / maddata)))
    normObj@data2mad <- normObj@data2mad + mean(mediandata)
    normObj
}

performCyclicLoessNormalization <- function(normObj) {
    if (DEBUG_PRINTS_ON) { print("Function: performCyclicLoessNormalization") }
    
    normObj@data2loess <- normalizeCyclicLoess(normObj@data2log2, method="fast")
    normObj
}

performGlobalRLRNormalization <- function(normObj) {
    if (DEBUG_PRINTS_ON) { print("Function: performGlobalLRLNormalization") }
    

    mediandata <- apply(normObj@data2log2, 1, "median", na.rm=T)
    isFirstSample <- TRUE

    for (j in 1:ncol(normObj@data2log2)) {

        # print(sprintf("Index: %i", j))

        LRfit <- rlm(as.matrix(normObj@data2log2[, j])~mediandata, na.action=na.exclude)
        Coeffs <- LRfit$coefficients
        coeffs2 <- Coeffs[2]
        coeffs1 <- Coeffs[1]

        if (isFirstSample) {
            value_to_assign <- (normObj@data2log2[, j] - coeffs1) / coeffs2
            globalfittedRLR <- (normObj@data2log2[, j] - coeffs1) / coeffs2
            isFirstSample <- FALSE
        }
        else {
            globalfittedRLR <- cbind(globalfittedRLR, (normObj@data2log2[, j] - coeffs1) / coeffs2)
        }
    }

    normObj@globalfittedRLR <- globalfittedRLR
    colnames(normObj@globalfittedRLR) <- colnames(normObj@data2log2)

    normObj
}

performReplicateBasedNormalizations <- function(normObj, sampleReplicateGroups, filterrawdata) {
    if (DEBUG_PRINTS_ON) { print("Function: performReplicateBasedNormalizations") }
    
    firstIndices <- getFirstIndicesInVector(sampleReplicateGroups, reverse=F)
    lastIndices <- getFirstIndicesInVector(sampleReplicateGroups, reverse=T)
    stopifnot(length(firstIndices) == length(lastIndices))
    
    normObj@fittedLR <- matrix(, nrow=nrow(filterrawdata), ncol=0)
    normObj@data2limloess <- matrix(, nrow=nrow(filterrawdata), ncol=0)
    normObj@data2vsnrep <- matrix(, nrow=nrow(filterrawdata), ncol=0)
    
    for (sampleIndex in 1:length(firstIndices)) {
        
        startColIndex <- firstIndices[sampleIndex]
        endColIndex <- lastIndices[sampleIndex]
        
        # Median based LR normalization
        mediandata <- apply(normObj@data2log2[, startColIndex:endColIndex], 1, "median", na.rm=T)
        
        for (currentCol in startColIndex:endColIndex) {
            
            LRfit <- rlm(as.matrix(normObj@data2log2[, currentCol])~mediandata, na.action=na.exclude)
            LRcoeffs <- LRfit$coefficients
            normObj@fittedLR <- cbind(normObj@fittedLR, (normObj@data2log2[, currentCol] - LRcoeffs[1]) / LRcoeffs[2])
        }
        
        normObj@data2limloess <- cbind(normObj@data2limloess, normalizeCyclicLoess(normObj@data2log2[, startColIndex:endColIndex], method="fast"))
        normObj@data2vsnrep <- cbind(normObj@data2vsnrep, justvsn(as.matrix(filterrawdata[, startColIndex:endColIndex])))
    }
    
    normObj
}

## Retrieve index vector with first or last occurences of each element
getFirstIndicesInVector <- function(targetVector, reverse=F) {
    
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


