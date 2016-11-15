DEBUG_PRINTS_ON = T

normMethods <- function(nds, currentjob, jobdir) {

    
    
    
    
    
    
    
    
    
    
    
    nr <- NormalyzerResults(nds=nds)

    nr <- performBasicNormalizations(nr)
    
    houseKeepingFlag <- nrow(nds@normfinderFilterRawData) < 1000

    if (houseKeepingFlag) {
        Hkvar <- normfinder(nds)
        nr <- calculateHKdataForNormObj(nr, Hkvar)
        # nds <- calculateHKdataForNormObj(nds, nds@inputHeaderValues, nds@sampleReplicateGroups, Hkvar, nds@filterrawdata)
    }
    
    nr <- basicMetricNormalizations(nr)
    
    methodlist <- getMethodList(houseKeepingFlag, nds)
    methodnames <- getMethodNames(houseKeepingFlag)

    print(nds)
    
    filterrawdata <- nr@nds@filterrawdata
    
    
    
    
    
    
    
    
    # Perform other norm. if the dataset is not small  
    if (nrow(nds@filterrawdata) > 50) {
        
        nr <- performVSNNormalization(nr)
        
        nr <- performQuantileNormalization(nr)
        nr <- performSMADNormalization(nr)
        nr <- performCyclicLoessNormalization(nr)
        nr <- performGlobalRLRNormalization(nr)
        
        ## NORMALIZATION within REPLICATES RLR, VSN and Loess
        nr <- performReplicateBasedNormalizations(nr)

        # --- This entire part should be automated as part of result object and removed
        if (houseKeepingFlag) {
            methodlist <- list(nr@data2log2, nr@data2limloess, nr@fittedLR, nr@data2vsnrep, nr@data2loess, 
                               nr@globalfittedRLR, nr@data2vsn, nr@data2GI, nr@data2med, nr@data2mean, 
                               nr@data2ctrlog, nr@data2quantile)
            methodnames <- c("Log2", "Loess-R", "RLR-R", "VSN-R", "Loess-G", "RLR-G", "VSN-G", "TI-G", "MedI-G", "AI-G", "NF-G", "Quantile")
        } 
        else {
            methodlist <- list(nr@data2log2, nr@data2limloess, nr@fittedLR, nr@data2vsnrep, nr@data2loess, 
                               nr@globalfittedRLR, nr@data2vsn, nr@data2GI, nr@data2med, nr@data2mean, nr@data2quantile)
            methodnames <- c("Log2", "Loess-R", "RLR-R", "VSN-R", "Loess-G", "RLR-G", "VSN-G", "TI-G", "MedI-G", "AI-G", "Quantile")
        }
        # -----------------
        
    }

    # Create tmp dir for the current job
    for (sampleIndex in 1:length(methodlist)) {
        write.table(file=paste(jobdir, "/", methodnames[sampleIndex], "-normalized.txt", sep=""), 
                    cbind(nds@rawData[-(1:2), (1:(length(nds@inputHeaderValues) - length(nds@sampleReplicateGroups)))], 
                          methodlist[[sampleIndex]]), sep="\t", row.names=F, col.names=nds@rawData[2,], quote=F)
    }
    
    if (houseKeepingFlag) {
        write.table(file=paste(jobdir, "/housekeeping-variables.txt", sep=""), Hkvar, sep="\t", row.names=F, col.names=nds@rawData[2,], quote=F)
    }
    
    write.table(file=paste(jobdir, "/submitted_rawdata.txt", sep=""), 
                cbind(nds@rawData[-(1:2), (1:(length(nds@inputHeaderValues) - length(nds@sampleReplicateGroups)))], nds@filterrawdata), sep="\t", row.names=F,
                col.names=nds@rawData[2,], quote=F)
    methodlist <- list(methodlist, methodnames, nds@rawData, nds@filterrawdata, nds@sampleReplicateGroups, houseKeepingFlag)
    
    # sink(NULL)

    return(methodlist)
}

# performBasicNormalizations <- function(normObj, filterrawdata) {
#     if (DEBUG_PRINTS_ON) { print("Function: performBasicNormalizations") }
#     
#     normObj@data2log2 <- log2(filterrawdata)
#     normObj@data2GI <- matrix(nrow=nrow(filterrawdata), ncol=ncol(filterrawdata), byrow=T)
#     normObj@data2ctr <- matrix(nrow=nrow(filterrawdata), ncol=ncol(filterrawdata), byrow=T)
#     normObj@data2med <- matrix(nrow=nrow(filterrawdata), ncol=ncol(filterrawdata), byrow=T) 
#     normObj@data2mean <- matrix(nrow=nrow(filterrawdata), ncol=ncol(filterrawdata), byrow=T)
#     
#     normObj
# }

calculateHKdataForNormObj <- function(nr, Hkvar) {
    # calculateHKdataForNormObj <- function(normObj, getEDdata, filterED, Hkvar, filterrawdata) {
        
    if (DEBUG_PRINTS_ON) { print("Function: calculateHKdataForNormObj") }
    
    # nds <- calculateHKdataForNormObj(nds, nds@inputHeaderValues, nds@sampleReplicateGroups, Hkvar, nds@filterrawdata)
    
    nds <- nr@nds
    
    getEDdata <- nds@inputHeaderValues
    filterED <- nds@sampleReplicateGroups
    filterrawdata <- nds@filterrawdata
    
    Hkvartemp <- as.matrix(Hkvar[, -(1:(length(getEDdata) - length(filterED)))])
    class(Hkvartemp) <- "numeric"
    colmedianctr <- apply(Hkvartemp, 2, FUN="mean")
    
    for(i in 1:nrow(filterrawdata)) {
        nr@data2ctr[i,] <- unlist(sapply(1:ncol(filterrawdata), 
                                              function(zd) { (filterrawdata[i, zd] / colmedianctr[zd]) * (mean(colmedianctr)) }))
    }
    
    nr@data2ctrlog <- log2(nr@data2ctr)
    colnames(nr@data2ctrlog) <- colnames(nds@filterrawdata)
    # colnames(normObj@data2ctrlog) <- colnames(normObj@data2log2)
    
    nr
}

## Normalizing using average sample sum, sample mean and sample median
basicMetricNormalizations <- function(nr) {

    nds <- nr@nds
    
    avgcolsum <- median(nds@colsum)
    for (i in 1:nrow(nds@filterrawdata)) {
        nr@data2GI[i, ] <- unlist(sapply(1:ncol(nds@filterrawdata), function(zd) { (nds@filterrawdata[i, zd] / nds@colsum[zd]) * (avgcolsum) }))
        nr@data2med[i, ] <- unlist(sapply(1:ncol(nds@filterrawdata), function(zd) { (nds@filterrawdata[i, zd] / nds@medofdata[zd]) * mean(nds@medofdata) }))
        nr@data2mean[i, ] <- unlist(sapply(1:ncol(nds@filterrawdata), function(zd) { (nds@filterrawdata[i, zd] / nds@meanofdata[zd]) * mean(nds@meanofdata) }))
    } 

    # print(head(normObj@data2GI))
    # print(head(normObj@data2med))
    # print(head(normObj@data2mean))
    
    nr@data2GI <- log2(nr@data2GI)
    nr@data2med <- log2(nr@data2med)
    nr@data2mean <- log2(nr@data2mean)
    colnames(nr@data2GI) <- colnames(nr@data2log2)
    colnames(nr@data2med) <- colnames(nr@data2log2)
    colnames(nr@data2mean) <- colnames(nr@data2log2)
    
    nr
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

performVSNNormalization <- function(nr) {
    
    filterrawdata <- nr@nds@filterrawdata
    nr@data2vsn <- justvsn(filterrawdata)
    nr
}

performQuantileNormalization <- function(nr) {
    
    stopifnot(!is.null(nr@data2log2))
    
    nr@data2quantile <- normalize.quantiles(nr@data2log2, copy=T)
    nr
}

performSMADNormalization <- function(nr) {

    mediandata <- apply(nr@data2log2, 2, "median", na.rm=T)
    maddata <- apply(nr@data2log2, 2, function(x) mad(x, na.rm=T))
    nr@data2mad <- t(apply(nr@data2log2, 1, function(x) ((x - mediandata) / maddata)))
    nr@data2mad <- nr@data2mad + mean(mediandata)
    nr
}

performCyclicLoessNormalization <- function(nr) {

    nr@data2loess <- normalizeCyclicLoess(nr@data2log2, method="fast")
    nr
}

performGlobalRLRNormalization <- function(nr) {

    mediandata <- apply(nr@data2log2, 1, "median", na.rm=T)
    isFirstSample <- TRUE

    for (j in 1:ncol(nr@data2log2)) {

        # print(sprintf("Index: %i", j))

        lrFit <- rlm(as.matrix(nr@data2log2[, j])~mediandata, na.action=na.exclude)
        coeffs <- lrFit$coefficients
        coeffs2 <- coeffs[2]
        coeffs1 <- coeffs[1]

        if (isFirstSample) {
            value_to_assign <- (nr@data2log2[, j] - coeffs1) / coeffs2
            globalfittedRLR <- (nr@data2log2[, j] - coeffs1) / coeffs2
            isFirstSample <- FALSE
        }
        else {
            globalfittedRLR <- cbind(globalfittedRLR, (nr@data2log2[, j] - coeffs1) / coeffs2)
        }
    }

    nr@globalfittedRLR <- globalfittedRLR
    colnames(nr@globalfittedRLR) <- colnames(nr@data2log2)

    nr
}

performReplicateBasedNormalizations <- function(nr) {
    # performReplicateBasedNormalizations <- function(normObj, sampleReplicateGroups, filterrawdata) {
        
    # nds@sampleReplicateGroups, nds@filterrawdata
        
    nds <- nr@nds
    
    sampleReplicateGroups <- nds@sampleReplicateGroups
    filterrawdata <- nds@filterrawdata

    
    firstIndices <- getFirstIndicesInVector(sampleReplicateGroups, reverse=F)
    lastIndices <- getFirstIndicesInVector(sampleReplicateGroups, reverse=T)
    
    stopifnot(length(firstIndices) == length(lastIndices))
    
    nr@fittedLR <- matrix(, nrow=nrow(filterrawdata), ncol=0)
    nr@data2limloess <- matrix(, nrow=nrow(filterrawdata), ncol=0)
    nr@data2vsnrep <- matrix(, nrow=nrow(filterrawdata), ncol=0)
    
    for (sampleIndex in 1:length(firstIndices)) {
        
        startColIndex <- firstIndices[sampleIndex]
        endColIndex <- lastIndices[sampleIndex]
        
        # Median based LR normalization
        mediandata <- apply(nr@data2log2[, startColIndex:endColIndex], 1, "median", na.rm=T)
        
        for (currentCol in startColIndex:endColIndex) {
            
            LRfit <- rlm(as.matrix(nr@data2log2[, currentCol])~mediandata, na.action=na.exclude)
            LRcoeffs <- LRfit$coefficients
            nr@fittedLR <- cbind(nr@fittedLR, (nr@data2log2[, currentCol] - LRcoeffs[1]) / LRcoeffs[2])
        }
        
        nr@data2limloess <- cbind(nr@data2limloess, normalizeCyclicLoess(nr@data2log2[, startColIndex:endColIndex], method="fast"))
        nr@data2vsnrep <- cbind(nr@data2vsnrep, justvsn(as.matrix(filterrawdata[, startColIndex:endColIndex])))
    }

    colnames(nr@fittedLR) <- colnames(nr@data2log2)
    colnames(nr@data2limloess) <- colnames(nr@data2log2)
    colnames(nr@data2vsnrep) <- colnames(nr@data2log2)
        
    nr
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


