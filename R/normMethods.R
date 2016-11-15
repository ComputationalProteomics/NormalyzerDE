DEBUG_PRINTS_ON = T

normMethods <- function(nds, currentjob, jobdir) {

    nr <- generateNormalyzerResultsObject(nds)
    nr <- performNormalizations(nr)
    
    methodlist <- getMethodList(!is.null(nr@houseKeepingVars), nds)
    methodnames <- getMethodNames(!is.null(nr@houseKeepingVars))

    filterrawdata <- nr@nds@filterrawdata
    
    # Perform other norm. if the dataset is not small  
    if (nrow(nds@filterrawdata) > 50) {
        
        # nr <- performVSNNormalization(nr)
        # nr <- performQuantileNormalization(nr)
        # nr <- performSMADNormalization(nr)
        # nr <- performCyclicLoessNormalization(nr)
        # nr <- performGlobalRLRNormalization(nr)
        # 
        # ## NORMALIZATION within REPLICATES RLR, VSN and Loess
        # nr <- performReplicateBasedNormalizations(nr)

        # --- This entire part should be automated as part of result object and removed
        if (!is.null(nr@houseKeepingVars)) {
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
        
        print(nds@rawData[2,])
        print(methodlist[[sampleIndex]])
        print(methodnames[sampleIndex])
        
        write.table(file=paste(jobdir, "/", methodnames[sampleIndex], "-normalized.txt", sep=""), 
                    cbind(nds@rawData[-(1:2), (1:(length(nds@inputHeaderValues) - length(nds@sampleReplicateGroups)))], 
                          methodlist[[sampleIndex]]), sep="\t", row.names=F, col.names=nds@rawData[2,], quote=F)
    }
    
    if (!is.null(nr@houseKeepingVars)) {
        write.table(file=paste(jobdir, "/housekeeping-variables.txt", sep=""), nr@houseKeepingVars, sep="\t", row.names=F, col.names=nds@rawData[2,], quote=F)
    }
    
    write.table(file=paste(jobdir, "/submitted_rawdata.txt", sep=""), 
                cbind(nds@rawData[-(1:2), (1:(length(nds@inputHeaderValues) - length(nds@sampleReplicateGroups)))], nds@filterrawdata), sep="\t", row.names=F,
                col.names=nds@rawData[2,], quote=F)
    methodlist <- list(methodlist, methodnames, nds@rawData, nds@filterrawdata, nds@sampleReplicateGroups, !is.null(nr@houseKeepingVars))
    
    return(methodlist)
}

generateNormalyzerResultsObject <- function(nds) {
    nr <- NormalyzerResults(nds=nds)
    nr <- initializeResultsObject(nr)
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

# performVSNNormalization <- function(nr) {
#     
#     filterrawdata <- nr@nds@filterrawdata
#     nr@data2vsn <- justvsn(filterrawdata)
#     nr
# }

# performQuantileNormalization <- function(nr) {
#     
#     stopifnot(!is.null(nr@data2log2))
#     
#     nr@data2quantile <- normalize.quantiles(nr@data2log2, copy=T)
#     nr
# }

# performSMADNormalization <- function(nr) {
# 
#     mediandata <- apply(nr@data2log2, 2, "median", na.rm=T)
#     maddata <- apply(nr@data2log2, 2, function(x) mad(x, na.rm=T))
#     nr@data2mad <- t(apply(nr@data2log2, 1, function(x) ((x - mediandata) / maddata)))
#     nr@data2mad <- nr@data2mad + mean(mediandata)
#     nr
# }

# performCyclicLoessNormalization <- function(nr) {
# 
#     nr@data2loess <- normalizeCyclicLoess(nr@data2log2, method="fast")
#     nr
# }

# performGlobalRLRNormalization <- function(nr) {
# 
#     mediandata <- apply(nr@data2log2, 1, "median", na.rm=T)
#     isFirstSample <- TRUE
# 
#     for (j in 1:ncol(nr@data2log2)) {
# 
#         # print(sprintf("Index: %i", j))
# 
#         lrFit <- rlm(as.matrix(nr@data2log2[, j])~mediandata, na.action=na.exclude)
#         coeffs <- lrFit$coefficients
#         coeffs2 <- coeffs[2]
#         coeffs1 <- coeffs[1]
# 
#         if (isFirstSample) {
#             value_to_assign <- (nr@data2log2[, j] - coeffs1) / coeffs2
#             globalfittedRLR <- (nr@data2log2[, j] - coeffs1) / coeffs2
#             isFirstSample <- FALSE
#         }
#         else {
#             globalfittedRLR <- cbind(globalfittedRLR, (nr@data2log2[, j] - coeffs1) / coeffs2)
#         }
#     }
# 
#     nr@globalfittedRLR <- globalfittedRLR
#     colnames(nr@globalfittedRLR) <- colnames(nr@data2log2)
# 
#     nr
# }

# performReplicateBasedNormalizations <- function(nr) {
#     # performReplicateBasedNormalizations <- function(normObj, sampleReplicateGroups, filterrawdata) {
#         
#     # nds@sampleReplicateGroups, nds@filterrawdata
#         
#     nds <- nr@nds
#     
#     sampleReplicateGroups <- nds@sampleReplicateGroups
#     filterrawdata <- nds@filterrawdata
# 
#     
#     firstIndices <- getFirstIndicesInVector(sampleReplicateGroups, reverse=F)
#     lastIndices <- getFirstIndicesInVector(sampleReplicateGroups, reverse=T)
#     
#     stopifnot(length(firstIndices) == length(lastIndices))
#     
#     nr@fittedLR <- matrix(, nrow=nrow(filterrawdata), ncol=0)
#     nr@data2limloess <- matrix(, nrow=nrow(filterrawdata), ncol=0)
#     nr@data2vsnrep <- matrix(, nrow=nrow(filterrawdata), ncol=0)
#     
#     for (sampleIndex in 1:length(firstIndices)) {
#         
#         startColIndex <- firstIndices[sampleIndex]
#         endColIndex <- lastIndices[sampleIndex]
#         
#         # Median based LR normalization
#         mediandata <- apply(nr@data2log2[, startColIndex:endColIndex], 1, "median", na.rm=T)
#         
#         for (currentCol in startColIndex:endColIndex) {
#             
#             LRfit <- rlm(as.matrix(nr@data2log2[, currentCol])~mediandata, na.action=na.exclude)
#             LRcoeffs <- LRfit$coefficients
#             nr@fittedLR <- cbind(nr@fittedLR, (nr@data2log2[, currentCol] - LRcoeffs[1]) / LRcoeffs[2])
#         }
#         
#         nr@data2limloess <- cbind(nr@data2limloess, normalizeCyclicLoess(nr@data2log2[, startColIndex:endColIndex], method="fast"))
#         nr@data2vsnrep <- cbind(nr@data2vsnrep, justvsn(as.matrix(filterrawdata[, startColIndex:endColIndex])))
#     }
# 
#     colnames(nr@fittedLR) <- colnames(nr@data2log2)
#     colnames(nr@data2limloess) <- colnames(nr@data2log2)
#     colnames(nr@data2vsnrep) <- colnames(nr@data2log2)
#         
#     nr
# }

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


