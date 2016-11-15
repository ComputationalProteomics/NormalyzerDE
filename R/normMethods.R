DEBUG_PRINTS_ON = T

normMethods <- function(nds, currentjob, jobdir) {

    nr <- generateNormalyzerResultsObject(nds)
    nr <- performNormalizations(nr)
    
    
    
    # methodlist <- getMethodList(!is.null(nr@houseKeepingVars), nds)
    # methodnames <- getMethodNames(!is.null(nr@houseKeepingVars))

    methodnames <- getMethodNames(nr)
    slotNameList <- getSlotNameList(nr)
    
    # print(methodnames)
    # print(methodlist)
    # 
    # stop("")
    
    # Perform other norm. if the dataset is not small  
    # if (nrow(nds@filterrawdata) > 50) {
    #     
    #     # --- This entire part should be automated as part of result object and removed
    #     if (!is.null(nr@houseKeepingVars)) {
    #         methodlist <- list(nr@data2log2, nr@data2limloess, nr@fittedLR, nr@data2vsnrep, nr@data2loess, 
    #                            nr@globalfittedRLR, nr@data2vsn, nr@data2GI, nr@data2med, nr@data2mean, 
    #                            nr@data2ctrlog, nr@data2quantile)
    #         methodnames <- c("Log2", "Loess-R", "RLR-R", "VSN-R", "Loess-G", "RLR-G", "VSN-G", "TI-G", "MedI-G", "AI-G", "NF-G", "Quantile")
    #     } 
    #     else {
    #         methodlist <- list(nr@data2log2, nr@data2limloess, nr@fittedLR, nr@data2vsnrep, nr@data2loess, 
    #                            nr@globalfittedRLR, nr@data2vsn, nr@data2GI, nr@data2med, nr@data2mean, nr@data2quantile)
    #         methodnames <- c("Log2", "Loess-R", "RLR-R", "VSN-R", "Loess-G", "RLR-G", "VSN-G", "TI-G", "MedI-G", "AI-G", "Quantile")
    #     }
    #     # -----------------
    #     
    # }

    for (sampleIndex in 1:length(methodnames)) {
        
        print(slotNameList)
        print(slotNameList[[sampleIndex]])
        
        # write.table(file=paste(jobdir, "/", methodnames[sampleIndex], "-normalized.txt", sep=""), 
        #             cbind(nds@rawData[-(1:2), (1:(length(nds@inputHeaderValues) - length(nds@sampleReplicateGroups)))], 
        #                   methodlist[[sampleIndex]]), sep="\t", row.names=F, col.names=nds@rawData[2,], quote=F)
        
        write.table(file=paste(jobdir, "/", methodnames[sampleIndex], "-normalized.txt", sep=""),
                    cbind(nds@rawData[-(1:2), (1:(length(nds@inputHeaderValues) - length(nds@sampleReplicateGroups)))],
                          slot(nr, slotNameList[[sampleIndex]]), sep="\t", row.names=F, col.names=nds@rawData[2,], quote=F))
    }
    
    if (!all(is.null(nr@houseKeepingVars))) {
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


