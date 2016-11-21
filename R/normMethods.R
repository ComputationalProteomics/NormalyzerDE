normMethods <- function(nds, currentjob, jobdir) {

    nr <- generateNormalyzerResultsObject(nds)
    nr <- performNormalizations(nr)
    
    methodnames <- getMethodNames(nr)
    methodlist <- getNormalizationMatrices(nr)

    for (sampleIndex in 1:length(methodnames)) {
        
        write.table(file=paste(jobdir, "/", methodnames[sampleIndex], "-normalized.txt", sep=""),
                    cbind(nds@rawData[-(1:2), (1:(length(nds@inputHeaderValues) - length(nds@sampleReplicateGroups)))],
                          methodlist[[sampleIndex]]), sep="\t", row.names=F, col.names=nds@rawData[2,], quote=F)
    }
    
    if (!all(is.na(nr@houseKeepingVars))) {
        write.table(file=paste(jobdir, "/housekeeping-variables.txt", sep=""), nr@houseKeepingVars, sep="\t", row.names=F, col.names=nds@rawData[2,], quote=F)
    }
    
    write.table(file=paste(jobdir, "/submitted_rawdata.txt", sep=""), 
                cbind(nds@rawData[-(1:2), (1:(length(nds@inputHeaderValues) - length(nds@sampleReplicateGroups)))], nds@filterrawdata), sep="\t", row.names=F,
                col.names=nds@rawData[2,], quote=F)

    return(nr)
}

generateNormalyzerResultsObject <- function(nds) {
    nr <- NormalyzerResults(nds=nds)
    nr <- initializeResultsObject(nr)
    nr
}

getMethodList <- function(HKflag, normObj) {

    if (HKflag) {
        methodlist <- list(normObj@data2log2, normObj@data2GI, normObj@data2med, normObj@data2mean, normObj@data2ctrlog)
    }
    else {
        methodlist <- list(normObj@data2log2, normObj@data2GI, normObj@data2med, normObj@data2mean)
    }
    methodlist
}

getMethodNames <- function(houseKeepingFlag) {

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


