#' Perform normalizations on Normalyzer dataset
#' 
#' @param nds Normalyzer dataset object.
#' @param currentjob Name of the ongoing processing run.
#' @param jobdir Path to the output directory
#' @return Returns Normalyzer results object with performed analyzes assigned
#'  as attributes
#' @export
normMethods <- function(nds, currentjob, jobdir) {

    nr <- generateNormalyzerResultsObject(nds)
    nr <- performNormalizations(nr)
    
    methodnames <- getUsedMethodNames(nr)
    methodlist <- getNormalizationMatrices(nr)

    for (sampleIndex in 1:length(methodnames)) {
        
        utils::write.table(file=paste(jobdir, "/", 
                                      methodnames[sampleIndex], 
                                      "-normalized.txt", sep=""),
                    cbind(nds@rawData[-(1:2), 1:(length(nds@inputHeaderValues) - length(nds@sampleReplicateGroups))],
                          methodlist[[sampleIndex]]), sep="\t", row.names=FALSE, col.names=nds@rawData[2,], quote=FALSE)
    }
    
    if (!all(is.na(nr@houseKeepingVars))) {
        utils::write.table(file=paste(jobdir, "/housekeeping-variables.txt", sep=""), nr@houseKeepingVars, sep="\t", row.names=FALSE, col.names=nds@rawData[2,], quote=FALSE)
    }
    
    utils::write.table(file=paste(jobdir, "/submitted_rawdata.txt", sep=""), 
                cbind(nds@rawData[-(1:2), (1:(length(nds@inputHeaderValues) - length(nds@sampleReplicateGroups)))], nds@filterrawdata), sep="\t", row.names=FALSE,
                col.names=nds@rawData[2,], quote=FALSE)

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

#' Retrieve vector with tags for used global normalization methods
#' 
#' @param houseKeepingFlag Boolean telling whether house-keeing normalization 
#' is used
#' @return Vector with string names for normalization tags
getMethodNames <- function(houseKeepingFlag) {

    if (houseKeepingFlag) {
        methodnames <- c("Log2", "TI-G", "MedI-G", "AI-G", "NF-G")
    }
    else {
        methodnames <- c("Log2", "TI-G", "MedI-G", "AI-G")
    }
    methodnames
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


