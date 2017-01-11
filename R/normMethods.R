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
    
    # print(head(annotationColumns))
    # stop("")
    
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


