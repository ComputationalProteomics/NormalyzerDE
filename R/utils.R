#' Create empty directory for run if not already present
#' 
#' @param jobName Name of the run.
#' @param outputDir Path to directory where to create the output directory.
#' @return Path to newly created directory.
#' @export
#' @examples
#' setupJobDir("job_name", "path/to/outdir")
setupJobDir <- function(jobName, outputDir) {
    
    if (is.null(outputDir)) {
        jobDir <- paste(getwd(), "/", jobName[1], sep="")
    }
    else {
        jobDir <- paste(outputDir, "/", jobName[1], sep="")
    }
    createDirectory(jobDir)
    
    jobDir
}

#' Create directory, or return error if already present
#' 
#' @param targetPath Path where to attempt to create directory
#' @return None
createDirectory <- function(targetPath) {
    
    if (file.exists(targetPath)) {
        errorHandler <- "Directory already exists"
        class(errorHandler) <- "try-error"
        if (inherits(errorHandler, "try-error")) {
            return(errorHandler)
        }
        stop("Directory already exists")
    } 
    else {
        dir.create(targetPath, recursive=TRUE)
    }
}

#' Get number of seconds between two Sys.time() objects
#' 
#' @param start Start-time object
#' @param end End-time object
#' @return None
elapsedSecondsBetweenSystimes <- function(start, end) {
    
    startSecond <- strtoi(format(start, "%s"))
    endSecond <- strtoi(format(end, "%s"))
    elapsed <- end - start
    elapsed
}


getIndexList <- function(targetVector) {
    
    indexList <- list()
    uniqVals <- unique(targetVector)
    for (val in uniqVals) {
        indexList[[toString(val)]] <- which(targetVector == val)
    }
    indexList
}


#' Get contrast vector (TRUE/FALSE-values) indicating whether all samples
#' specified in replicateHeader have at least one non-NA value
#' 
#' @param dataMatrix Matrix with expression values for entities in replicate 
#'  samples.
#' @param replicateHeader Header showing how samples in matrix are replicated.
#' @param minCount Minimum number of non-NA required for sample.
#' @return Contrast vector
getLowRepCountFilterContrast <- function(dataMatrix, replicateHeader, minCount=1) {
    
    replicatesHaveData <- rep(TRUE, nrow(dataMatrix))
    indexList <- getIndexList(replicateHeader)
    
    for (sampleIndex in 1:length(names(indexList))) {
        
        repVal <- names(indexList)[sampleIndex]
        cols <- indexList[[repVal]]
        
        nbrNAperReplicate <- rowSums(is.na(dataMatrix[, cols, drop=FALSE]))
        nbrReplicates <- length(cols)
        nbrNonNA <- nbrReplicates - nbrNAperReplicate
        replicatesHaveData <- (nbrNonNA >= minCount & replicatesHaveData)
    }
    
    replicatesHaveData
}


#' Get contrast vector (TRUE/FALSE-values) indicating whether at least half
#' of values in row have non-NA values
#' 
#' @param dataMatrix Matrix with expression values for entities in replicate 
#'  samples.
#' @param replicateHeader Header showing how samples in matrix are replicated.
#' @return Contrast vector
getLowNALinesContrast <- function(dataMatrix, replicateHeader) {
    
    nbsNAperLine <- rowSums(is.na(dataMatrix))
    enoughNAContrast <- (nbsNAperLine < ncol(dataMatrix) / 2)
    
    enoughNAContrast
}


#' Get contrast vector (TRUE/FALSE-values) indicating whether both at least
#' half values are present, and each sample has at least one non-NA value
#' 
#' @param dataMatrix Matrix with expression values for entities in replicate 
#'  samples.
#' @param replicateHeader Header showing how samples in matrix are replicated.
#' @param varFilterFrac Variance filtering fraction.
#' @param minCount Minimum number of required values present in samples.
#' @return Contrast vector
getRowNAFilterContrast <- function(dataMatrix, 
                                   replicateHeader, 
                                   varFilterFrac=NULL,
                                   minCount=1) {
    
    lowNALinesContrast <- getLowNALinesContrast(dataMatrix, replicateHeader)
    samplesHaveValuesContrast <- getLowRepCountFilterContrast(dataMatrix, replicateHeader, minCount=minCount)
    
    if (!is.null(varFilterFrac)) {
        samplesPassVarThres <- getVarFilteredContrast(dataMatrix, varFilterFrac)
    }
    else {
        samplesPassVarThres <- rep(TRUE, nrow(dataMatrix))
    }
    
    return (lowNALinesContrast & samplesHaveValuesContrast & samplesPassVarThres)
}


getVarFilteredContrast <- function(logMatrix, varFilterFrac) {
  
    featureVars <- apply(logMatrix, 1, function(x) {stats::var(x, na.rm=TRUE)})
    fracThres <- sort(featureVars)[(length(featureVars) - 1)  * varFilterFrac + 1]
    featureVars > fracThres
}


