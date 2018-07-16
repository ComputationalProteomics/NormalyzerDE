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

#' Generate a random test dataset with features, sample values and retention times
#' 
#' @param nSamples Number of samples
#' @param nFeatures Number of features
#' @param rtMin Minimum retention time
#' @param rtMax Maximum retention time
#' @param mean Mean value for sample intensities
#' @param sd Standard deviation for sample intensities
#' @return Test dataset
#' @export
#' @examples
#' df <- setupTestData(6, 20)
#' df <- setupTestData(6, 20, mean=15, sd=1)
setupTestData <- function(nSamples, nFeatures, rtMin=40, rtMax=80, mean=20, sd=4) {

    featureNames <- paste0("feature_", seq(1, nFeatures))
    set.seed(37)
    sampleData <- matrix(rnorm(nSamples * nFeatures, mean, sd), nFeatures, nSamples)
    rtData <- runif(nFeatures, rtMin, rtMax)

    df <- data.frame(feature=featureNames, 
                     RT=rtData, 
                     as.data.frame(sampleData))

    colnames(df) <- c("feature", "RT", paste0("S", seq(1, nSamples)))
    df
}


#' Generate design matrix
#' 
#' @param nSamples Number of samples
#' @return Design matrix
#' @export
#' @examples
#' df <- setupTestDesign(6)
setupTestDesign <- function(nSamples) {
    
    df <- data.frame(sample=rep(NA, nSamples), condition=rep(NA, nSamples))
    for (i in 1:nSamples) {
        row <- c(sample=paste0("S", i), condition=as.character(as.factor(i %% 2)))
        df[i, ] <- row
    }
    df
}





