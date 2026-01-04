#' Create empty directory for run
#' 
#' Creates a directory at provided path named to the jobname.
#' 
#' @param jobName Name of the run.
#' @param outputDir Path to directory where to create the output directory.
#' @return Path to newly created directory.
#' @export
#' @examples
#' setupJobDir("job_name", "path/to/outdir")
setupJobDir <- function(jobName, outputDir) {
    
    sanitizedJobName <- sanitizeJobName(jobName)

    if (is.null(outputDir)) {
        jobDir <- file.path(getwd(), sanitizedJobName)
    }
    else {
        jobDir <- file.path(outputDir, sanitizedJobName)
    }
    createDirectory(jobDir)
    
    jobDir
}

sanitizeJobName <- function(jobName) {
    
    if (length(jobName) < 1 || is.null(jobName) || is.na(jobName[1])) {
        stop("Invalid jobName: must be a non-empty character value")
    }
    
    sanitizedJobName <- gsub("[^[:alnum:]_-]", "_", as.character(jobName[1]))
    sanitizedJobName <- basename(sanitizedJobName)
    
    if (!nzchar(sanitizedJobName)) {
        stop("Invalid jobName: must contain at least one character after sanitization")
    }
    
    sanitizedJobName
}

#' Get dataframe with raw data column sorted on replicates
#' 
#' @param rawDataOnly Dataframe with unparsed input data matrix.
#' @param groups Vector containing condition levels.
#' @return rawData sorted on replicate
#' @keywords internal
getReplicateSortedData <- function(rawDataOnly, groups) {

    indexList <- getIndexList(groups)
    orderedIndices <- unlist(indexList[sort(names(indexList))])

    rawDataOnly[, orderedIndices, drop=FALSE]
}

#' Create directory, or return error if already present
#' 
#' @param targetPath Path where to attempt to create directory
#' @return None
#' @keywords internal
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
#' @keywords internal
elapsedSecondsBetweenSystimes <- function(start, end) {
    
    startSecond <- strtoi(format(start, "%s"))
    endSecond <- strtoi(format(end, "%s"))
    elapsed <- end - start
    elapsed
}

#' Return list containing vector positions of values in string
#' 
#' @param targetVector 
#' @return indexList List where key is condition level and values are indices
#'   for the condition
#' @keywords internal
getIndexList <- function(targetVector) {
    
    indexList <- list()
    uniqVals <- unique(targetVector)
    for (val in uniqVals) {
        indexList[[toString(val)]] <- which(targetVector == val)
    }
    indexList
}

#' Get contrast vector (TRUE/FALSE-values) indicating whether both at least
#' half values are present, and each sample has at least one non-NA value
#' 
#' @param dataMatrix Matrix with expression values for entities in replicate 
#'  samples.
#' @param replicateHeader Header showing how samples in matrix are replicated.
#' @param minCount Minimum number of required values present in samples.
#' @return Contrast vector
#' @keywords internal
getRowNAFilterContrast <- function(dataMatrix, 
                                   replicateHeader, 
                                   minCount=1) {
    
    replicatesHaveData <- rep(TRUE, nrow(dataMatrix))
    indexList <- getIndexList(replicateHeader)
    
    for (sampleIndex in seq_along(names(indexList))) {
        
        repVal <- names(indexList)[sampleIndex]
        cols <- indexList[[repVal]]
        
        nbrNAperReplicate <- rowSums(is.na(dataMatrix[, cols, drop=FALSE]))
        nbrReplicates <- length(cols)
        nbrNonNA <- nbrReplicates - nbrNAperReplicate
        replicatesHaveData <- (nbrNonNA >= minCount & replicatesHaveData)
    }
    
    replicatesHaveData
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
#' @keywords internal
setupTestData <- function(nSamples, nFeatures, rtMin=40, rtMax=80, mean=20, sd=4) {

    featureNames <- paste0("feature_", seq(1, nFeatures))
    sampleData <- matrix(stats::rnorm(nSamples * nFeatures, mean, sd), nFeatures, nSamples)
    rtData <- stats::runif(nFeatures, rtMin, rtMax)

    df <- data.frame(feature=featureNames, 
                     RT=rtData, 
                     as.data.frame(sampleData))

    colnames(df) <- c("feature", "RT", paste0("S", seq(1, nSamples)))
    df
}

#' General function for calculating percentage difference of average column
#' means in matrix
#' 
#' @param targetMat Matrix for which column means should be compared
#' @return percDiffVector Vector with percentage difference, where first element
#'   always will be 100
#' @keywords internal
calculatePercentageAvgDiffInMat <- function(targetMat) {
    
    calculatePercDiff <- function (sampleIndex, mat) {
        mean(mat[, sampleIndex]) * 100 / mean(mat[, 1])
    }
    
    percDiffVector <- vapply(
        seq_len(ncol(targetMat)), 
        calculatePercDiff,
        0,
        mat=targetMat)
    
    percDiffVector
}

#' Filter rows with lower than given number of replicates for any condition
#' 
#' @param df Dataframe with expression data to filter
#' @param groups Condition groups header
#' @param leastRep Minimum number of replicates in each group
#'   to retain
#' @return collDesignDf Reduced design matrix
#' @keywords internal
filterLowRep <- function(df, groups, leastRep = 2) {
  
    allReplicatesHaveValuesContrast <- function(row, groups, minCount) {
        names(row) <- groups
        repCounts <- table(names(stats::na.omit(row)))
        length(repCounts) == length(unique(groups)) && min(repCounts) >= minCount || minCount == 0
    }
    
    rowMeetThresContrast <- apply(
        df, 
        1,
        allReplicatesHaveValuesContrast, 
        groups = groups, 
        minCount = leastRep)
    
    filteredDf <- df[rowMeetThresContrast, , drop=FALSE]
    filteredDf
}

#' Impute one value in groups with only NAs if any group has more than fraction non-NA values
#' 
#' @param df Dataframe with expression data
#' @param groups Condition groups header
#' @param minFraction Minimum fraction with values in one group before imputing low value in other groups
#' @return imputedDf Imputed data frame
#' @keywords internal
imputeGroupValues <- function(df, groups, minFraction = 0.75) {
  
    if (!any(!is.na(df))) {
        return(df)
    }

    minValue <- min(df, na.rm = TRUE)
    groupIndices <- split(seq_along(groups), groups)
    groupTotalCounts <- vapply(groupIndices, length, integer(1))
    inputRowNames <- if (is.data.frame(df)) attr(df, "row.names") else rownames(df)
  
    imputeRow <- function(row, groupIndices, groupTotalCounts, minFraction, minVal) {
        groupNonNaCounts <- vapply(
            groupIndices,
            function(indices) sum(!is.na(row[indices])),
            integer(1)
        )

        if (any(groupNonNaCounts / groupTotalCounts >= minFraction) && any(groupNonNaCounts == 0)) {
            missingGroups <- names(groupNonNaCounts)[groupNonNaCounts == 0]
            for (groupLabel in missingGroups) {
                row[groupIndices[[groupLabel]][1]] <- minVal
            }
        }

        row
    }

    imputed <- apply(
        df, 
        1,
        imputeRow, 
        groupIndices = groupIndices,
        groupTotalCounts = groupTotalCounts,
        minFraction = minFraction,
        minVal = minValue
    )

    imputedMat <- t(imputed)
    colnames(imputedMat) <- colnames(df)
    if (!is.null(inputRowNames)) rownames(imputedMat) <- rownames(df)

    if (is.data.frame(df)) {
        imputedDf <- as.data.frame(imputedMat)
        if (!is.null(inputRowNames)) attr(imputedDf, "row.names") <- inputRowNames
        return(imputedDf)
    }
    imputedMat
}
