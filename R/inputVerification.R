#' Verify that input data is in correct format, and if so, return a generated
#'  Normalyzer data object from that input data
#' 
#' @param inputPath Path to file with raw input to Normalyzer.
#' @param jobName Name of ongoing run.
#' @return Normalyzer data object representing verified input data.
#' @export
getVerifiedNormalyzerObject <- function(inputPath, 
                                        jobName, 
                                        designMatrixPath=NULL,
                                        threshold=15, 
                                        omitSamples=FALSE,
                                        requireReplicates=TRUE,
                                        sampleCol="sample",
                                        groupCol="group",
                                        zeroToNA=FALSE,
                                        inputFormat="default") {

    if (inputFormat == "default") {
        rawData <- loadRawDataFromFile(inputPath)
    } else if (inputFormat == "proteios") {
        rawData <- proteoisToNormalyzer(inputPath)
    } else if (inputFormat == "maxquantpep") {
        rawData <- maxQuantToNormalyzer(inputPath, protLevel=FALSE)
    } else if (inputFormat == "maxquantprot") {
        rawData <- maxQuantToNormalyzer(inputPath, protLevel=TRUE)
    } else {
        valids <- c("default", "proteios", "maxquantpep", "maxquantprot")
        stop(paste("Unknown inputFormat:", inputFormat, "valids are:", paste(valids, collapse=", ")))
    }
    
    if (is.null(designMatrixPath)) {
        print("No design matrix providing, inferring it from legacy header")
        designMatrix <- legacyNormalyzerToDesign(inputPath)
        rawData <- rawData[-1,]
    }
    else {
        designMatrix <- read.table(designMatrixPath, sep="\t", stringsAsFactors=F, header=T)
    }
    
    if (zeroToNA) {
        rawData[rawData == "0"] <- NA
    }
    
    groups <- as.numeric(as.factor(designMatrix[, groupCol]))
    fullMatrix <- rawData[-1,]
    colnames(fullMatrix) <- rawData[1,]
    
    verifyDesignMatrix(fullMatrix, designMatrix, sampleCol=sampleCol)
    dataMatrix <- fullMatrix[, designMatrix[, sampleCol]]
    
    verifyValidNumbers(dataMatrix, groups)
    
    repSortedRawData <- getReplicateSortedData(dataMatrix, groups)
    processedRawData <- preprocessData(repSortedRawData)
    
    # processedRawData <- preprocessData(dataMatrix)

    lowCountSampleFiltered <- getLowCountSampleFiltered(processedRawData, 
                                                        groups,
                                                        threshold=threshold, 
                                                        stopIfTooFew=!omitSamples)
    
    # If no samples left after omitting, stop
    verifyMultipleSamplesPresent(lowCountSampleFiltered, groups, requireReplicates=requireReplicates)
    validateSampleReplication(lowCountSampleFiltered, groups, requireReplicates=requireReplicates)
    
    nds <- generateNormalyzerDataset(rawData, jobName, designMatrix, sampleCol, groupCol)
    nds
}

#' Try reading raw Normalyzer matrix from provided filepath
#' 
#' @param inputPath Path to Normalyzer data.
#' @return Table containing raw data from input file.
loadRawDataFromFile <- function(inputPath) {
    
    tryCatch(
        rawData <- as.matrix(utils::read.table(inputPath, 
                                               header=FALSE, 
                                               sep="\t", 
                                               stringsAsFactors=FALSE,
                                               quote="",
                                               comment.char="")),
        error=function(e) {
            print(paste0("Error encountered for input file:", inputPath, ", error: ", e))
            stop("Please provide a valid input file.")
        },
        
        warning=function(w) {
            print(paste0("Provided input file (", inputPath, ") not found:"))
            stop("Please provide a valid input file.")
        }
    )
    
    rawData
}


#' Verify that input fields conform to the expected formats
#' 
#' @param rawData Dataframe with Normalyzer input data.
#' @return Parsed rawdata where 0 values are replaced with NA
verifyValidNumbers <- function(rawDataOnly, groups) {
    
    # Fields expected to contain numbers in decimal or scientific notation, or containing NA or null
    validPatterns <- c("\\d+(\\.\\d+)?", "NA", "null", "\\d+(\\.\\d+)?[eE]([\\+\\-])?\\d+$")
    
    regexPattern <- sprintf("^(%s)$", paste(validPatterns, collapse="|"))
    matches <- grep(regexPattern, rawDataOnly, perl=TRUE, ignore.case=TRUE)
    nonMatches <- na.omit(rawDataOnly[-matches])
    
    # browser()
    
    if (length(nonMatches) > 0) {
        error_string <- paste(
            "Invalid values encountered in input data.",
            "Only valid data is numeric and NA- or na-fields",
            "Invalid fields: ",
            paste(unique(nonMatches), collapse=" "),
            "Aborting...",
            sep="\n"
        )
        
        stop(error_string)
    }
    
    zeroRegexPattern <- c("^0$")
    zeroMatches <- grep(zeroRegexPattern, rawDataOnly, perl=TRUE)
    if (length(zeroMatches) > 0) {
        error_string <- paste(
            "Encountered zeroes in data. Must be replaced with NA before processing.",
            "This can be done automatically setting the zeroToNA-flag option to TRUE.",
            sep="\n"
        )
        stop(error_string)
    }

    print("Input data checked. All fields are valid.")
}

#' Verify that design matrix setup matches the data matrix
#' 
#' @param rawData Dataframe with Normalyzer input data
#' @param designMatrix Dataframe with design setup for Normalyzer
#' 
#' @return None
verifyDesignMatrix <- function(fullMatrix, designMatrix, sampleCol="sample") {

    designColnames <- designMatrix[, sampleCol]
    
    if (!all(designColnames %in% colnames(fullMatrix))) {
        errorString <- paste(
            "Not all columns present in design matrix are present in data matrix",
            "The following element was not found:",
            paste(setdiff(designColnames, colnames(fullMatrix)), collapse=" ")
        )
        stop(errorString)
    }

    dataMatrix <- fullMatrix[, designMatrix[, sampleCol]]
    dataColumns <- dataMatrix[, designColnames]
    
    if (length(designColnames) != ncol(dataColumns)) {
        errorString <- paste(
            "Number of samples does not match number of selected columns",
            "Found number of columns:",
            ncol(dataColumns),
            "Expected number of columns:",
            length(designColnames),
            "Are all columns in the design matrix present in the data matrix?"
        )
        stop(errorString)
    }
    
    if (length(unique(designColnames)) != length(designColnames)) {
        errorString <- paste(
            "Sample labels must be unique, found: ",
            length(unique(designColnames)),
            "unique, expected:",
            length(designColnames)
        )
        stop(errorString)
    }
}


#' Get dataframe with raw data column sorted on replicates
#' 
#' @param rawData Dataframe with unparsed Normalyzer input data.
#' @return rawData sorted on replicate
getReplicateSortedData <- function(rawDataOnly, groups) {

    factorLevels <- sort(as.numeric(unique(groups)))
    tempDf <- NULL

    for (i in 1:length(factorLevels)) {
        tempDf <- cbind(tempDf, rawDataOnly[, which(groups == factorLevels[as.numeric(i)]), drop=FALSE])
    }

    tempDf
}


#' Replace 0 values with NA in input data
#' 
#' @param rawData Dataframe with Normalyzer input data.
#' @return Parsed rawdata where 0 values are replaced with NA
preprocessData <- function(dataMatrix) {

    dataMatrix[which(dataMatrix == 0)] <- NA

    dataMatrix
}

#' Verify that samples contain at least a lowest number of values
#' 
#' @param dfWithNAs Dataframe with processed Normalyzer input data.
#'        Zero values are expected to have been replaced with NAs.
#' @param threshold Lowest number of allowed values in a column.
#' @return None
getLowCountSampleFiltered <- function(dataMatrix, groups, threshold=15, stopIfTooFew=TRUE) {
    
    sampleIndices <- which(as.numeric(groups) > 0)
    numberOfValues <- vector(length=length(sampleIndices), mode="numeric")
    
    for (i in 1:length(sampleIndices)) {
        
        sampleIndex <- sampleIndices[i]
        numberOfValues[i] <- length(na.omit(dataMatrix[, sampleIndex]))
    }
    
    notPassingThreshold <- which(numberOfValues < threshold)
    
    if (length(notPassingThreshold) == length(numberOfValues)) {
        error_string <- paste(
            "None of the samples had enough valid non-NA values",
            "Found number of non-NA values:",
            paste(numberOfValues[notPassingThreshold], collapse=" "),
            "Current threshold: ",
            threshold,
            "You could try lowering the threshold by adjusting \"sampleAbundThres\"",
            "Be aware that this will likely lead to downstream crashes.",
            sep="\n"
        )
        stop(error_string)
    }
    else if (length(notPassingThreshold) > 0) {
        error_string <- paste(
            "Following samples does not contain enough non-NA values:",
            paste(sampleIndices[notPassingThreshold], collapse=" "),
            "Found number of values:",
            paste(numberOfValues[notPassingThreshold], collapse=" "),
            "Current threshold:",
            threshold,
            "You can force processing without this sample by specifying the",
            "option \"omitLowAbundSamples\"",
            sep="\n")
        
        if (stopIfTooFew) {
            stop(error_string)
        }
        else {
            warning(error_string)
        }
    }
    
    if (length(notPassingThreshold > 0)) {
        naSamplesOmittedDf <- dataMatrix[, -sampleIndices[notPassingThreshold]]
    }
    else {
        dataMatrix
    }
}


#' Check whether all samples have replicates
#' 
#' @param processedDf Prepared Normalyzer dataframe.
#' @param requireReplicates By default stops processing if not all samples
#'  have replicates
#' @return None
validateSampleReplication <- function(dataMatrix, groups, requireReplicates=TRUE) {
    
    headerCounts <- table(groups)
    nonReplicatedSamples <- names(headerCounts[which(headerCounts == 1)])

    if (length(nonReplicatedSamples) > 0) {
        
        error_string <- paste(
            "Following samples does not have replicates:",
            paste(nonReplicatedSamples, collapse=" "),
            "By default this is not allowed.",
            "You can force limited processing of non-replicated data by",
            "specifying the \"requireReplicates\" option",
            sep="\n")
        
        if (requireReplicates) {
            stop(error_string)
        }
        else {
            warning(error_string)
        }
    }
    else {
        print("Sample replication check: All samples have replicates")
    }
}


#' Check whether more than one sample is present
#' 
#' @param processedDf Prepared Normalyzer dataframe.
#' @return None
verifyMultipleSamplesPresent <- function(dataMatrix, groups, requireReplicates=TRUE) {
    
    samples <- groups[which(as.numeric(groups) > 0)]
    distinctSamples <- unique(samples)

    if (length(samples) < 2) {
        error_string <- paste(
            "At least two samples are required to run Normalyzer",
            "Here, we found:",
            paste(samples, collapse=" "),
            sep="\n")
        
        stop(error_string)
    }
            
    if (length(distinctSamples) == 1) {
        
        error_string <- paste(
            "Found less than two distinct samples. Following was found:",
            paste(distinctSamples, collapse=" "),
            "For full processing two samples are required",
            "You can force limited processing by turning of the", 
            "\"requireReplicates\" option",
            sep="\n")
        
        if (requireReplicates) {
            stop(error_string)
        }
        else {
            warning(error_string)
        }
    }
    else if (length(distinctSamples) == 0) {
        stop("No replicate groups found. Double check your input file,
             and check that your data haven't been filtered out in preceeding
             input validation steps.")
    }
    else {
        print("Sample check: More than one sample group found")
    }
}


#' Setup Normalyzer dataset from given raw data
#' 
#' @param rawData Dataframe with unparsed Normalyzer input data.
#' @param jobName Name of ongoing run.
#' @return Normalyzer data object representing loaded data.
generateNormalyzerDataset <- function(fullRawMatrix, jobName, designMatrix, sampleNameCol, groupNameCol) {
    
    rawData <- fullRawMatrix[-1,]
    colnames(rawData) <- fullRawMatrix[1,]
    
    nds <- NormalyzerDataset(jobName=jobName, rawData=rawData, designMatrix=designMatrix, sampleNameCol=sampleNameCol, groupNameCol=groupNameCol)
    nds <- setupValues(nds)
    nds
}

