#' Load raw data into dataframe
#' 
#' @param dataPath File path to design matrix.
#' @param inputFormat If input is given in standard NormalyzerDE format, Proteios format
#'   or in MaxQuant protein or peptide format
#' @param zeroToNA Automatically convert zeroes to NA values
#' @return rawData Raw data loaded into data frame
#' @examples \dontrun{
#' df <- loadData("data.tsv")
#' @export
#' }
loadData <- function(dataPath, inputFormat="default", zeroToNA=FALSE) {
    
    if (inputFormat == "default") {
        rawData <- loadRawDataFromFile(dataPath)
    } else if (inputFormat == "proteios") {
        rawData <- proteiosToNormalyzer(dataPath)
    } else if (inputFormat == "maxquantpep") {
        rawData <- maxQuantToNormalyzer(dataPath, protLevel=FALSE)
    } else if (inputFormat == "maxquantprot") {
        rawData <- maxQuantToNormalyzer(dataPath, protLevel=TRUE)
    } else {
        valids <- c("default", "proteios", "maxquantpep", "maxquantprot")
        stop(paste("Unknown inputFormat:", inputFormat, "valids are:", paste(valids, collapse=", ")))
    }
    
    if (zeroToNA) {
        rawData[rawData == "0"] <- NA
    }
    
    rawData
}

#' Load raw design into dataframe
#' 
#' @param designPath File path to design matrix.
#' @param sampleCol Column name for column containing sample names.
#' @param groupCol Column name for column containing condition levels.
#' @return designMatrix Design data loaded into data frame
#' @examples \dontrun{
#' df <- loadDesign("design.tsv")
#' }
#' @export
loadDesign <- function(designPath, sampleCol="sample", groupCol="group") {
    designMatrix <- utils::read.table(designPath, sep="\t", stringsAsFactors=FALSE, header=TRUE, comment.char="")
    designMatrix[, sampleCol] <- as.character(designMatrix[, sampleCol])
    designMatrix[, groupCol] <- as.factor(as.character(designMatrix[, groupCol]))
    designMatrix
}

#' Verify that input data is in correct format, and if so, return a generated
#'  Normalyzer data object from that input data
#' 
#' @param jobName Name of ongoing run.
#' @param designMatrix Data frame containing design data.
#' @param rawData Data frame containing raw data.
#' @param threshold Minimum number of features.
#' @param omitSamples Automatically omit invalid samples from analysis.
#' @param requireReplicates Require there to be at least to samples per
#'        condition
#' @param sampleCol Name of the column containing samples in design matrix.
#' @param groupCol Name of column containing conditions in design matrix.
#' @param quiet Don't print output messages during processing 
#' @return Normalyzer data object representing verified input data.
#' @export
#' @examples
#' data(example_data)
#' data(example_design)
#' normObj <- getVerifiedNormalyzerObject("job_name", example_design, example_data)
getVerifiedNormalyzerObject <- function(
        jobName, 
        designMatrix, 
        rawData,
        threshold=15, 
        omitSamples=FALSE,
        requireReplicates=TRUE,
        sampleCol="sample",
        groupCol="group",
        quiet=FALSE
    ) {

    if (!groupCol %in% colnames(designMatrix)) {
        stop(paste0("Given groupCol: '", groupCol, "' was not present among design matrix columns"))
    }

    groups <- as.numeric(as.factor(designMatrix[, groupCol]))
    samples <- as.character(designMatrix[, sampleCol])
    
    fullMatrix <- rawData[-1,]
    colnames(fullMatrix) <- rawData[1,]

    verifyDesignMatrix(fullMatrix, designMatrix, sampleCol=sampleCol)
    dataMatrix <- fullMatrix[, samples]
    annotationMatrix <- fullMatrix[, -which(colnames(fullMatrix) %in% samples)]

    verifyValidNumbers(dataMatrix, groups, quiet=quiet)
    
    repSortedRawData <- getReplicateSortedData(dataMatrix, groups)
    processedRawData <- preprocessData(repSortedRawData)
    
    lowCountSampleFiltered <- getLowCountSampleFiltered(processedRawData, 
                                                        groups,
                                                        threshold=threshold, 
                                                        stopIfTooFew=!omitSamples)
    
    designMatrix <- designMatrix[which(designMatrix$sample %in% colnames(lowCountSampleFiltered)), ]
    
    # If no samples left after omitting, stop
    verifyMultipleSamplesPresent(lowCountSampleFiltered, groups, requireReplicates=requireReplicates, quiet=quiet)
    validateSampleReplication(lowCountSampleFiltered, groups, requireReplicates=requireReplicates, quiet=quiet)
    
    nds <- generateNormalyzerDataset(jobName, designMatrix, cbind(annotationMatrix, lowCountSampleFiltered), sampleCol, groupCol, quiet=quiet)
    nds
}

#' Try reading raw Normalyzer matrix from provided filepath
#' 
#' @param inputPath Path to Normalyzer data.
#' @return Table containing raw data from input file.
#' @keywords internal
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
#' @param rawDataOnly Dataframe with input data.
#' @param groups Condition levels for comparisons.
#' @return Parsed rawdata where 0 values are replaced with NA
#' @keywords internal
verifyValidNumbers <- function(rawDataOnly, groups, quiet=FALSE) {
    
    # Fields expected to contain numbers in decimal or scientific notation, or containing NA or null
    validPatterns <- c("\\d+(\\.\\d+)?", "NA", "\"NA\"", "null", "\\d+(\\.\\d+)?[eE]([\\+\\-])?\\d+$")
    
    regexPattern <- sprintf("^(%s)$", paste(validPatterns, collapse="|"))
    matches <- grep(regexPattern, rawDataOnly, perl=TRUE, ignore.case=TRUE)
    nonMatches <- stats::na.omit(rawDataOnly[-matches])
    
    if (length(nonMatches) > 0) {
        errorString <- paste(
            "Invalid values encountered in input data.",
            "Only valid data is numeric and NA- or na-fields",
            "Invalid fields: ",
            paste(unique(nonMatches), collapse=" "),
            "Aborting...",
            sep="\n"
        )
        
        stop(errorString)
    }
    
    zeroRegexPattern <- c("^0$")
    zeroMatches <- grep(zeroRegexPattern, rawDataOnly, perl=TRUE)
    if (length(zeroMatches) > 0) {
        errorString <- paste(
            "Encountered zeroes in data. Must be replaced with NA before processing.",
            "This can be done automatically setting the zeroToNA-flag option to TRUE.",
            sep="\n"
        )
        stop(errorString)
    }

    if (!quiet) print("Input data checked. All fields are valid.")
}

#' Verify that design matrix setup matches the data matrix
#' 
#' @param fullMatrix Dataframe with input data.
#' @param designMatrix Dataframe with design setup.
#' @param sampleCol Column in design matrix containing sample IDs.
#' 
#' @return None
#' @keywords internal
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
#' @param rawDataOnly Dataframe with unparsed input data matrix.
#' @param groups Vector containing condition levels.
#' @return rawData sorted on replicate
#' @keywords internal
getReplicateSortedData <- function(rawDataOnly, groups) {

    factorLevels <- sort(as.numeric(unique(groups)))
    tempDf <- NULL

    for (i in seq_len(length(factorLevels))) {
        tempDf <- cbind(tempDf, rawDataOnly[, which(groups == factorLevels[as.numeric(i)]), drop=FALSE])
    }

    tempDf
}


#' Replace 0 values with NA in input data
#' 
#' @param dataMatrix Matrix with raw data.
#' @return Parsed rawdata where 0 values are replaced with NA
#' @keywords internal
preprocessData <- function(dataMatrix) {

    dataMatrix[which(dataMatrix == 0)] <- NA
    dataMatrix
}

#' Verify that samples contain at least a lowest number of values
#' 
#' @param dataMatrix Dataframe with processed input data.
#' @param groups Vector containing condition levels.
#' @param threshold Lowest number of allowed values in a column.
#' @param stopIfTooFew Abort run if lower than threshold number of values in
#'        column
#' @return None
#' @keywords internal
getLowCountSampleFiltered <- function(dataMatrix, groups, threshold=15, stopIfTooFew=TRUE) {
    
    sampleIndices <- seq_len(length(groups))
    numberOfValues <- colSums(!is.na(dataMatrix))
    notPassingThreshold <- which(numberOfValues < threshold)
    
    if (length(notPassingThreshold) == length(numberOfValues)) {
        errorString <- paste(
            "None of the samples had enough valid non-NA values",
            "Found number of non-NA values:",
            paste(numberOfValues[notPassingThreshold], collapse=" "),
            "Current threshold: ",
            threshold,
            "You could try lowering the threshold by adjusting \"sampleAbundThres\"",
            "Be aware that this will likely lead to downstream crashes.",
            sep="\n"
        )
        stop(errorString)
    }
    else if (length(notPassingThreshold) > 0) {
        errorString <- paste(
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
            stop(errorString)
        }
        else {
            warning(errorString)
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
#' @param dataMatrix Prepared matrix containing expression data.
#' @param groups Vector containing condition levels
#' @param requireReplicates By default stops processing if not all samples
#'  have replicates
#' @return None
#' @keywords internal
validateSampleReplication <- function(dataMatrix, groups, requireReplicates=TRUE, quiet=FALSE) {
    
    headerCounts <- table(groups)
    nonReplicatedSamples <- names(headerCounts[which(headerCounts == 1)])

    if (length(nonReplicatedSamples) > 0) {
        
        errorString <- paste(
            "Following samples does not have replicates:",
            paste(nonReplicatedSamples, collapse=" "),
            "By default this is not allowed.",
            "You can force limited processing of non-replicated data by",
            "specifying the \"requireReplicates\" option",
            sep="\n")
        
        if (requireReplicates) {
            stop(errorString)
        }
        else {
            warning(errorString)
        }
    }
    else {
        if (!quiet) print("Sample replication check: All samples have replicates")
    }
}


#' Check whether more than one sample is present
#' 
#' @param dataMatrix Prepared dataframe.
#' @param groups Vector containing condition levels
#' @param requireReplicates By default stops processing if not all samples
#'  have replicates
#' @return None
#' @keywords internal
verifyMultipleSamplesPresent <- function(dataMatrix, groups, requireReplicates=TRUE, quiet=FALSE) {
    
    samples <- groups[which(as.numeric(groups) > 0)]
    distinctSamples <- unique(samples)

    if (length(samples) < 2) {
        errorString <- paste(
            "At least two samples are required to run Normalyzer",
            "Here, we found:",
            paste(samples, collapse=" "),
            sep="\n")
        
        stop(errorString)
    }
            
    if (length(distinctSamples) == 1) {
        
        errorString <- paste(
            "Found less than two distinct samples. Following was found:",
            paste(distinctSamples, collapse=" "),
            "For full processing two samples are required",
            "You can force limited processing by turning of the", 
            "\"requireReplicates\" option",
            sep="\n")
        
        if (requireReplicates) {
            stop(errorString)
        }
        else {
            warning(errorString)
        }
    }
    else if (length(distinctSamples) == 0) {
        stop("No replicate groups found. Double check your input file,
             and check that your data haven't been filtered out in preceeding
             input validation steps.")
    }
    else {
        if (!quiet) print("Sample check: More than one sample group found")
    }
}


#' Setup Normalyzer dataset from given raw data
#' 
#' @param jobName Name of ongoing run.
#' @param designMatrix Dataframe containing condition matrix.
#' @param fullRawMatrix Dataframe with unparsed input data.
#' @param sampleNameCol Name of column in design matrix containing sample names.
#' @param groupNameCol Name of column in design matrix contaning conditions.
#' @return Data object representing loaded data.
#' @keywords internal
generateNormalyzerDataset <- function(jobName, designMatrix, fullRawMatrix, sampleNameCol, groupNameCol, quiet=FALSE) {
    
    nds <- NormalyzerDataset(jobName=jobName, rawData=fullRawMatrix, designMatrix=designMatrix, sampleNameCol=sampleNameCol, groupNameCol=groupNameCol)
    nds <- setupValues(nds, quiet=quiet)
    nds
}

