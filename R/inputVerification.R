#' Verify that input data is in correct format, and if so, return a generated
#'  Normalyzer data object from that input data
#' 
#' @param inputPath Path to file with raw input to Normalyzer.
#' @param jobName Name of ongoing run.
#' @return Normalyzer data object representing verified input data.
#' @export
getVerifiedNormalyzerObject <- function(inputPath, jobName, threshold=15, 
                                        omitSamples=FALSE,
                                        requireReplicates=TRUE) {

    rawData <- loadRawDataFromFile(inputPath)
    
    verifyValidNumbers(rawData)
    
    repSortedRawData <- getReplicateSortedData(rawData)
    processedRawData <- preprocessData(repSortedRawData)
    
    lowCountSampleFiltered <- getLowCountSampleFiltered(processedRawData, 
                                                        threshold=threshold, 
                                                        stopIfTooFew=!omitSamples)
    
    # If no samples left after omitting, stop
    verifyMultipleSamplesPresent(lowCountSampleFiltered, requireReplicates=requireReplicates)
    validateSampleReplication(lowCountSampleFiltered, requireReplicates=requireReplicates)
    
    nds <- generateNormalyzerDataset(lowCountSampleFiltered, jobName)
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
            print(paste0("Provided input file (", inputPath, ") not found:"))
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
verifyValidNumbers <- function(normalyzerDfAll) {
    
    normalyzerDf <- normalyzerDfAll[, which(as.numeric(normalyzerDfAll[1,]) > 0), 
                                    drop=FALSE]
    
    validPatterns <- c("\\d+(\\.\\d+)?", "NA", "null", "\\d+\\.\\d+e\\d+$")
    
    rawData <- normalyzerDf[-1:-2,]
    
    regexPattern <- sprintf("^(%s)$", paste(validPatterns, collapse="|"))
    matches <- grep(regexPattern, rawData, perl=TRUE, ignore.case=TRUE)
    nonMatches <- na.omit(rawData[-matches])
    
    if (length(nonMatches > 0)) {
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
    
    print("Input data checked. All fields are valid.")
}


#' Get dataframe with raw data column sorted on replicates
#' 
#' @param rawData Dataframe with unparsed Normalyzer input data.
#' @return rawData sorted on replicate
getReplicateSortedData <- function(rawData) {

    factor_levels <- sort(as.numeric(unique(rawData[1,])))
    temp_df <- NULL

    for (i in 1:length(factor_levels)) {
        temp_df <- cbind(temp_df, rawData[, which(rawData[1,] == factor_levels[as.numeric(i)]), drop=FALSE])
    }

    temp_df
}


#' Replace 0 values with NA in input data
#' 
#' @param rawData Dataframe with Normalyzer input data.
#' @return Parsed rawdata where 0 values are replaced with NA
preprocessData <- function(normalyzerDf) {

    dataMatrix <- normalyzerDf[-1:-2, ]
    dataMatrix[which(dataMatrix == 0)] <- NA
    processedDf <- rbind(normalyzerDf[1:2, ], dataMatrix)

    processedDf
}

#' Verify that samples contain at least a lowest number of values
#' 
#' @param dfWithNAs Dataframe with processed Normalyzer input data.
#'        Zero values are expected to have been replaced with NAs.
#' @param threshold Lowest number of allowed values in a column.
#' @return None
getLowCountSampleFiltered <- function(dfWithNAs, threshold=15, stopIfTooFew=TRUE) {
    
    rawData <- dfWithNAs[-1:-2,]
    header <- dfWithNAs[1,]
    sampleIndices <- which(as.numeric(header) > 0)
    
    numberOfValues <- vector(length=length(sampleIndices), mode="numeric")
    
    for (i in 1:length(sampleIndices)) {
        
        sampleIndex <- sampleIndices[i]
        numberOfValues[i] <- length(na.omit(rawData[,sampleIndex]))
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
        naSamplesOmittedDf <- dfWithNAs[, -sampleIndices[notPassingThreshold]]
    }
    else {
        dfWithNAs
    }
}


#' Check whether all samples have replicates
#' 
#' @param processedDf Prepared Normalyzer dataframe.
#' @param requireReplicates By default stops processing if not all samples
#'  have replicates
#' @return None
validateSampleReplication <- function(processedDf, requireReplicates=TRUE) {
    
    nonReplicatedSamples <- getNonReplicatedFromDf(processedDf)
    
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
verifyMultipleSamplesPresent <- function(processedDf, requireReplicates=TRUE) {
    
    header <- processedDf[1,]
    samples <- header[which(as.numeric(header) > 0)]
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
generateNormalyzerDataset <- function(rawData, jobName) {
    
    nds <- NormalyzerDataset(jobName=jobName, rawData=rawData)
    nds <- setupValues(nds)
    nds
}

