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
    # verifySampleReplication(repSortedRawData)  # This is the way to do it!
    processedRawData <- preprocessData(repSortedRawData)
    
    
    lowCountSampleFiltered <- getLowCountSampleFiltered(processedRawData, 
                                                        threshold=threshold, 
                                                        stopIfTooFew=!omitSamples)
    

    # If no samples left after omitting, stop
    
    # Verify replicates - Control whether to crash or progress
    validateSampleReplication(lowCountSampleFiltered, 
                              requireReplicates=requireReplicates)
    
    

    # nds <- generateNormalyzerDataset(processedRawData, jobName)
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
    
    normalyzerDf <- normalyzerDfAll[, which(normalyzerDfAll[1,] != "0"), 
                                    drop=FALSE]
    
    validPatterns <- c("\\d+(\\.\\d+)?", "NA")
    
    rawData <- normalyzerDf[-1:-2,]
    
    regexPattern <- sprintf("^(%s)$", paste(validPatterns, collapse="|"))
    matches <- grep(regexPattern, rawData, perl=TRUE, ignore.case=TRUE)
    nonMatches <- rawData[-matches]
    
    if (length(na.omit(nonMatches) > 0)) {
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

    temp_df <- NULL
    temp_df <- as.factor(rawData[1, ])
    factor_levels <- levels(temp_df)
    temp_df <- NULL
    
    for (i in 1:length(factor_levels)) {
        temp_df <- cbind(temp_df, rawData[, which(rawData[1, ] == factor_levels[as.numeric(i)])])
    }
    
    rawData <- temp_df
    rawData
}

# #' Perform validations that input data is in correct format, and abort
# #'  otherwise
# #' @param rawData Dataframe with unparsed Normalyzer input data.
# #' @return None
# #' @export
# verifySampleReplication <- function(rawData) {
# 
#     replicateLine <- rawData[1, ]
#     replicateLineUnique <- unique(replicateLine)
#     
#     for (i in 1:length(replicateLineUnique)) {
#         if (replicateLineUnique[i] != 0) {
#             if (length(grep(replicateLineUnique[i], replicateLine)) < 2) {
#                 
#                 error_object <- paste("Number of replicates are less than 2 for the group ", replicateLineUnique[i], sep="")
#                 class(error_object) <- "try-error"
#                 if (inherits(error_object, "try-error")) {
#                     return(error_object)
#                 }
#                 stop(paste("Number of replicates are less than 2 for the group ", replicateLineUnique[i], sep=""))
#             }
#         }
#     }
# }

#' Replace 0 values with NA in input data
#' 
#' @param rawData Dataframe with Normalyzer input data.
#' @return Parsed rawdata where 0 values are replaced with NA
preprocessData <- function(normalyzerDf) {

    # TODO: Is this problematic? What if you have zero as sample name?

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
    sampleIndices <- which(header != "0")
    
    numberOfValues <- vector(length=length(sampleIndices), mode="numeric")
    
    for (i in 1:length(sampleIndices)) {
        
        sampleIndex <- sampleIndices[i]
        numberOfValues[i] <- length(na.omit(rawData[,sampleIndex]))
        print(paste("Number found: ", numberOfValues[i]))
    }
    
    notPassingThreshold <- which(numberOfValues < threshold)
    
    if (length(notPassingThreshold) > 0) {
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
        dfWithNAs[, -sampleIndices[notPassingThreshold]]
    }
    else {
        dfWithNAs
    }
}


#' Check whether all samples have replicates
#' 
#' @param processedDf Prepared Normalyzer dataframe.
#' @return None
validateSampleReplication <- function(processedDf, requireReplicates=TRUE) {
    
    # header <- processedDf[1,]
    # sampleValues <- header[which(header != "0")]
    # headerCounts <- table(sampleValues)
    # 
    # nonReplicatedSamples <- names(headerCounts[which(headerCounts == 1)])
    
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

