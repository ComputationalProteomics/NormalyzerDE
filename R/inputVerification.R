#' Verify that input data is in correct format, and if so, return a generated
#'  Normalyzer data object from that input data
#' 
#' @param inputPath Path to file with raw input to Normalyzer.
#' @param jobName Name of ongoing run.
#' @return Normalyzer data object representing verified input data.
#' @export
#' @examples
#' getVerifiedNormalyzerObjectFromFile("path/to/mydata.csv", "my_run")
getVerifiedNormalyzerObjectFromFile <- function(inputPath, jobName) {

    # rawData <- retrieveRawData(inputPath)
    rawData <- as.matrix(read.table(inputPath, header=F, sep="\t", stringsAsFactors=F, quote=""))
    
    rawData <- getReplicateSortedData(rawData)
    verifySampleReplication(rawData)  # This is the way to do it!
    rawData <- preprocessData(rawData)
    
    nds <- generateNormalyzerDataset(rawData, jobName)

    # Verification step: Check whether all colsum values are > 0
    
    nds
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

#' Perform validations that input data is in correct format, and abort
#'  otherwise
#' 
#' @param rawData Dataframe with unparsed Normalyzer input data.
#' @return Vector with replicate line
#' @export
verifySampleReplication <- function(rawData) {

    replicateLine <- rawData[1, ]
    replicateLineUnique <- unique(replicateLine)
    
    for (i in 1:length(replicateLineUnique)) {
        if (replicateLineUnique[i] != 0) {
            if (length(grep(replicateLineUnique[i], replicateLine)) < 2) {
                
                error_object <- paste("Number of replicates are less than 2 for the group ", replicateLineUnique[i], sep="")
                class(error_object) <- "try-error"
                if (inherits(error_object, "try-error")) {
                    return(error_object)
                }
                stop(paste("Number of replicates are less than 2 for the group ", replicateLineUnique[i], sep=""))
            }
        }
    }
    replicateLine
}

#' Replace 0 values with NA in input data
#' 
#' @param rawData Dataframe with Normalyzer input data.
#' @return Parsed rawdata where 0 values are replaced with NA
preprocessData <- function(rawData) {

    rep0 <- rawData[-1,]
    rep0[which(rep0 == 0)] <- NA
    rawData <- rbind(rawData[1, ], rep0)
    
    rawData
}



