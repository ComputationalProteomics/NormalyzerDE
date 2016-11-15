
## Input: File path to raw data
## Output: Initialized normalization object
## Performs verification steps on input data, and provides informative error
## messages if failed to 
getVerifiedNormalyzerObjectFromFile <- function(inputFile, jobName) {

    rawData <- retrieveRawData(inputFile)
    
    rawData <- getReplicateSortedData(rawData)
    verifySampleReplication(rawData)  # This is the way to do it!
    rawData <- preprocessData(rawData)
    
    nds <- generateNormalyzerDataset(rawData, jobName)

    # Verification step: Check whether all colsum values are > 0
    
    nds
}

generateNormalyzerDataset <- function(rawData, jobName) {

    nds <- NormalyzerDataset(jobName=jobName, rawData=rawData)
    nds <- setupValues(nds)
    # nds <- setupFilterRawData(nds)
    
    nds
}

retrieveRawData <- function(datafile) {

    if (class(datafile) == "character") {
        rawData <- as.matrix((read.table(datafile, header=F, sep="\t", stringsAsFactors=F, quote="")))    
    } 
    else if (class(datafile) == "data.frame") {
        rawData <- as.matrix(datafile)
    } 
    else if (class(datafile) == "matrix") {
        rawData <- datafile  
    }
    rawData
}

getReplicateSortedData <- function(rawData) {

    b <- NULL
    b <- as.factor(rawData[1,])
    l <- levels(b)
    b <- NULL
    for (i in 1:length(l)) {
        b <- cbind(b, rawData[, which(rawData[1,] == l[as.numeric(i)])])
    }
    rawData <- b
    rawData
}

## Evaluate whether the input format is in correct format
## Abort processing if this isn't the case
verifySampleReplication <- function(rawData) {

    checkrep <- rawData[1,]
    
    repunique <- unique(checkrep)
    
    for (i in 1:length(repunique)) {
        if (repunique[i] != 0) {
            if (length(grep(repunique[i], checkrep)) < 2) {
                abc <- paste("Number of replicates are less than 2 for the group ", repunique[i], sep="")
                class(abc) <- "try-error"
                if (inherits(abc, "try-error")) {
                    return(abc)
                }
                stop(paste("Number of replicates are less than 2 for the group ", repunique[i], sep=""))
            }
        }
    }
    checkrep
}

## Replaces zeroes with NAs
preprocessData <- function(rawData) {

    rep0 <- rawData[-1,]
    rep0[which(rep0==0)] <- NA
    rawData <- rbind(rawData[1,], rep0)
    
    rawData
}



