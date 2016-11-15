
## Input: File path to raw data
## Output: Initialized normalization object
## Performs verification steps on input data, and provides informative error
## messages if failed to 
getVerifiedNormalyzerObjectFromFile <- function(inputFile) {

    rawData <- retrieveRawData(inputFile)
    
    
    
    rawData <- getReplicateSortedData(rawData)
    verifySampleReplication(rawData)  # This is the way to do it!
    rawData <- preprocessData(rawData)
    
    normObj <- generateNormalyzerDataset(rawData)
    
    normObj <- setupFilterRawData(normObj, rawData, normObj@inputHeaderValues, normObj@sampleReplicateGroups)
    normObj <- setupNormfinderFilterRawData(normObj, rawData, normObj@inputHeaderValues)

    print(normObj)

    normObj
}

generateNormalyzerDataset <- function(rawData) {

    normObj <- NormalyzerDataset()
    normObj@rawData <- rawData
    normObj <- setupValues(normObj)
    

    normObj
}

retrieveRawData <- function(datafile) {

    if (class(datafile) == "character") {
        getrawdata <- as.matrix((read.table(datafile, header=F, sep="\t", stringsAsFactors=F, quote="")))    
    } 
    else if (class(datafile) == "data.frame") {
        getrawdata <- as.matrix(datafile)
    } 
    else if (class(datafile) == "matrix") {
        getrawdata <- datafile  
    }
    getrawdata
}

getReplicateSortedData <- function(getrawdata) {

    b <- NULL
    b <- as.factor(getrawdata[1,])
    l <- levels(b)
    b <- NULL
    for (i in 1:length(l)) {
        b <- cbind(b, getrawdata[, which(getrawdata[1,] == l[as.numeric(i)])])
    }
    getrawdata <- b
    getrawdata
}

## Evaluate whether the input format is in correct format
## Abort processing if this isn't the case
verifySampleReplication <- function(getrawdata) {

    checkrep <- getrawdata[1,]
    
    # print(paste("Inside parseforerrors, checkrep: ", checkrep))
    
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
preprocessData <- function(getrawdata) {

    rep0 <- getrawdata[-1,]
    rep0[which(rep0==0)] <- NA
    getrawdata <- rbind(getrawdata[1,], rep0)
    
    getrawdata
}

setupFilterRawData <- function(normObj, rawData, fullSampleHeader, filteredSampleHeader) {

    filterrawdata <- rawData[, -(1:(length(fullSampleHeader) - length(filteredSampleHeader)))]
    colnames(filterrawdata) <- rawData[2, -(1:(length(fullSampleHeader) - length(filteredSampleHeader)))]
    filterrawdata <- (as.matrix((filterrawdata[-(1:2), ])))
    class(filterrawdata) <- "numeric"
    
    normObj@filterrawdata <- filterrawdata
    
    normObj
}

setupNormfinderFilterRawData <- function(normObj, rawData, fullSampleHeader) {

    countna <- rowSums(!is.na(rawData[, which(fullSampleHeader > 0)]))
    normObj@normfinderFilterRawData <- rawData[countna >= (1 * ncol(rawData[, which(fullSampleHeader > 0)])), ]
    normObj
}