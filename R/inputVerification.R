
## Input: File path to raw data
## Output: Initialized normalization object
## Performs verification steps on input data, and provides informative error
## messages if failed to 
getVerifiedNormalyzerObjectFromFile <- function(inputFile) {

    rawData <- retrieveRawData(inputFile)
    # jobDir <- setupJobDir(outputDir, jobName)

    rawData <- getReplicateSortedData(rawData)
    verifySampleReplication(rawData)
    
    rawData <- preprocessData(rawData)
    
    normObj <- NormalyzerObject()
    inputHeaderValues <- rawData[1,]

    normObj@rawData <- rawData
    
    sampleReplicateGroups <- as.numeric(inputHeaderValues[-which(inputHeaderValues < 1)])
    normObj <- setupFilterRawData(normObj, rawData, inputHeaderValues, sampleReplicateGroups)
    
    normObj <- setupNormfinderFilterRawData(normObj, rawData, inputHeaderValues)
    
    normObj@inputHeaderValues <- inputHeaderValues
    normObj@sampleReplicateGroups <- sampleReplicateGroups
    
    # TODO: Implement empty column checking
    
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

# setupJobDir <- function(outputDir, jobName) {
#     
#     if (is.null(outputDir)) {
#         jobDir <- paste(getwd(), "/", jobName[1], sep="")
#     }
#     else {
#         jobDir <- paste(outputDir, "/", jobName[1], sep="")
#     }
#     createDirectory(jobDir)
#     jobDir
# }

# createDirectory <- function(targetPath) {
# 
#     if (file.exists(targetPath)) {
#         abc <- "Directory already exists"
#         class(abc) <- "try-error"
#         if (inherits(abc, "try-error")) {
#             return(abc)
#         }
#         stop("Directory already exists")
#     } 
#     else {
#         dir.create(targetPath)
#     }
# }

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