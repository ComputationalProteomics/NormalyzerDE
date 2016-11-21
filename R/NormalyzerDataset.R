NormalyzerDataset <- setClass("NormalyzerDataset",
                             slots = c(
                                 jobName = "character",
                                 
                                 rawData = "matrix",
                                 filterrawdata = "matrix",
                                 normfinderFilterRawData = "matrix",

                                 inputHeaderValues = "character",
                                 sampleReplicateGroups = "numeric",
                                 
                                 colsum = "numeric",
                                 medofdata = "numeric",
                                 meanofdata = "numeric"

                                 # data2log2 = "matrix",
                                 # data2GI = "matrix",
                                 # data2ctr = "matrix",
                                 # data2med = "matrix",
                                 # data2mean = "matrix",
                                 # data2ctrlog = "matrix",
                                 # data2vsn = "matrix",
                                 # data2quantile = "matrix",
                                 # data2mad = "matrix",
                                 # data2loess = "matrix",
                                 # data2vsnrep = "matrix",
                                 # data2limloess = "matrix",
                                 # 
                                 # globalfittedRLR = "matrix",
                                 # fittedLR = "matrix"
                             ),
                             prototype=prototype(jobName=NULL, rawData=NULL))

setGeneric(name="setupValues", function(nds) standardGeneric("setupValues"))
setMethod("setupValues", "NormalyzerDataset",
          function(nds) {
              nds@inputHeaderValues <- nds@rawData[1,]
              nds@sampleReplicateGroups <- as.numeric(nds@inputHeaderValues[-which(nds@inputHeaderValues < 1)])
              
              nds <- setupFilterRawData(nds)
              nds <- setupNormfinderFilterRawData(nds)
              
              nds@colsum <- colSums(nds@filterrawdata, na.rm=T)
              nds@medofdata <- apply(nds@filterrawdata, 2, FUN="median", na.rm=T)
              nds@meanofdata <- apply(nds@filterrawdata, 2, FUN="mean", na.rm=T)
              
              nds
          })

setGeneric(name="setupFilterRawData", function(nds) standardGeneric("setupFilterRawData"))
setMethod("setupFilterRawData", "NormalyzerDataset",
          function(nds) {
            print("DEBUG: setupFilterRawData")

            stopifnot(!is.null(nds@rawData))

            fullSampleHeader <- nds@inputHeaderValues
            filteredSampleHeader <- nds@sampleReplicateGroups
            
            filterrawdata <- nds@rawData[, -(1:(length(fullSampleHeader) - length(filteredSampleHeader)))]
            colnames(filterrawdata) <- nds@rawData[2, -(1:(length(fullSampleHeader) - length(filteredSampleHeader)))]
            filterrawdata <- (as.matrix((filterrawdata[-(1:2), ])))
            class(filterrawdata) <- "numeric"

            nds@filterrawdata <- filterrawdata

            nds
        })

setGeneric(name="setupNormfinderFilterRawData", function(nds) standardGeneric("setupNormfinderFilterRawData"))
setMethod("setupNormfinderFilterRawData", "NormalyzerDataset",
          function(nds) {

              stopifnot(!is.null(nds@rawData))

              fullSampleHeader <- nds@inputHeaderValues
              filteredSampleHeader <- nds@sampleReplicateGroups

              countna <- rowSums(!is.na(nds@rawData[, which(fullSampleHeader > 0)]))
              nds@normfinderFilterRawData <- nds@rawData[countna >= (1 * ncol(nds@rawData[, which(fullSampleHeader > 0)])), ]

              nds
          })


## --- Templates --- ##

## TestShow
## Template experimental method
# setGeneric(name="testshow", def=function(object){})
# setMethod(f="testshow", signature="NormalyzerDataset",
#     definition=function(object) {
#         cat("This is a class-method testprint!")
#     }
# )

# setGeneric(name="setupValues<-", function(x) standardGeneric("setupValues<-"))
# setMethod("setupValues", "NormalyzerObject",
#     function(x) {
#         x@inputHeaderValues <- x@rawData[1,]
#         x@sampleReplicateGroups <- as.numeric(x@inputHeaderValues[-which(x@inputHeaderValues < 1)])
#         
#         x
#     }
# )

# setGeneric(name="inputHeaderValues", function(x) standardGeneric("inputHeaderValues"))
# setMethod("inputHeaderValues", "NormalyzerDataset", function(x) x@inputHeaderValues)
# 
# setGeneric(name="inputHeaderValues<-", function(x, value) standardGeneric("inputHeaderValues<-"))
# setReplaceMethod("inputHeaderValues", "NormalyzerDataset", 
#     function(x, value) {
#         x@inputHeaderValues <- value; x
#     }
# )
# 
# setGeneric(name="sampleReplicateGroups", function(x) standardGeneric("sampleReplicateGroups"))
# setMethod("sampleReplicateGroups", "NormalyzerDataset", function(x) x@sampleReplicateGroups)
# 
# setGeneric(name="sampleReplicateGroups<-", function(x, value) standardGeneric("sampleReplicateGroups<-"))
# setReplaceMethod("sampleReplicateGroups", "NormalyzerDataset", function(x, value) x@sampleReplicateGroups <- value)

# setMethod(f="inputHeaderValues", signature="NormalyzerObject", definition=function(x) x@rawData[1,])
# 
# setGeneric(name="sampleReplicateGroups", function(x) standardGeneric("sampleReplicateGroups"))
# setMethod("sampleReplicateGroups", "NormalyzerObject", 
#     function(x) {
#         as.numeric(x@inputHeaderValues[-which(inputHeaderValues < 1)])
#     }
# )





