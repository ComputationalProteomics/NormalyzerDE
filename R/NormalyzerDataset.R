NormalyzerDataset <- setClass("NormalyzerDataset",
                             slots = c(
                                 rawData = "matrix",
                                 filterrawdata = "matrix",
                                 normfinderFilterRawData = "matrix",

                                 inputHeaderValues = "character",
                                 sampleReplicateGroups = "numeric",

                                 data2log2 = "matrix",
                                 data2GI = "matrix",
                                 data2ctr = "matrix",
                                 data2med = "matrix",
                                 data2mean = "matrix",
                                 data2ctrlog = "matrix",
                                 data2vsn = "matrix",
                                 data2quantile = "matrix",
                                 data2mad = "matrix",
                                 data2loess = "matrix",
                                 data2vsnrep = "matrix",
                                 data2limloess = "matrix",

                                 globalfittedRLR = "matrix",
                                 fittedLR = "matrix"
                             ),
                             prototype=prototype())

setGeneric(name="setupValues", function(x) standardGeneric("setupValues"))
setMethod("setupValues", "NormalyzerDataset",
          function(x) {
              x@inputHeaderValues <- x@rawData[1,]
              x@sampleReplicateGroups <- as.numeric(x@inputHeaderValues[-which(x@inputHeaderValues < 1)])
              x
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





