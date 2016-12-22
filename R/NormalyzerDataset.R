#' S4 class to represent dataset information
#' 
#' @slot jobName Name of the job represented by the dataset
#' @slot rawData Matrix with raw values
#' @slot filterrawdata Reduced raw data matrix where low abundance rows are 
#'  removed
#' @slot normfinderFilterRawData Reduced raw data matrix used for Normfinder
#'  normalization
#' @slot inputHeaderValues Vector with sample descriptions
#' @slot sampleReplicateGroups Vector with sample replicate information
#' @slot colsum Vector with sum of values for each column
#' @slot medofdata Vector with median values for each column
#' @slot meanofdata Vector with mean values for each column
#' @export

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
                              ),
                              prototype=prototype(jobName=NULL, rawData=NULL))

#' Initialize values for dataset
#'
#' @param nds Normalyzer dataset.
#' @return None
#' @rdname setupValues
#' @export
setGeneric(name="setupValues", function(nds) standardGeneric("setupValues"))

#' @rdname setupValues
setMethod("setupValues", "NormalyzerDataset",
          function(nds) {
              nds@inputHeaderValues <- nds@rawData[1, ]
              sampleRepStr <- nds@inputHeaderValues[-which(nds@inputHeaderValues < 1)]
              nds@sampleReplicateGroups <- as.numeric(sampleRepStr)
              
              nds <- setupFilterRawData(nds)
              nds <- setupNormfinderFilterRawData(nds)
              
              nds@colsum <- colSums(nds@filterrawdata, na.rm=TRUE)
              nds@medofdata <- apply(nds@filterrawdata, 2, 
                                     FUN="median", na.rm=TRUE)
              nds@meanofdata <- apply(nds@filterrawdata, 2, 
                                      FUN="mean", na.rm=TRUE)
              
              nds
          }
)

#' Setup filtered raw data slot
#'
#' @param nds Normalyzer dataset.
#' @return None
#' @rdname setupFilterRawData
#' @export
setGeneric(name="setupFilterRawData", 
           function(nds) standardGeneric("setupFilterRawData"))

#' @rdname setupFilterRawData
setMethod("setupFilterRawData", "NormalyzerDataset",
          function(nds) {
              print("DEBUG: setupFilterRawData")
              
              stopifnot(!is.null(nds@rawData))
              
              
              fullSampleHead <- nds@inputHeaderValues
              filteredSampleHead <- nds@sampleReplicateGroups
              
              nbrNonRepInHeader <- length(fullSampleHead) - length(filteredSampleHead)
              
              annotTrimMatrix <- nds@rawData[, -(1:nbrNonRepInHeader)]
              
              colnames(annotTrimMatrix) <- nds@rawData[2, -(1:nbrNonRepInHeader)]
              filterRawData <- as.matrix(annotTrimMatrix[-(1:2), ])
              class(filterRawData) <- "numeric"
              
              nds@filterrawdata <- filterRawData
              
              nds
          }
)

#' Setup Normfinder filtered raw data
#'
#' @param nds Normalyzer dataset.
#' @return None
#' @rdname setupNormfinderFilterRawData
#' @export
setGeneric(name="setupNormfinderFilterRawData", 
           function(nds) standardGeneric("setupNormfinderFilterRawData"))

#' @rdname setupNormfinderFilterRawData
setMethod("setupNormfinderFilterRawData", "NormalyzerDataset",
          function(nds) {
              
              stopifnot(!is.null(nds@rawData))
              
              fullSampleHeader <- nds@inputHeaderValues
              filteredSampleHeader <- nds@sampleReplicateGroups
              
              countVal <- rowSums(!is.na(nds@rawData[, which(fullSampleHeader > 0)]))
              
              nbrReplicates <- 1 * ncol(nds@rawData[, which(fullSampleHeader > 0)])
              nds@normfinderFilterRawData <- nds@rawData[countVal >= nbrReplicates, ]
              
              nds
          }
)


