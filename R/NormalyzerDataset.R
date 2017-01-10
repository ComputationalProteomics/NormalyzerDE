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
                                  meanofdata = "numeric",
                                  
                                  retentionTimes = "numeric",
                                  annotationValues = "matrix",
                                  
                                  singleReplicateRun = "logical"
                              ),
                              prototype=prototype(jobName=NULL, 
                                                  rawData=NULL))

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
              
              
              
              # Detect single replicate, and assign related logical
              nonReplicatedSamples <- getNonReplicatedFromDf(nds@rawData)
              if (length(nonReplicatedSamples) > 0) {
                  nds@singleReplicateRun <- TRUE
                  print(paste("Non replicated samples in dataset:",
                        paste(nonReplicatedSamples, collapse=" "),
                        "performing limited single-replicate run",
                        sep="\n"))
              }
              else {
                  nds@singleReplicateRun <- FALSE
              }
              
              # Detect single sample group
              header <- nds@rawData[1,]
              distinctSamples <- unique(header[which(as.numeric(header) > 0)])
              if (length(distinctSamples) == 1) {
                  nds@singleReplicateRun <- TRUE
                  print(paste("Only one replicate group present.",
                              paste("Group: ", distinctSamples[1]),
                              "Proceeding with limited processing",
                              sep="\n"))
              }
              else if (length(distinctSamples) == 0) {
                  stop("No sample replicate groups found - Aborting.")
              }
              
              # Setup
              nds@inputHeaderValues <- nds@rawData[1, ]
              sampleRepStr <- nds@inputHeaderValues[which(as.numeric(nds@inputHeaderValues) > 0)]
              nds@sampleReplicateGroups <- as.numeric(sampleRepStr)
              
              nds@annotationValues <- nds@rawData[-1:-2, which(as.numeric(inputHeaderValues) < 1)]
              
              print(head(nds@rawData))
              print(inputHeaderValues)
              print(as.numeric(inputHeaderValues) < 1)
              print(head(nds@annotationValues))
              stop("")

              nds <- setupFilterRawData(nds)
              if (!nds@singleReplicateRun) {
                nds <- setupNormfinderFilterRawData(nds)
              }
              
              nds@colsum <- colSums(nds@filterrawdata, na.rm=TRUE)
              nds@medofdata <- apply(nds@filterrawdata, 2, FUN="median", na.rm=TRUE)
              nds@meanofdata <- apply(nds@filterrawdata, 2, FUN="mean", na.rm=TRUE)
              
              # If present - extract RT column
              rt_annotations <- nds@inputHeaderValues[which(inputHeaderValues == "-1")]
              
              if (length(rt_annotations) > 1) {
                  error_message <- paste(
                      "Only able to handle single RT column (marked with -1)",
                      "Please set unwanted RT columns to 0")
                  stop(error_message)
              }
              else if (length(rt_annotations) == 1) {
                  print("RT annotation column found")
                  rtValues <- nds@rawData[, which(inputHeaderValues == "-1")][-1:-2]
                  
                  # print(as.numeric(rtValues))
                  
                  nds@retentionTimes <- as.numeric(rtValues)
              }
              
              print(head(nds@rawData))
              print(head(nds@annotationValues))
              stop("")
              
              
              
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

              stopifnot(!is.null(nds@rawData))
              
              fullSampleHead <- nds@inputHeaderValues
              sampleDataWithHead <- nds@rawData[, which(as.numeric(fullSampleHead) > 0), drop=FALSE]
              
              filterRawData <- as.matrix(sampleDataWithHead[-(1:2),, drop=FALSE ])
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
                            
              nds@normfinderFilterRawData <- nds@rawData[countVal >= nbrReplicates,,drop=FALSE]

              nds
          }
)


