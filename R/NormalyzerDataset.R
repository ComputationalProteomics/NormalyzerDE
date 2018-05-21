#' S4 class to represent dataset information
#' 
#' @slot jobName Name of the job represented by the dataset.
#' @slot rawData Matrix with raw values.
#' @slot sampleNameCol Name column for sample.
#' @slot groupNameCol Name column for groups.
#' @slot designMatrix Data frame containing design.
#' @slot sampleNames Vector containing sample names.
#' @slot filterrawdata Reduced raw data matrix where low abundance rows are 
#'  removed - TODO: CHECK, IS THIS CURRENTLY WORKING?
#' @slot sampleReplicateGroups Vector with sample replicate information
#' @slot colsum Vector with sum of values for each column
#' @slot medofdata Vector with median values for each column
#' @slot meanofdata Vector with mean values for each column
#' @slot annotationValues Annotation part of original dataframe.
#' @slot retentionTimes Vector of retention time values.
#' @slot singleReplicateRun Conditional whether run is single replicate.
#' @export
NormalyzerDataset <- setClass("NormalyzerDataset",
                              representation(
                                  
                                  jobName = "character",
                                  rawData = "matrix",
                                  
                                  sampleNameCol = "character",
                                  groupNameCol = "character",
                                  
                                  designMatrix = "data.frame",
                                  sampleNames = "character",

                                  filterrawdata = "matrix",
                                  sampleReplicateGroups = "numeric",
                                  
                                  colsum = "numeric",
                                  medofdata = "numeric",
                                  meanofdata = "numeric",
                                  
                                  annotationValues = "matrix",
                                  retentionTimes = "numeric",
                                                                    
                                  singleReplicateRun = "logical"
                              ),
                              prototype=prototype(jobName=NULL, 
                                                  rawData=NULL,
                                                  designMatrix=NULL,
                                                  sampleNameCol=NULL,
                                                  groupNameCol=NULL))

#' Initialize values for dataset
#'
#' @param nds Normalyzer dataset
#' @return None
#' @rdname setupValues
#' @keywords internal
setGeneric(name="setupValues", function(nds) standardGeneric("setupValues"))

#' @rdname setupValues
setMethod("setupValues", "NormalyzerDataset",
          function(nds) {
              
              nds@sampleReplicateGroups <- as.factor(nds@designMatrix[, nds@groupNameCol])
              nds@sampleNames <- nds@designMatrix[, nds@sampleNameCol]

              singleReplicateRun <- detectSingleReplicate(nds)
              singletonSamplePresent <- detectSingletonSample(nds)
              
              if (singleReplicateRun || singletonSamplePresent) {
                  nds@singleReplicateRun <- TRUE
              }
              else {
                  nds@singleReplicateRun <- FALSE
              }
              
              nds <- setupBasicValues(nds)
              nds <- setupRTColumn(nds)
              
              nds
          }
)


#' Detect single replicate, and assign related logical
#'
#' @param nds Normalyzer dataset
#' @return bool on whether sample contains only one sample group
#' @rdname detectSingleReplicate
#' @keywords internal
setGeneric(name="detectSingleReplicate", 
           function(nds) standardGeneric("detectSingleReplicate"))

#' @rdname detectSingleReplicate
setMethod("detectSingleReplicate", "NormalyzerDataset",
          function(nds) {
              
              headerCounts <- table(nds@sampleReplicateGroups)
              nonReplicatedSamples <- names(headerCounts[which(headerCounts == 1)])
              
              if (length(nonReplicatedSamples) > 0) {
                  singleReplicateRun <- TRUE
                  print(paste("Non replicated samples in dataset:",
                              paste(nonReplicatedSamples, collapse=" "),
                              "performing limited single-replicate run"))
              }
              else {
                  singleReplicateRun <- FALSE
              }
              
              singleReplicateRun
          }
)

#' Detect single sample group 
#'
#' @param nds Normalyzer dataset.
#' @return None
#' @rdname detectSingletonSample
#' @keywords internal
setGeneric(name="detectSingletonSample", 
           function(nds) standardGeneric("detectSingletonSample"))

#' @rdname detectSingletonSample
setMethod("detectSingletonSample", "NormalyzerDataset",
          function(nds) {
              
              groups <- nds@sampleReplicateGroups
              distinctSamples <- unique(groups)
              singletonSamplePresent <- FALSE
              
              if (length(distinctSamples) == 1) {
                  singletonSamplePresent <- TRUE
                  print(paste("Only one replicate group present.",
                              paste("Group: ", distinctSamples[1]),
                              "Proceeding with limited processing",
                              sep="\n"))
              }
              else if (length(distinctSamples) == 0) {
                  stop("No sample replicate groups found - Aborting.")
              }
              
              singletonSamplePresent
          }
)

#' Setup basic values for dataset
#'
#' @param nds Normalyzer dataset
#' @return Updated dataset object
#' @rdname setupBasicValues
#' @keywords internal
setGeneric(name="setupBasicValues", 
           function(nds) standardGeneric("setupBasicValues"))

#' @rdname setupBasicValues
setMethod("setupBasicValues", "NormalyzerDataset",
          function(nds) {
              
              nds@sampleReplicateGroups <- as.numeric(nds@sampleReplicateGroups)
              nds@annotationValues <- nds@rawData[, !(colnames(nds@rawData) %in% nds@sampleNames), drop=FALSE]
              nds <- setupFilterRawData(nds)

              nds@colsum <- colSums(nds@filterrawdata, na.rm=TRUE)
              nds@medofdata <- apply(nds@filterrawdata, 2, FUN="median", na.rm=TRUE)
              nds@meanofdata <- apply(nds@filterrawdata, 2, FUN="mean", na.rm=TRUE)
              
              nds
          }
)

#' Check for and set up RT column if single -1 annotation present
#'
#' @param nds Normalyzer dataset
#' @return Updated dataset object
#' @rdname setupRTColumn
#' @keywords internal
setGeneric(name="setupRTColumn", 
           function(nds) standardGeneric("setupRTColumn"))

#' @rdname setupRTColumn
setMethod("setupRTColumn", "NormalyzerDataset",
          function(nds) {
              
              rtColumns <- grep("\\bRT\\b", colnames(nds@rawData))
              
              if (length(rtColumns) > 1) {
                  error_message <- paste(
                      "Only able to handle single RT column (name containing RT standing by itself)",
                      "Please change name or remove unwanted RT columns")
                  stop(error_message)
              }
              else if (length(rtColumns) == 1) {
                  
                  print(paste0("RT annotation column found (", rtColumns, ")"))
                  
                  rtValues <- nds@rawData[, rtColumns]
                  nds@retentionTimes <- as.numeric(rtValues)
              }
              
              nds
          }
)


#' Setup filtered raw data slot
#'
#' @param nds Normalyzer dataset.
#' @return None
#' @rdname setupFilterRawData
#' @keywords internal
setGeneric(name="setupFilterRawData", 
           function(nds) standardGeneric("setupFilterRawData"))

#' @rdname setupFilterRawData
setMethod("setupFilterRawData", "NormalyzerDataset",
          function(nds) {

              stopifnot(!is.null(nds@rawData))
              
              nds@filterrawdata <- nds@rawData[, nds@sampleNames]
              class(nds@filterrawdata) <- "numeric"

              nds
          }
)

