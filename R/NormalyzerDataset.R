#' Represents raw input data together with basic annotation information
#' 
#' Takes a job name, a data matrix, a design matrix as well as specification
#' of the group and sample columns in the design matrix. Provides the
#' basic representation of a dataset in the NormalyzerDE normalization part.
#' 
#' @slot jobName Name of the job represented by the dataset.
#' @slot rawData Matrix with raw values.
#' @slot sampleNameCol Name column for sample.
#' @slot groupNameCol Name column for groups.
#' @slot designMatrix Data frame containing design.
#' @slot sampleNames Vector containing sample names.
#' @slot filterrawdata Reduced raw data matrix where low abundance rows are 
#'  removed
#' @slot sampleReplicateGroups Vector with sample replicate information
#' @slot samplesGroupsWithReplicates Vector with replicated sample replicate information
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
                                  samplesGroupsWithReplicates = "numeric",
                                  
                                  annotationValues = "matrix",
                                  retentionTimes = "numeric",
                                                                    
                                  singleReplicateRun = "logical"
                              ),
                              prototype=prototype(jobName=NULL, 
                                                  rawData=NULL,
                                                  designMatrix=NULL,
                                                  sampleNameCol=NULL,
                                                  groupNameCol=NULL))

# Getters
setGeneric("jobName", function(object) { standardGeneric("jobName") })
setMethod("jobName", signature(object="NormalyzerDataset"), 
          function(object) { slot(object, "jobName") })

setGeneric("rawData", function(object) { standardGeneric("rawData") })
setMethod("rawData", signature(object="NormalyzerDataset"), 
          function(object) { slot(object, "rawData") })

setGeneric("sampleNameCol", function(object) { standardGeneric("sampleNameCol") })
setMethod("sampleNameCol", signature(object="NormalyzerDataset"), 
          function(object) { slot(object, "sampleNameCol") })

setGeneric("groupNameCol", function(object) { standardGeneric("groupNameCol") })
setMethod("groupNameCol", signature(object="NormalyzerDataset"), 
          function(object) { slot(object, "groupNameCol") })

setGeneric("designMatrix", function(object) { standardGeneric("designMatrix") })
setMethod("designMatrix", signature(object="NormalyzerDataset"), 
          function(object) { slot(object, "designMatrix") })

setGeneric("sampleNames", function(object) { standardGeneric("sampleNames") })
setMethod("sampleNames", signature(object="NormalyzerDataset"), 
          function(object) { slot(object, "sampleNames") })
setGeneric("sampleNames<-", function(object, value) { standardGeneric("sampleNames<-") })
setReplaceMethod("sampleNames", signature(object="NormalyzerDataset"), 
                 function(object, value) { 
                     slot(object, "sampleNames") <- value
                     validObject(object)
                     object
                     })

setGeneric("filterrawdata", function(object) { standardGeneric("filterrawdata") })
setMethod("filterrawdata", signature(object="NormalyzerDataset"), 
          function(object) { slot(object, "filterrawdata") })
setGeneric("filterrawdata<-", function(object, value) { standardGeneric("filterrawdata<-") })
setReplaceMethod("filterrawdata", signature(object="NormalyzerDataset"), 
                 function(object, value) { 
                     slot(object, "filterrawdata") <- value
                     validObject(object)
                     object
                 })

setGeneric("sampleReplicateGroups", function(object) { standardGeneric("sampleReplicateGroups") })
setMethod("sampleReplicateGroups", signature(object="NormalyzerDataset"), 
          function(object) { slot(object, "sampleReplicateGroups") })
setGeneric("sampleReplicateGroups<-", function(object, value) { standardGeneric("sampleReplicateGroups<-") })
setReplaceMethod("sampleReplicateGroups", signature(object="NormalyzerDataset"), 
                 function(object, value) { 
                     slot(object, "sampleReplicateGroups") <- value
                     validObject(object)
                     object
                 })

setGeneric("samplesGroupsWithReplicates", function(object) { standardGeneric("samplesGroupsWithReplicates") })
setMethod("samplesGroupsWithReplicates", signature(object="NormalyzerDataset"), 
          function(object) { slot(object, "samplesGroupsWithReplicates") })
setGeneric("samplesGroupsWithReplicates<-", function(object, value) { standardGeneric("samplesGroupsWithReplicates<-") })
setReplaceMethod("samplesGroupsWithReplicates", signature(object="NormalyzerDataset"), 
                 function(object, value) { 
                     slot(object, "samplesGroupsWithReplicates") <- value
                     validObject(object)
                     object
                 })

setGeneric("annotationValues", function(object) { standardGeneric("annotationValues") })
setMethod("annotationValues", signature(object="NormalyzerDataset"), 
          function(object) { slot(object, "annotationValues") })
setGeneric("annotationValues<-", function(object, value) { standardGeneric("annotationValues<-") })
setReplaceMethod("annotationValues", signature(object="NormalyzerDataset"), 
                 function(object, value) { 
                     slot(object, "annotationValues") <- value
                     validObject(object)
                     object
                 })

setGeneric("retentionTimes", function(object) { standardGeneric("retentionTimes") })
setMethod("retentionTimes", signature(object="NormalyzerDataset"), 
          function(object) { slot(object, "retentionTimes") })

setGeneric("singleReplicateRun", function(object) { standardGeneric("singleReplicateRun") })
setMethod("singleReplicateRun", signature(object="NormalyzerDataset"), 
          function(object) { slot(object, "singleReplicateRun") })
setGeneric("singleReplicateRun<-", function(object, value) { standardGeneric("singleReplicateRun<-") })
setReplaceMethod("singleReplicateRun", signature(object="NormalyzerDataset"), 
                 function(object, value) { 
                     slot(object, "singleReplicateRun") <- value
                     validObject(object)
                     object
                 })

#' Initialize values for dataset
#'
#' @param nds Normalyzer dataset
#' @return None
#' @rdname setupValues
#' @keywords internal
setGeneric(name="setupValues", function(nds, quiet) standardGeneric("setupValues"))

#' @rdname setupValues
setMethod("setupValues", "NormalyzerDataset",
          function(nds, quiet=FALSE) {
              
              sampleReplicateGroups(nds) <- as.numeric(as.factor(designMatrix(nds)[, groupNameCol(nds)]))
              samplesGroupsWithReplicates(nds) <- as.numeric(
                  names(table(sampleReplicateGroups(nds))[which(table(sampleReplicateGroups(nds)) > 1)]))
              sampleNames(nds) <- as.character(designMatrix(nds)[, sampleNameCol(nds)])

              singleReplicateRun <- detectSingleReplicate(nds)
              singletonSamplePresent <- detectSingletonSample(nds)
              
              if (singleReplicateRun || singletonSamplePresent) {
                  singleReplicateRun(nds) <- TRUE
              }
              else {
                  singleReplicateRun(nds) <- FALSE
              }
              
              nds <- setupBasicValues(nds)
              nds <- setupRTColumn(nds, quiet=quiet)
              
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
              
              headerCounts <- table(sampleReplicateGroups(nds))
              nonReplicatedSamples <- names(headerCounts[which(headerCounts == 1)])
              
              if (length(nonReplicatedSamples) > 0) {
                  singleReplicateRun <- TRUE
                  message("Non replicated samples in dataset: ",
                          paste(nonReplicatedSamples, collapse=" "),
                          " performing limited single-replicate run")
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
              
              groups <- sampleReplicateGroups(nds)
              distinctSamples <- unique(groups)
              singletonSamplePresent <- FALSE
              
              if (length(distinctSamples) == 1) {
                  singletonSamplePresent <- TRUE
                  message("Only one replicate group present. ",
                          paste("Group: ", distinctSamples[1]),
                          " Proceeding with limited processing")
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
              
              annotationValues(nds) <- rawData(nds)[, !(colnames(rawData(nds)) %in% sampleNames(nds)), drop=FALSE]
              nds <- setupFilterRawData(nds)

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
           function(nds, quiet) standardGeneric("setupRTColumn"))

#' @rdname setupRTColumn
setMethod("setupRTColumn", "NormalyzerDataset",
          function(nds, quiet=FALSE) {
              
              rtColumns <- grep("\\bRT\\b", colnames(rawData(nds)))
              
              if (length(rtColumns) > 1) {
                  errorMessage <- paste(
                      "Only able to handle single RT column (name containing RT standing by itself)",
                      "Please change name or remove unwanted RT columns")
                  stop(errorMessage)
              }
              else if (length(rtColumns) == 1) {
                  
                  if (!quiet) message("RT annotation column found (", rtColumns, ")")
                  
                  rtValues <- rawData(nds)[, rtColumns]
                  retentionTimes(nds) <- as.numeric(rtValues)
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

              stopifnot(!is.null(rawData(nds)))
              
              filterrawdata(nds) <- rawData(nds)[, sampleNames(nds)]
              class(filterrawdata(nds)) <- "numeric"

              nds
          }
)

