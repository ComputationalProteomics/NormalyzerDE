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
                                  
                                  isTinyRun = "logical",
                                  singleReplicateRun = "logical"
                              ))

#' Constructor for NormalyzerDataset
#' 
#' @param jobName Name of the NormalyzerDE processing run
#' @param designMatrix Matrix containing sample conditions
#' @param rawData Matrix containing raw input data
#' @param annotationData Matrix containing annotation information for each
#'   input feature. Is expected to contain the same number of rows as the data
#'   but can contain any number of features
#' @param sampleNameCol Name of column in design matrix containing sample
#'   information
#' @param groupNameCol Name of column in design matrix containing condition
#'   information
#' @param tinyRunThres If fewer features than this is present in the input
#'   a limited run will be performed to avoid some steps requiring a more
#'   extensive number of features.
#' @param quiet If set to TRUE no information messages will be printed
#' @return nds Generated NormalyzerDataset instance
#' @keywords internal
NormalyzerDataset <- function(jobName, designMatrix, rawData, annotationData,
                                sampleNameCol, groupNameCol, tinyRunThres=50, quiet=FALSE) {
              
              sampleReplicateGroups <- as.numeric(
                  as.factor(designMatrix[, groupNameCol]))
              samplesGroupsWithReplicates <- as.numeric(
                  names(table(sampleReplicateGroups)[table(sampleReplicateGroups) > 1]))
              sampleNames <- as.character(designMatrix[, sampleNameCol])
              
              # annotationValues <- rawData[, !(colnames(rawData) %in% sampleNames), drop=FALSE]
              filterrawdata <- rawData[, sampleNames]
              class(filterrawdata) <- "numeric"
              
              if (nrow(filterrawdata) < tinyRunThres) {
                  isTinyRun <- TRUE
                  if (!quiet) {
                      message("Number of features (", nrow(filterrawdata), 
                              ") is less than the theshold for normal run (", 
                              tinyRunThres, "), a limited run is performed ",
                              "\nExact settings can be adjusted using the 'tinyRunThreshold' setting")
                  }
              }
              else {
                  isTinyRun <- FALSE
              }
              
              nds <- new(
                  "NormalyzerDataset",
                  jobName=jobName,
                  designMatrix=designMatrix,
                  rawData=rawData,
                  sampleNameCol=sampleNameCol,
                  groupNameCol=groupNameCol,
                  filterrawdata=filterrawdata,
                  sampleReplicateGroups=sampleReplicateGroups,
                  samplesGroupsWithReplicates=samplesGroupsWithReplicates,
                  annotationValues=annotationData,
                  isTinyRun=isTinyRun
              )
              
              rtColumn <- getRTColumn(annotationData, quiet=quiet)
              if (!is.null(rtColumn)) {
                  retentionTimes(nds) <- rtColumn
              }
              
              singleReplicateRun(nds) <- checkSingleReplicateRun(nds, quiet=quiet)
              
              nds
          }

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

setGeneric("filterrawdata", function(object) { standardGeneric("filterrawdata") })
setMethod("filterrawdata", signature(object="NormalyzerDataset"), 
          function(object) { slot(object, "filterrawdata") })

setGeneric("sampleReplicateGroups", function(object) { standardGeneric("sampleReplicateGroups") })
setMethod("sampleReplicateGroups", signature(object="NormalyzerDataset"), 
          function(object) { slot(object, "sampleReplicateGroups") })

setGeneric("samplesGroupsWithReplicates", function(object) { standardGeneric("samplesGroupsWithReplicates") })
setMethod("samplesGroupsWithReplicates", signature(object="NormalyzerDataset"), 
          function(object) { slot(object, "samplesGroupsWithReplicates") })

setGeneric("annotationValues", function(object) { standardGeneric("annotationValues") })
setMethod("annotationValues", signature(object="NormalyzerDataset"), 
          function(object) { slot(object, "annotationValues") })

setGeneric("retentionTimes", function(object) { standardGeneric("retentionTimes") })
setMethod("retentionTimes", signature(object="NormalyzerDataset"), 
          function(object) { slot(object, "retentionTimes") })
setGeneric("retentionTimes<-", function(object, value) { standardGeneric("retentionTimes<-") })
setReplaceMethod("retentionTimes", signature(object="NormalyzerDataset"), 
                 function(object, value) { 
                     slot(object, "retentionTimes") <- value
                     validObject(object)
                     object
                 })

setGeneric("isTinyRun", function(object) { standardGeneric("isTinyRun") })
setMethod("isTinyRun", signature(object="NormalyzerDataset"), 
          function(object) { slot(object, "isTinyRun") })

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

#' Detect single replicate, and assign related logical
#'
#' @param nds Normalyzer dataset
#' @param quiet Don't give non-error output
#' @return bool on whether sample contains only one sample group
#' @rdname detectSingleReplicate
#' @keywords internal
setGeneric(name="detectSingleReplicate", 
           function(nds, quiet) standardGeneric("detectSingleReplicate"))

#' @rdname detectSingleReplicate
setMethod("detectSingleReplicate", "NormalyzerDataset",
          function(nds, quiet=FALSE) {
              
              headerCounts <- table(sampleReplicateGroups(nds))
              nonReplicatedSamples <- names(headerCounts[headerCounts == 1])
              
              if (length(nonReplicatedSamples) > 0) {
                  singleReplicateRun <- TRUE
                  
                  if (!quiet) {
                      message("Non replicated samples in dataset: ",
                              paste(nonReplicatedSamples, collapse=" "),
                              " performing limited single-replicate run\n")
                  }
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
#' @param quiet Only print error messages
#' @return None
#' @rdname detectSingletonSample
#' @keywords internal
setGeneric(name="detectSingletonSample", 
           function(nds, quiet) standardGeneric("detectSingletonSample"))

#' @rdname detectSingletonSample
setMethod("detectSingletonSample", "NormalyzerDataset",
          function(nds, quiet=FALSE) {
              
              groups <- sampleReplicateGroups(nds)
              distinctSamples <- unique(groups)
              singletonSamplePresent <- FALSE
              
              if (length(distinctSamples) == 1) {
                  singletonSamplePresent <- TRUE
                  if (!quiet) {
                      message("Only one replicate group present. ",
                              paste("Group: ", distinctSamples[1]),
                              " Proceeding with limited processing\n")
                  }
              }
              else if (length(distinctSamples) == 0) {
                  stop("No sample replicate groups found - Aborting.\n")
              }
              
              singletonSamplePresent
          }
)

setGeneric(name="checkSingleReplicateRun", 
           function(nds, quiet) standardGeneric("checkSingleReplicateRun"))
setMethod("checkSingleReplicateRun", "NormalyzerDataset",
          function(nds, quiet=FALSE) {
              
              singleReplicateRun <- detectSingleReplicate(nds, quiet=quiet)
              singletonSamplePresent <- detectSingletonSample(nds, quiet=quiet)
              
              if (singleReplicateRun || singletonSamplePresent) {
                  singleReplicateRun <- TRUE
              }
              else {
                  singleReplicateRun <- FALSE
              }
              singleReplicateRun
          })

getRTColumn <- function(annotData, quiet=FALSE) {
    rtColumns <- grep("\\bRT\\b", colnames(annotData), ignore.case = TRUE)
    
    if (length(rtColumns) > 1) {
        errorMessage <- paste(
            "Only able to handle single RT column (name containing RT standing by itself) ",
            "Please change name or remove unwanted RT columns\n")
        stop(errorMessage)
    }
    else if (length(rtColumns) == 1) {
        
        if (!quiet) message("RT annotation column found (", rtColumns, ")")
        
        rtValues <- annotData[, rtColumns]
        retentionTimes <- as.numeric(rtValues)
        return(retentionTimes)
    }
    else {
        if (!quiet) message("No RT column found, skipping RT processing\n")
        return(NULL)
    }
}

