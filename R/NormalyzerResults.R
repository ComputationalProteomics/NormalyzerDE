#' Representation of the results from performing normalization over a dataset
#' 
#' It is linked to a NormalyzerDataset instance representing the raw data
#' which has been processed. After performing evaluation it also links to
#' an instance of NormalyzerEvaluationResults representing the results
#' from the evaluation.
#' 
#' @slot normalizations SummarizedExperiment object containing calculated
#'   normalization results
#' @slot nds Normalyzer dataset representing run data
#' @slot ner Normalyzer evaluation results
#' @slot furtherNormalizationMinThreshold Minimum number of features threshold 
#' for running extended normalizations
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @export
NormalyzerResults <- setClass("NormalyzerResults",
                              slots=c(
                                  nds = "NormalyzerDataset",
                                  ner = "NormalyzerEvaluationResults",
                                  furtherNormalizationMinThreshold="numeric",
                                  normalizations = "SummarizedExperiment"
                              ),
                              prototype=prototype(
                                  nds=NULL, 
                                  furtherNormalizationMinThreshold=100))


setGeneric("normalizations", function(object) { standardGeneric("normalizations") })
setMethod("normalizations", signature(object="NormalyzerResults"), 
          function(object) { slot(object, "normalizations") })
setGeneric("normalizations<-", function(object, value) { standardGeneric("normalizations<-") })
setReplaceMethod("normalizations", signature(object="NormalyzerResults"), 
                 function(object, value) { 
                     slot(object, "normalizations") <- value
                     validObject(object)
                     object
                 })

setGeneric("nds", function(object) { standardGeneric("nds") })
setMethod("nds", signature(object="NormalyzerResults"), 
          function(object) { slot(object, "nds") })
setGeneric("nds<-", function(object, value) { standardGeneric("nds<-") })
setReplaceMethod("nds", signature(object="NormalyzerResults"), 
                 function(object, value) { 
                     slot(object, "nds") <- value
                     validObject(object)
                     object
                 })

setGeneric("ner", function(object) { standardGeneric("ner") })
setMethod("ner", signature(object="NormalyzerResults"), 
          function(object) { slot(object, "ner") })
setGeneric("ner<-", function(object, value) { standardGeneric("ner<-") })
setReplaceMethod("ner", signature(object="NormalyzerResults"), 
                 function(object, value) { 
                     slot(object, "ner") <- value
                     validObject(object)
                     object
                 })

setGeneric("furtherNormalizationMinThreshold", function(object) { standardGeneric("furtherNormalizationMinThreshold") })
setMethod("furtherNormalizationMinThreshold", signature(object="NormalyzerResults"), 
          function(object) { slot(object, "furtherNormalizationMinThreshold") })
setGeneric("furtherNormalizationMinThreshold<-", function(object, value) { standardGeneric("furtherNormalizationMinThreshold<-") })
setReplaceMethod("furtherNormalizationMinThreshold", signature(object="NormalyzerResults"), 
                 function(object, value) { 
                     slot(object, "furtherNormalizationMinThreshold") <- value
                     validObject(object)
                     object
                 })

#' Initialize Normalyzer results object
#'
#' @param nr Normalyzer results object.
#' @return None
#' @rdname initializeResultsObject
#' @keywords internal
setGeneric(name="initializeResultsObject",
           function(nr) standardGeneric("initializeResultsObject"))

#' @rdname initializeResultsObject
setMethod("initializeResultsObject", "NormalyzerResults",
          function(nr) {
              
              nds <- nds(nr)
              normalizations(nr) <- SummarizedExperiment::SummarizedExperiment(
                  assays=list(log2=log2(filterrawdata(nds)))
              )
              
              nr
          }
)

setGeneric(name="addNormalization",
           function(nr, name, normalization) standardGeneric("addNormalization"))
setMethod("addNormalization", "NormalyzerResults",
          function(nr, name, normalization) {
              
              normObj <- normalizations(nr)
              SummarizedExperiment::assays(normObj)[[name]] <- normalization
              normalizations(nr) <- normObj
              nr
          }
)

#' Main function for executing normalizations
#'
#' @param nr Normalyzer results object.
#' @param forceAll Ignore dataset size limits and run all normalizations
#'  (only meant for testing purposes)
#' @param rtNorm Perform retention time based normalizations
#' 
#' @param rtStepSizeMinutes Retention time normalization window size.
#' @param rtWindowMinCount Minimum number of datapoints in each retention-time
#'   segment.
#' @param rtWindowShifts Number of layered retention time normalized windows.
#' @param rtWindowMergeMethod Merge approach for layered retention time windows.
#' 
#' @return nr NormalyzerDE results object
#' @rdname performNormalizations
#' @keywords internal
setGeneric(name="performNormalizations",
           function(nr, forceAll, rtNorm, rtStepSizeMinutes, rtWindowMinCount, 
                    rtWindowShifts, rtWindowMergeMethod, quiet) 
               standardGeneric("performNormalizations"))
#' @rdname performNormalizations
setMethod("performNormalizations", "NormalyzerResults",
          function(nr, forceAll=FALSE, rtNorm=FALSE, rtStepSizeMinutes=1, 
                   rtWindowMinCount=100, rtWindowShifts=1, 
                   rtWindowMergeMethod="median") {
              
              nds <- nds(nr)
              nr <- basicMetricNormalizations(nr)
              rtColPresent <- length(retentionTimes(nds)) > 0
              
              nr <- addNormalization(
                  nr, 
                  "VSN", 
                  performVSNNormalization(filterrawdata(nds)))
              nr <- addNormalization(
                  nr, 
                  "SMAD", 
                  performSMADNormalization(filterrawdata(nds)))
              nr <- addNormalization(
                  nr, 
                  "CycLoess", 
                  performCyclicLoessNormalization(filterrawdata(nds)))
              nr <- addNormalization(
                  nr, 
                  "RLR", 
                  performGlobalRLRNormalization(filterrawdata(nds)))
              
              if (rtNorm) {
                  if (rtColPresent) {
                      nr <- performRTNormalizations(
                          nr, 
                          stepSizeMinutes=rtStepSizeMinutes,
                          minWindowSize=rtWindowMinCount,
                          windowShifts=rtWindowShifts,
                          mergeMethod=rtWindowMergeMethod)
                  }
                  else {
                      warning("No RT column specified (column named 'RT')." 
                              ,"Skipping RT normalization.")
                  }
              }
              nr
          }
)


#' Generate basic metrics normalizations
#'
#' @param nr Normalyzer results object.
#' @return None
#' @rdname basicMetricNormalizations
#' @keywords internal
setGeneric(name="basicMetricNormalizations",
           function(nr) standardGeneric("basicMetricNormalizations"))

#' @rdname basicMetricNormalizations
setMethod("basicMetricNormalizations", "NormalyzerResults",
          function(nr) {
              
              nds <- nds(nr)

              nr <- addNormalization(
                  nr, 
                  "GI", 
                  globalIntensityNormalization(filterrawdata(nds)))
              
              nr <- addNormalization(
                  nr, 
                  "median", 
                  medianNormalization(filterrawdata(nds)))
              
              nr <- addNormalization(
                  nr, 
                  "mean", 
                  meanNormalization(filterrawdata(nds)))
              
              nr
          }
)

#' Perform retention time normalizations
#'
#' @param nr Results object.
#' @param stepSizeMinutes Size of normalization windows in minutes retention time.
#' @param minWindowSize Minimum number of datapoints within each retention time window.
#' @param windowShifts Number of overlapping retention-time segmented normalizations.
#' @param mergeMethod Type of merge approach for overlapping windows.
#' @return nr NormalyzerDE results object
#' @rdname performRTNormalizations
#' @keywords internal
setGeneric(name="performRTNormalizations",
           function(nr, stepSizeMinutes, minWindowSize, windowShifts, mergeMethod) standardGeneric("performRTNormalizations"))

#' @rdname performRTNormalizations
setMethod("performRTNormalizations", "NormalyzerResults",
          function(nr, stepSizeMinutes, minWindowSize, windowShifts, mergeMethod) {
              
              nds <- nds(nr)
              
              smoothedRTMed <- getSmoothedRTNormalizedMatrix(
                  rawMatrix=filterrawdata(nds), 
                  retentionTimes=retentionTimes(nds), 
                  normMethod=medianNormalization, 
                  stepSizeMinutes=stepSizeMinutes,
                  windowMinCount=minWindowSize,
                  mergeMethod=mergeMethod,
                  windowShifts=windowShifts
              )
              
              smoothedRTMean <- getSmoothedRTNormalizedMatrix(
                  rawMatrix=filterrawdata(nds), 
                  retentionTimes=retentionTimes(nds), 
                  normMethod=meanNormalization, 
                  stepSizeMinutes=stepSizeMinutes,
                  windowMinCount=minWindowSize,
                  mergeMethod=mergeMethod,
                  windowShifts=windowShifts
              )
              
              smoothedRTLoess <- getSmoothedRTNormalizedMatrix(
                  rawMatrix=filterrawdata(nds), 
                  retentionTimes=retentionTimes(nds), 
                  normMethod=performCyclicLoessNormalization, 
                  stepSizeMinutes=stepSizeMinutes,
                  windowMinCount=minWindowSize,
                  mergeMethod=mergeMethod,
                  windowShifts=windowShifts
              )
              
              smoothedRTVSN <- getSmoothedRTNormalizedMatrix(
                  rawMatrix=filterrawdata(nds), 
                  retentionTimes=retentionTimes(nds), 
                  normMethod=performVSNNormalization, 
                  stepSizeMinutes=stepSizeMinutes,
                  windowMinCount=minWindowSize,
                  mergeMethod=mergeMethod,
                  windowShifts=windowShifts
              )
              
              nr <- addNormalization(
                  nr, 
                  "RT-median", 
                  smoothedRTMed)
              
              nr <- addNormalization(
                  nr, 
                  "RT-mean", 
                  smoothedRTMean)
              
              nr <- addNormalization(
                  nr, 
                  "RT-Loess", 
                  smoothedRTLoess)
              
              nr <- addNormalization(
                  nr, 
                  "RT-VSN", 
                  smoothedRTVSN)
              
              nr
          })

#' Get vector of labels for used methods
#'
#' @param nr Normalyzer results object.
#' @return None
#' @rdname getUsedMethodNames
#' @keywords internal
setGeneric(name="getUsedMethodNames",
           function(nr) standardGeneric("getUsedMethodNames"))

#' @rdname getUsedMethodNames
setMethod("getUsedMethodNames", "NormalyzerResults",
          function(nr) {

              matrices <- SummarizedExperiment::assays(normalizations(nr))
              names(matrices)
          })


#' Get list with normalization matrices
#'
#' @param nr Normalyzer results object.
#' @return None
#' @rdname getNormalizationMatrices
#' @keywords internal
setGeneric(name="getNormalizationMatrices",
           function(nr) standardGeneric("getNormalizationMatrices"))

#' @rdname getNormalizationMatrices
setMethod("getNormalizationMatrices", "NormalyzerResults",
          function(nr) {
              
              matrices <- SummarizedExperiment::assays(normalizations(nr))
              matrices
          })

