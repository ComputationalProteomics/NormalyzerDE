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
                                  normalizations = "list"
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
              
              rtColPresent <- length(retentionTimes(nds)) > 0

              normResults <- list()
              
              normResults[["log2"]] <- log2(filterrawdata(nds))
              normResults[["GI"]] <- globalIntensityNormalization(filterrawdata(nds))
              normResults[["median"]] <- medianNormalization(filterrawdata(nds))
              normResults[["mean"]] <- meanNormalization(filterrawdata(nds))
              normResults[["VSN"]] <- performVSNNormalization(filterrawdata(nds))
              normResults[["SMAD"]] <- performSMADNormalization(filterrawdata(nds))
              normResults[["CycLoess"]] <- performCyclicLoessNormalization(filterrawdata(nds))
              normResults[["RLR"]] <- performGlobalRLRNormalization(filterrawdata(nds))
              
              if (rtNorm && rtColPresent) {

                  normResults[["RT-median"]] <- getSmoothedRTNormalizedMatrix(
                      rawMatrix=filterrawdata(nds), 
                      retentionTimes=retentionTimes(nds), 
                      normMethod=medianNormalization, 
                      stepSizeMinutes=rtStepSizeMinutes,
                      windowMinCount=rtWindowMinCount,
                      mergeMethod=rtWindowMergeMethod,
                      windowShifts=rtWindowShifts
                  )
                  
                  normResults[["RT-mean"]] <- getSmoothedRTNormalizedMatrix(
                      rawMatrix=filterrawdata(nds), 
                      retentionTimes=retentionTimes(nds), 
                      normMethod=meanNormalization, 
                      stepSizeMinutes=rtStepSizeMinutes,
                      windowMinCount=rtWindowMinCount,
                      mergeMethod=rtWindowMergeMethod,
                      windowShifts=rtWindowShifts
                  )
                  
                  normResults[["RT-Loess"]] <- getSmoothedRTNormalizedMatrix(
                      rawMatrix=filterrawdata(nds), 
                      retentionTimes=retentionTimes(nds), 
                      normMethod=performCyclicLoessNormalization, 
                      stepSizeMinutes=rtStepSizeMinutes,
                      windowMinCount=rtWindowMinCount,
                      mergeMethod=rtWindowMergeMethod,
                      windowShifts=rtWindowShifts
                  )
                  
                  normResults[["RT-VSN"]] <- getSmoothedRTNormalizedMatrix(
                      rawMatrix=filterrawdata(nds), 
                      retentionTimes=retentionTimes(nds), 
                      normMethod=performVSNNormalization, 
                      stepSizeMinutes=rtStepSizeMinutes,
                      windowMinCount=rtWindowMinCount,
                      mergeMethod=rtWindowMergeMethod,
                      windowShifts=rtWindowShifts
                  )
              }
              else {
                  warning("No RT column specified (column named 'RT') or option not specified"
                          ,"Skipping RT normalization.")
              }
              
              normalizations(nr) <- normResults
              
              nr
          }
)


