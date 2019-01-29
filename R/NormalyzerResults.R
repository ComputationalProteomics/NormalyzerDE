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
#' for running extended normalizations
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @export
NormalyzerResults <- setClass("NormalyzerResults",
                              slots=c(
                                  nds = "NormalyzerDataset",
                                  ner = "NormalyzerEvaluationResults",
                                  normalizations = "list"
                              ))

#' Constructor for NormalyzerResults
#' 
#' @param nds NormalyzerDataset object
#' @return nr Prepared NormalyzerResults object
#' @export
#' @examples
#' data(example_summarized_experiment)
#' normObj <- getVerifiedNormalyzerObject("job_name", example_summarized_experiment)
#' emptyNormResults <- NormalyzerResults(normObj)
NormalyzerResults <- function(nds) {
    nr <- new("NormalyzerResults", nds=nds)
    nr
}

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
#' @param noLogTransform Prevent log-transforming input
#' @param quiet Don't show regular output messages
#' 
#' @return nr NormalyzerDE results object
#' @rdname performNormalizations
#' @keywords internal
setGeneric(name="performNormalizations",
           function(nr, forceAll=FALSE, rtNorm=FALSE, rtStepSizeMinutes=1, 
                    rtWindowMinCount=100, rtWindowShifts=1, 
                    rtWindowMergeMethod="median", noLogTransform=FALSE, quiet=FALSE) 
               standardGeneric("performNormalizations"))
#' @rdname performNormalizations
setMethod("performNormalizations", "NormalyzerResults",
          function(nr, forceAll=FALSE, rtNorm=FALSE, rtStepSizeMinutes=1, 
                   rtWindowMinCount=100, rtWindowShifts=1, 
                   rtWindowMergeMethod="median", noLogTransform=FALSE, quiet=FALSE) {
              
              nds <- nds(nr)
              rtColPresent <- length(retentionTimes(nds)) > 0
              normResults <- list()
              
              if (!noLogTransform) {
                  normResults[["log2"]] <- log2(filterrawdata(nds))
                  
                  if (!isTinyRun(nds)) {
                      normResults[["VSN"]] <- performVSNNormalization(filterrawdata(nds))
                  }
                  else {
                      if (!quiet) {
                          message(
                              "Skipping VSN normalization due to small number of features, ",
                              nrow(filterrawdata(nds)), " features found"
                          )
                      }
                  }
              }
              else {
                  normResults[["log2"]] <- filterrawdata(nds)
                  if (!quiet) {
                      message("VSN normalization assumes non log-transformed data, as the option ",
                              "noLogTransform is specified it is assumed to not need log2-transformation ",
                              "and thus the VSN normalization is skipped\n")
                  }
              }
              
              normResults[["GI"]] <- globalIntensityNormalization(filterrawdata(nds), noLogTransform=noLogTransform)
              normResults[["median"]] <- medianNormalization(filterrawdata(nds), noLogTransform=noLogTransform)
              normResults[["mean"]] <- meanNormalization(filterrawdata(nds), noLogTransform=noLogTransform)
              normResults[["Quantile"]] <- performQuantileNormalization(filterrawdata(nds), noLogTransform=noLogTransform)
              # normResults[["SMAD"]] <- performSMADNormalization(filterrawdata(nds), noLogTransform=noLogTransform)
              normResults[["CycLoess"]] <- performCyclicLoessNormalization(filterrawdata(nds), noLogTransform=noLogTransform)
              normResults[["RLR"]] <- performGlobalRLRNormalization(filterrawdata(nds), noLogTransform=noLogTransform)
              
              enoughDataForRT <- nrow(filterrawdata(nds)) >= rtWindowMinCount
              
              if (!enoughDataForRT) {
                  if (!quiet) {
                      warning("Number of features in data matrix is smaller than minimum normalization window (set by option 'rtWindowMinCount') ",
                              "for RT normalization - skipping RT normalizations.\n")
                  }
              }
              else if (rtNorm && rtColPresent) {

                  normResults[["RT-median"]] <- getSmoothedRTNormalizedMatrix(
                      rawMatrix=filterrawdata(nds), 
                      retentionTimes=retentionTimes(nds), 
                      normMethod=medianNormalization, 
                      stepSizeMinutes=rtStepSizeMinutes,
                      windowMinCount=rtWindowMinCount,
                      mergeMethod=rtWindowMergeMethod,
                      windowShifts=rtWindowShifts, 
                      noLogTransform=noLogTransform
                  )
                  
                  normResults[["RT-mean"]] <- getSmoothedRTNormalizedMatrix(
                      rawMatrix=filterrawdata(nds), 
                      retentionTimes=retentionTimes(nds), 
                      normMethod=meanNormalization, 
                      stepSizeMinutes=rtStepSizeMinutes,
                      windowMinCount=rtWindowMinCount,
                      mergeMethod=rtWindowMergeMethod,
                      windowShifts=rtWindowShifts,
                      noLogTransform=noLogTransform
                  )
                  
                  normResults[["RT-Loess"]] <- getSmoothedRTNormalizedMatrix(
                      rawMatrix=filterrawdata(nds), 
                      retentionTimes=retentionTimes(nds), 
                      normMethod=performCyclicLoessNormalization, 
                      stepSizeMinutes=rtStepSizeMinutes,
                      windowMinCount=rtWindowMinCount,
                      mergeMethod=rtWindowMergeMethod,
                      windowShifts=rtWindowShifts,
                      noLogTransform=noLogTransform
                  )
                  
                  if (!noLogTransform) {
                      normResults[["RT-VSN"]] <- getSmoothedRTNormalizedMatrix(
                          rawMatrix=filterrawdata(nds), 
                          retentionTimes=retentionTimes(nds), 
                          normMethod=performVSNNormalization, 
                          stepSizeMinutes=rtStepSizeMinutes,
                          windowMinCount=rtWindowMinCount,
                          mergeMethod=rtWindowMergeMethod,
                          windowShifts=rtWindowShifts,
                          noLogTransform=noLogTransform
                      )
                  }
                  else {
                      if (!quiet) {
                          message("Skipping RT-VSN, only available for non log-transformed data\n")
                      }
                  }
              }
              else {
                  if (!quiet) {
                      message("No RT column specified (column named 'RT') or option not specified",
                              " Skipping RT normalization.\n")
                  }
              }
              
              normalizations(nr) <- normResults
              
              nr
          }
)


