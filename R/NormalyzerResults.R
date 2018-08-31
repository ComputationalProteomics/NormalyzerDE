slotNamesMap <- c("log2Matrix"="Log2",
                  "loessMatrix"="Loess-G",
                  "vsnMatrix"="VSN-G",
                  "GIMatrix"="AI-G",
                  "medMatrix"="MedI-G",
                  "meanMatrix"="MeanI-G",
                  "quantileMatrix"="Quantile",
                  "rtMedMatrix"="RT-med",
                  "rtMeanMatrix"="RT-mean",
                  "rtLoessMatrix"="RT-Loess",
                  "rtVSNMatrix"="RT-VSN")

#' Representation of the results from performing normalization over a dataset
#' 
#' It is linked to a NormalyzerDataset instance representing the raw data
#' which has been processed. After performing evaluation it also links to
#' an instance of NormalyzerEvaluationResults representing the results
#' from the evaluation.
#' 
#' @slot nds Normalyzer dataset representing run data
#' @slot ner Normalyzer evaluation results
#' @slot methodnames Short names for included normalization methods
#' @slot furtherNormalizationMinThreshold Min threshold for running extended 
#' normalizations
#' @slot log2Matrix Log2 transformed filtered raw data
#' @slot fittedLR Fitted Loess regression 
#' @slot loessMatrix Loess normalization
#' @slot globalfittedRLR Global fitted RLR normalization
#' @slot vsnMatrix Global VSN normalized data
#' @slot GIMatrix GI-normalized data
#' @slot medMatrix Median normalized data
#' @slot meanMatrix Mean normalized data
#' @slot quantileMatrix Quantile normalized data
#' @slot rtMedMatrix Retention time normalized for median
#' @slot rtMeanMatrix Retention time normalized for Mean
#' @slot rtLoessMatrix Retention time normalized for Loess
#' @slot rtVSNMatrix Retention time normalized for VSN
#' @slot ctrMatrix CTR normalized data
#' @slot madMatrix MAD normalized data
#' @export
NormalyzerResults <- setClass("NormalyzerResults",
                              slots=c(
                                  nds = "NormalyzerDataset",
                                  ner = "NormalyzerEvaluationResults",
                                  
                                  methodnames = "character",
                                  furtherNormalizationMinThreshold="numeric",
                                  
                                  normalizations = "SummarizedExperiment",
                                  
                                  log2Matrix = "matrix",
                                  fittedLR = "matrix",
                                  loessMatrix = "matrix",
                                  globalfittedRLR = "matrix",
                                  vsnMatrix = "matrix",
                                  GIMatrix = "matrix",
                                  medMatrix = "matrix",
                                  meanMatrix = "matrix",
                                  quantileMatrix = "matrix",
                                  
                                  rtMedMatrix = "matrix",
                                  rtMeanMatrix = "matrix",
                                  rtLoessMatrix = "matrix",
                                  rtVSNMatrix = "matrix",
                                  
                                  ctrMatrix = "matrix",
                                  madMatrix = "matrix"
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

setGeneric("methodnames", function(object) { standardGeneric("methodnames") })
setMethod("methodnames", signature(object="NormalyzerResults"), 
          function(object) { slot(object, "methodnames") })
setGeneric("methodnames<-", function(object, value) { standardGeneric("methodnames<-") })
setReplaceMethod("methodnames", signature(object="NormalyzerResults"), 
                 function(object, value) { 
                     slot(object, "methodnames") <- value
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

setGeneric("log2Matrix", function(object) { standardGeneric("log2Matrix") })
setMethod("log2Matrix", signature(object="NormalyzerResults"), 
          function(object) { slot(object, "log2Matrix") })
setGeneric("log2Matrix<-", function(object, value) { standardGeneric("log2Matrix<-") })
setReplaceMethod("log2Matrix", signature(object="NormalyzerResults"), 
                 function(object, value) { 
                     slot(object, "log2Matrix") <- value
                     validObject(object)
                     object
                 })

setGeneric("fittedLR", function(object) { standardGeneric("fittedLR") })
setMethod("fittedLR", signature(object="NormalyzerResults"), 
          function(object) { slot(object, "fittedLR") })
setGeneric("fittedLR<-", function(object, value) { standardGeneric("fittedLR<-") })
setReplaceMethod("fittedLR", signature(object="NormalyzerResults"), 
                 function(object, value) { 
                     slot(object, "fittedLR") <- value
                     validObject(object)
                     object
                 })

setGeneric("loessMatrix", function(object) { standardGeneric("loessMatrix") })
setMethod("loessMatrix", signature(object="NormalyzerResults"), 
          function(object) { slot(object, "loessMatrix") })
setGeneric("loessMatrix<-", function(object, value) { standardGeneric("loessMatrix<-") })
setReplaceMethod("loessMatrix", signature(object="NormalyzerResults"), 
                 function(object, value) { 
                     slot(object, "loessMatrix") <- value
                     validObject(object)
                     object
                 })

setGeneric("globalfittedRLR", function(object) { standardGeneric("globalfittedRLR") })
setMethod("globalfittedRLR", signature(object="NormalyzerResults"), 
          function(object) { slot(object, "globalfittedRLR") })
setGeneric("globalfittedRLR<-", function(object, value) { standardGeneric("globalfittedRLR<-") })
setReplaceMethod("globalfittedRLR", signature(object="NormalyzerResults"), 
                 function(object, value) { 
                     slot(object, "globalfittedRLR") <- value
                     validObject(object)
                     object
                 })

setGeneric("vsnMatrix", function(object) { standardGeneric("vsnMatrix") })
setMethod("vsnMatrix", signature(object="NormalyzerResults"), 
          function(object) { slot(object, "vsnMatrix") })
setGeneric("vsnMatrix<-", function(object, value) { standardGeneric("vsnMatrix<-") })
setReplaceMethod("vsnMatrix", signature(object="NormalyzerResults"), 
                 function(object, value) { 
                     slot(object, "vsnMatrix") <- value
                     validObject(object)
                     object
                 })

setGeneric("GIMatrix", function(object) { standardGeneric("GIMatrix") })
setMethod("GIMatrix", signature(object="NormalyzerResults"), 
          function(object) { slot(object, "GIMatrix") })
setGeneric("GIMatrix<-", function(object, value) { standardGeneric("GIMatrix<-") })
setReplaceMethod("GIMatrix", signature(object="NormalyzerResults"), 
                 function(object, value) { 
                     slot(object, "GIMatrix") <- value
                     validObject(object)
                     object
                 })

setGeneric("medMatrix", function(object) { standardGeneric("medMatrix") })
setMethod("medMatrix", signature(object="NormalyzerResults"), 
          function(object) { slot(object, "medMatrix") })
setGeneric("medMatrix<-", function(object, value) { standardGeneric("medMatrix<-") })
setReplaceMethod("medMatrix", signature(object="NormalyzerResults"), 
                 function(object, value) { 
                     slot(object, "medMatrix") <- value
                     validObject(object)
                     object
                 })

setGeneric("meanMatrix", function(object) { standardGeneric("meanMatrix") })
setMethod("meanMatrix", signature(object="NormalyzerResults"), 
          function(object) { slot(object, "meanMatrix") })
setGeneric("meanMatrix<-", function(object, value) { standardGeneric("meanMatrix<-") })
setReplaceMethod("meanMatrix", signature(object="NormalyzerResults"), 
                 function(object, value) { 
                     slot(object, "meanMatrix") <- value
                     validObject(object)
                     object
                 })

setGeneric("quantileMatrix", function(object) { standardGeneric("quantileMatrix") })
setMethod("quantileMatrix", signature(object="NormalyzerResults"), 
          function(object) { slot(object, "quantileMatrix") })
setGeneric("quantileMatrix<-", function(object, value) { standardGeneric("quantileMatrix<-") })
setReplaceMethod("quantileMatrix", signature(object="NormalyzerResults"), 
                 function(object, value) { 
                     slot(object, "quantileMatrix") <- value
                     validObject(object)
                     object
                 })

setGeneric("rtMedMatrix", function(object) { standardGeneric("rtMedMatrix") })
setMethod("rtMedMatrix", signature(object="NormalyzerResults"), 
          function(object) { slot(object, "rtMedMatrix") })
setGeneric("rtMedMatrix<-", function(object, value) { standardGeneric("rtMedMatrix<-") })
setReplaceMethod("rtMedMatrix", signature(object="NormalyzerResults"), 
                 function(object, value) { 
                     slot(object, "rtMedMatrix") <- value
                     validObject(object)
                     object
                 })

setGeneric("rtMeanMatrix", function(object) { standardGeneric("rtMeanMatrix") })
setMethod("rtMeanMatrix", signature(object="NormalyzerResults"), 
          function(object) { slot(object, "rtMeanMatrix") })
setGeneric("rtMeanMatrix<-", function(object, value) { standardGeneric("rtMeanMatrix<-") })
setReplaceMethod("rtMeanMatrix", signature(object="NormalyzerResults"), 
                 function(object, value) { 
                     slot(object, "rtMeanMatrix") <- value
                     validObject(object)
                     object
                 })

setGeneric("rtLoessMatrix", function(object) { standardGeneric("rtLoessMatrix") })
setMethod("rtLoessMatrix", signature(object="NormalyzerResults"), 
          function(object) { slot(object, "rtLoessMatrix") })
setGeneric("rtLoessMatrix<-", function(object, value) { standardGeneric("rtLoessMatrix<-") })
setReplaceMethod("rtLoessMatrix", signature(object="NormalyzerResults"), 
                 function(object, value) { 
                     slot(object, "rtLoessMatrix") <- value
                     validObject(object)
                     object
                 })

setGeneric("rtVSNMatrix", function(object) { standardGeneric("rtVSNMatrix") })
setMethod("rtVSNMatrix", signature(object="NormalyzerResults"), 
          function(object) { slot(object, "rtVSNMatrix") })
setGeneric("rtVSNMatrix<-", function(object, value) { standardGeneric("rtVSNMatrix<-") })
setReplaceMethod("rtVSNMatrix", signature(object="NormalyzerResults"), 
                 function(object, value) { 
                     slot(object, "rtVSNMatrix") <- value
                     validObject(object)
                     object
                 })

setGeneric("ctrMatrix", function(object) { standardGeneric("ctrMatrix") })
setMethod("ctrMatrix", signature(object="NormalyzerResults"), 
          function(object) { slot(object, "ctrMatrix") })
setGeneric("ctrMatrix<-", function(object, value) { standardGeneric("ctrMatrix<-") })
setReplaceMethod("ctrMatrix", signature(object="NormalyzerResults"), 
                 function(object, value) { 
                     slot(object, "ctrMatrix") <- value
                     validObject(object)
                     object
                 })

setGeneric("madMatrix", function(object) { standardGeneric("madMatrix") })
setMethod("madMatrix", signature(object="NormalyzerResults"), 
          function(object) { slot(object, "madMatrix") })
setGeneric("madMatrix<-", function(object, value) { standardGeneric("madMatrix<-") })
setReplaceMethod("madMatrix", signature(object="NormalyzerResults"), 
                 function(object, value) { 
                     slot(object, "madMatrix") <- value
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
              
              # browser()
              
              nds <- nds(nr)
              # log2Matrix(nr) <- log2(filterrawdata(nds))
              
              # browser()
              
              normalizations(nr) <- SummarizedExperiment::SummarizedExperiment(
                  assays=list(log2=log2(filterrawdata(nds)))
              )
              
              # message("Is this matrix used?")
              # ctrMatrix(nr) <- matrix(
              #     nrow=nrow(filterrawdata(nds)),
              #     ncol=ncol(filterrawdata(nds)), byrow=TRUE)
              nr
          }
)


setGeneric(name="addNormalization",
           function(nr, name, normalization) standardGeneric("addNormalization"))
setMethod("addNormalization", "NormalyzerResults",
          function(nr, name, normalization) {
              
              # browser()
              
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
              
              # vsnMatrix(nr) <- performVSNNormalization(filterrawdata(nds))
              # quantileMatrix(nr) <- performQuantileNormalization(filterrawdata(nds))
              # madMatrix(nr) <- performSMADNormalization(filterrawdata(nds))
              # loessMatrix(nr) <- performCyclicLoessNormalization(filterrawdata(nds))
              # globalfittedRLR(nr) <- performGlobalRLRNormalization(filterrawdata(nds))
              
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
              
              # quantileMatrix(nr) <- performQuantileNormalization(filterrawdata(nds))
              # madMatrix(nr) <- performSMADNormalization(filterrawdata(nds))
              # loessMatrix(nr) <- performCyclicLoessNormalization(filterrawdata(nds))
              # globalfittedRLR(nr) <- performGlobalRLRNormalization(filterrawdata(nds))
              
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
                      warning("No RT column specified (column named 'RT'). Skipping RT normalization.")
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
              # GIMatrix(nr) <- globalIntensityNormalization(filterrawdata(nds))
              # medMatrix(nr) <- medianNormalization(filterrawdata(nds))
              # meanMatrix(nr) <- meanNormalization(filterrawdata(nds))
              
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
              
              rtMedMatrix(nr) <- smoothedRTMed
              rtMeanMatrix(nr) <- smoothedRTMean
              rtLoessMatrix(nr) <- smoothedRTLoess
              rtVSNMatrix(nr) <- smoothedRTVSN
              
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
              
              # Returns as a list of matrices rather than SimpleList of
              # sub-SummarizedExperiments
              # lapply(matrices@listData, function(normMat) { normMat })
              names(matrices)
              
              # usedMethodNames <- c()
              # for (i in seq_len(length(slotNamesMap))) {
              #     name <- names(slotNamesMap)[i]
              #     fieldValue <- methods::slot(nr, name)
              # 
              #     if (!all(is.na(fieldValue))) {
              #         outputName <- slotNamesMap[name]
              #         usedMethodNames <- c(usedMethodNames, outputName)
              #     }
              # }
              # 
              # usedMethodNames
          })

#' #' Get list of names for slots
#' #'
#' #' @param nr Normalyzer results object.
#' #' @return None
#' #' @rdname getSlotNameList
#' #' @keywords internal
#' setGeneric(name="getSlotNameList",
#'            function(nr) standardGeneric("getSlotNameList"))
#' 
#' #' @rdname getSlotNameList
#' setMethod("getSlotNameList", "NormalyzerResults",
#'           function(nr) {
#'               methodDataList <- c()
#'               
#'               browser()
#'               
#'               for (i in seq_len(length(slotNamesMap))) {
#'                   slotName <- names(slotNamesMap)[i]
#'                   fieldValue <- methods::slot(nr, slotName)
#'                   
#'                   if (!all(is.na(fieldValue))) {
#'                       methodDataList <- c(methodDataList, slotName)
#'                   }
#'               }
#'               
#'               methodDataList
#'           })

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
              
              # browser()
              
              matrices <- SummarizedExperiment::assays(normalizations(nr))
              
              # Returns as a list of matrices rather than SimpleList of
              # sub-SummarizedExperiments
              # lapply(matrices@listData, function(normMat) { normMat })
              matrices
              
              # methodDataList <- list()
              # listCounter <- 1
              # for (i in seq_len(length(slotNamesMap))) {
              #     slotName <- names(slotNamesMap)[i]
              #     fieldValue <- methods::slot(nr, slotName)
              #     
              #     if (!all(is.na(fieldValue))) {
              #         methodDataList[[listCounter]] <- fieldValue
              #         listCounter <- listCounter + 1
              #     }
              # }
              # methodDataList
          })

