slotNamesMap <- c("data2log2"="Log2",
                  "data2loess"="Loess-G",
                  "data2vsn"="VSN-G",
                  "data2GI"="AI-G",
                  "data2med"="MedI-G",
                  "data2mean"="MeanI-G",
                  "data2quantile"="Quantile",
                  "data2rtMed"="RT-med",
                  "data2rtMean"="RT-mean",
                  "data2rtLoess"="RT-Loess",
                  "data2rtVSN"="RT-VSN")

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
#' @slot data2log2 Log2 of filtered raw data
#' @slot fittedLR Fitted Loess regression 
#' @slot data2loess Loess normalization
#' @slot globalfittedRLR Global fitted RLR normalization
#' @slot data2vsn Global VSN normalized data
#' @slot data2GI GI-normalized data
#' @slot data2med Median normalized data
#' @slot data2mean Mean normalized data
#' @slot data2quantile Quantile normalized data
#' @slot data2rtMed Retention time normalized for median
#' @slot data2rtMean Retention time normalized for Mean
#' @slot data2rtLoess Retention time normalized for Loess
#' @slot data2rtVSN Retention time normalized for VSN
#' @slot data2ctr CTR normalized data
#' @slot data2mad MAD normalized data
#' @export
NormalyzerResults <- setClass("NormalyzerResults",
                              slots=c(
                                  nds = "NormalyzerDataset",
                                  ner = "NormalyzerEvaluationResults",
                                  
                                  methodnames = "character",
                                  furtherNormalizationMinThreshold="numeric",
                                  
                                  data2log2 = "matrix",
                                  fittedLR = "matrix",
                                  data2loess = "matrix",
                                  globalfittedRLR = "matrix",
                                  data2vsn = "matrix",
                                  data2GI = "matrix",
                                  data2med = "matrix",
                                  data2mean = "matrix",
                                  data2quantile = "matrix",
                                  
                                  data2rtMed = "matrix",
                                  data2rtMean = "matrix",
                                  data2rtLoess = "matrix",
                                  data2rtVSN = "matrix",
                                  
                                  data2ctr = "matrix",
                                  data2mad = "matrix"
                              ),
                              prototype=prototype(
                                  nds=NULL, 
                                  furtherNormalizationMinThreshold=100))

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
              
              nds <- nr@nds
              nr@data2log2 <- log2(filterrawdata(nds))
              nr@data2ctr <- matrix(nrow=nrow(filterrawdata(nds)),
                                    ncol=ncol(filterrawdata(nds)), byrow=TRUE)
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
              
              nds <- nr@nds
              nr <- basicMetricNormalizations(nr)
              rtColPresent <- length(retentionTimes(nds)) > 0
              
              nr@data2vsn <- performVSNNormalization(filterrawdata(nds))
              nr@data2quantile <- performQuantileNormalization(filterrawdata(nds))
              nr@data2mad <- performSMADNormalization(filterrawdata(nds))
              nr@data2loess <- performCyclicLoessNormalization(filterrawdata(nds))
              nr@globalfittedRLR <- performGlobalRLRNormalization(filterrawdata(nds))
              
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
              
              nds <- nr@nds
              nr@data2GI <- globalIntensityNormalization(filterrawdata(nds))
              nr@data2med <- medianNormalization(filterrawdata(nds))
              nr@data2mean <- meanNormalization(filterrawdata(nds))
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
              
              nds <- nr@nds
              
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
              
              nr@data2rtMed <- smoothedRTMed
              nr@data2rtMean <- smoothedRTMean
              nr@data2rtLoess <- smoothedRTLoess
              nr@data2rtVSN <- smoothedRTVSN
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
              
              usedMethodNames <- c()
              for (i in seq_len(length(slotNamesMap))) {
                  name <- names(slotNamesMap)[i]
                  fieldValue <- methods::slot(nr, name)
                  
                  if (!all(is.na(fieldValue))) {
                      outputName <- slotNamesMap[name]
                      usedMethodNames <- c(usedMethodNames, outputName)
                  }
              }
              
              usedMethodNames
          })

#' Get list of names for slots
#'
#' @param nr Normalyzer results object.
#' @return None
#' @rdname getSlotNameList
#' @keywords internal
setGeneric(name="getSlotNameList",
           function(nr) standardGeneric("getSlotNameList"))

#' @rdname getSlotNameList
setMethod("getSlotNameList", "NormalyzerResults",
          function(nr) {
              methodDataList <- c()
              
              for (i in seq_len(length(slotNamesMap))) {
                  slotName <- names(slotNamesMap)[i]
                  fieldValue <- methods::slot(nr, slotName)
                  
                  if (!all(is.na(fieldValue))) {
                      methodDataList <- c(methodDataList, slotName)
                  }
              }
              
              methodDataList
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
              
              methodDataList <- list()
              listCounter <- 1
              for (i in seq_len(length(slotNamesMap))) {
                  slotName <- names(slotNamesMap)[i]
                  fieldValue <- methods::slot(nr, slotName)
                  
                  if (!all(is.na(fieldValue))) {
                      methodDataList[[listCounter]] <- fieldValue
                      listCounter <- listCounter + 1
                  }
              }
              methodDataList
          })









