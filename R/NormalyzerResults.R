slotNamesMap <- c("data2log2"="Log2",
               "data2loess"="Loess-G",
               "data2vsn"="VSN-G",
               "data2GI"="AI-G",
               "data2med"="MedI-G",
               "data2mean"="MeanI-G",
               "data2quantile"="Quantile",
               "data2rtMed"="RT-med",
               "data2rtMean"="RT-mean",
               "data2rtLoess"="RT-Loess")

# slotOutput <- c("data2log2"="Log2",
#                 "data2loess"="Loess-G",
#                 "data2vsn"="VSN-G",
#                 "data2rlr"="RLR-G",
#                 "data2GI"="AI-G",
#                 "data2med"="MedI-G",
#                 "data2mean"="MeanI-G",
#                 "data2quantile"="Quantile",
#                 "data2rtMed"="RT-med",
#                 "data2rtMean"="RT-mean",
#                 "data2rtLoess"="RT-Loess")

# outputNames <- c("Log2",
#                  "Loess-G",
#                  "VSN-G",
#                  "RLR-G",
#                  "TI-G",
#                  "MedI-G",
#                  "AI-G",
#                  "Quantile",
#                  "RT-mean",
#                  "RT-med",
#                  "RT-Loess")

#' S4 class to represent dataset information
#' 
#' @slot nds Normalyzer dataset representing run data.
#' @slot ner Normalyzer evaluation results.
#' @slot methodnames Short names for included normalization methods.
#' @slot furtherNormalizationMinThreshold Min threshold for running extended normalizations.
#' @slot data2log2 Log2 of filtered raw data.
#' @slot data2limloess Lim-Loess normalized raw data.
#' @slot fittedLR Fitted Loess regression 
#' @slot data2vsnrep Replicate based VSN normalization.
#' @slot data2loess Loess normalization.
#' @slot globalfittedRLR Global fitted RLR normalization
#' @slot data2vsn Global VSN normalized data.
#' @slot data2GI GI-normalized data
#' @slot data2med Median normalized data.
#' @slot data2mean Mean normalized data.
#' @slot data2quantile Quantile normalized data.
#' @slot data2rtMed Retention time normalized for median.
#' @slot data2rtMean Retention time normalized for Mean.
#' @slot data2rtLoess Retention time normalized for Loess.
#' @slot data2ctr CTR normalized data
#' @slot data2mad MAD normalized data
#' @export
NormalyzerResults <- setClass("NormalyzerResults",
                              slots=c(
                                  nds = "NormalyzerDataset",
                                  ner = "NormalizationEvaluationResults",

                                  methodnames = "character",
                                  furtherNormalizationMinThreshold="numeric",

                                  data2log2 = "matrix",
                                  data2limloess = "matrix",
                                  fittedLR = "matrix",
                                  data2vsnrep = "matrix",
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

                                  data2ctr = "matrix",
                                  data2mad = "matrix"
                              ),
                              prototype=prototype(nds=NULL, furtherNormalizationMinThreshold=100))

#' Initialize Normalyzer results object
#'
#' @param nr Normalyzer results object.
#' @return None
#' @rdname initializeResultsObject
setGeneric(name="initializeResultsObject",
           function(nr) standardGeneric("initializeResultsObject"))

#' @rdname initializeResultsObject
setMethod("initializeResultsObject", "NormalyzerResults",
          function(nr) {

              nds <- nr@nds
              nr@data2log2 <- log2(nds@filterrawdata)
              nr@data2ctr <- matrix(nrow=nrow(nds@filterrawdata),
                                    ncol=ncol(nds@filterrawdata), byrow=TRUE)
              nr
          }
)

#' Main function for executing normalizations
#'
#' @param nr Normalyzer results object.
#' @param forceAll Ignore dataset size limits and run all normalizations
#'  (only meant for testing purposes)
#' @param rtNorm Perform retention time based normalizations
#' @param rtWindow Retention time window size.
#' @return None
#' @rdname performNormalizations
setGeneric(name="performNormalizations",
           function(nr, forceAll, rtNorm, rtWindow) standardGeneric("performNormalizations"))

#' @rdname performNormalizations
setMethod("performNormalizations", "NormalyzerResults",
          function(nr, forceAll=FALSE, rtNorm=FALSE, rtWindow=0.1) {

              nds <- nr@nds
              nr <- basicMetricNormalizations(nr)
              rtColPresent <- length(nds@retentionTimes) > 0

                nr@data2vsn <- performVSNNormalization(nds@filterrawdata)
                nr@data2quantile <- performQuantileNormalization(nds@filterrawdata)
                nr@data2mad <- performSMADNormalization(nds@filterrawdata)
                nr@data2loess <- performCyclicLoessNormalization(nds@filterrawdata)
                nr@globalfittedRLR <- performGlobalRLRNormalization(nds@filterrawdata)

                  # if (!nds@singleReplicateRun && !rtColPresent) {
                  #     nr <- performReplicateBasedNormalizations(nr)
                  # }
                  # else if (nds@singleReplicateRun) {
                  #     print("Processing in single replicate mode, replicate based normalizations are omitted")
                  # }
              # }

              if (rtNorm) {
                  if (rtColPresent) {
                      nr <- performRTNormalizations(nr, rtWindow)
                  }
                  else {
                      print("No RT column specified (column named 'RT'). Skipping RT normalization.")
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
setGeneric(name="basicMetricNormalizations",
           function(nr) standardGeneric("basicMetricNormalizations"))

#' @rdname basicMetricNormalizations
setMethod("basicMetricNormalizations", "NormalyzerResults",
          function(nr) {

              nds <- nr@nds
              nr@data2GI <- globalIntensityNormalization(nds@filterrawdata)
              nr@data2med <- medianNormalization(nds@filterrawdata)
              nr@data2mean <- meanNormalization(nds@filterrawdata)
              nr
          }
)


#' Generate replicate based normalizations
#'
#' @param nr Normalyzer results object.
#' @return None
#' @rdname performReplicateBasedNormalizations
setGeneric(name="performReplicateBasedNormalizations",
           function(nr) standardGeneric("performReplicateBasedNormalizations"))

#' @rdname performReplicateBasedNormalizations
setMethod("performReplicateBasedNormalizations", "NormalyzerResults",
          function(nr) {

              ## NORMALIZATION within REPLICATES RLR, VSN and Loess

              nds <- nr@nds

              sampleReplicateGroups <- nds@sampleReplicateGroups
              filterrawdata <- nds@filterrawdata

              limLoessMatrix <- matrix(, nrow=nrow(filterrawdata), ncol=0)
              vsnMatrix <- matrix(, nrow=nrow(filterrawdata), ncol=0)
              fittedLRMatrix <- matrix(, nrow=nrow(filterrawdata), ncol=0)

              indexList <- getIndexList(sampleReplicateGroups)
              colOrder <- c()

              for (repVal in names(indexList)) {

                  repValNr <- strtoi(repVal)
                  cols <- indexList[[repVal]]
                  colOrder <- c(colOrder, cols)

                  rawDataWindow <- filterrawdata[, cols]
                  LRWindow <- performGlobalRLRNormalization(rawDataWindow)
                  fittedLRMatrix <- cbind(fittedLRMatrix, LRWindow)

                  loessDataWindow <- performCyclicLoessNormalization(rawDataWindow)
                  limLoessMatrix <- cbind(limLoessMatrix, loessDataWindow)

                  vsnDataWindow <- performVSNNormalization(rawDataWindow)
                  vsnMatrix <- cbind(vsnMatrix, vsnDataWindow)
              }

              origPos <- c()
              for (i in 1:length(colOrder)) {
                  origPos <- c(origPos, which(colOrder == i))
              }

              nr@fittedLR <- fittedLRMatrix[, origPos]
              nr@data2limloess <- limLoessMatrix[, origPos]
              nr@data2vsnrep <- vsnMatrix[, origPos]

              nr
          })

#' Perform retention time normalizations
#'
#' @param nr Results object.
#' @param stepSizeMinutes Size of normalization windows in minutes retention time.
#' @param overlapWindows Number of overlapping normalization windows.
#' @return None
#' @rdname performRTNormalizations
setGeneric(name="performRTNormalizations",
           function(nr, stepSizeMinutes, overlapWindows) standardGeneric("performRTNormalizations"))

#' @rdname performRTNormalizations
setMethod("performRTNormalizations", "NormalyzerResults",
          function(nr, stepSizeMinutes, overlapWindows=FALSE) {

              nds <- nr@nds

              smoothedRTMed <- getSmoothedRTNormalizedMatrix(nds@filterrawdata, nds@retentionTimes, medianNormalization, stepSizeMinutes)
              smoothedRTMean <- getSmoothedRTNormalizedMatrix(nds@filterrawdata, nds@retentionTimes, meanNormalization, stepSizeMinutes)
              smoothedRTLoess <- getSmoothedRTNormalizedMatrix(nds@filterrawdata, nds@retentionTimes, performCyclicLoessNormalization, stepSizeMinutes)

              nr@data2rtMed <- smoothedRTMed
              nr@data2rtMean <- smoothedRTMean
              nr@data2rtLoess <- smoothedRTLoess
              nr
          })

#' Get vector of labels for used methods
#'
#' @param nr Normalyzer results object.
#' @return None
#' @rdname getUsedMethodNames
setGeneric(name="getUsedMethodNames",
           function(nr) standardGeneric("getUsedMethodNames"))

#' @rdname getUsedMethodNames
setMethod("getUsedMethodNames", "NormalyzerResults",
          function(nr) {

              usedMethodNames <- c()
              for (i in 1:length(slotNamesMap)) {
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
setGeneric(name="getSlotNameList",
           function(nr) standardGeneric("getSlotNameList"))

#' @rdname getSlotNameList
setMethod("getSlotNameList", "NormalyzerResults",
          function(nr) {
              methodDataList <- c()

              for (i in 1:length(slotNamesMap)) {
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
setGeneric(name="getNormalizationMatrices",
           function(nr) standardGeneric("getNormalizationMatrices"))

#' @rdname getNormalizationMatrices
setMethod("getNormalizationMatrices", "NormalyzerResults",
          function(nr) {

              methodDataList <- list()
              listCounter <- 1
              for (i in 1:length(slotNamesMap)) {
                  slotName <- names(slotNamesMap)[i]
                  fieldValue <- methods::slot(nr, slotName)

                  if (!all(is.na(fieldValue))) {
                      methodDataList[[listCounter]] <- fieldValue
                      listCounter <- listCounter + 1
                  }
                  
                  # if (!all(is.na(fieldValue))) {
                  #   methodDataList[[listCounter]] <- fieldValue
                  #   listCounter <- listCounter + 1
                  # }
              }

              methodDataList
          })









