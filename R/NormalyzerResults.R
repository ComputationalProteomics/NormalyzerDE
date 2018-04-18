slotNames <- c("data2log2",
               "data2limloess",
               "fittedLR",
               "data2vsnrep",
               "data2loess",
               "globalfittedRLR",
               "data2vsn",
               "data2GI",
               "data2med",
               "data2mean",
               "data2normfinder",
               "data2quantile",
               "data2rtMed",
               "data2rtMean",
               "data2rtLoess")

outputNames <- c("Log2",
                 "Loess-R",
                 "RLR-R",
                 "VSN-R", 
                 "Loess-G",
                 "RLR-G",
                 "VSN-G",
                 "TI-G",
                 "MedI-G",
                 "AI-G",
                 "NF-G",
                 "Quantile",
                 "RT-mean",
                 "RT-med",
                 "RT-Loess")

#' S4 class to represent dataset information
#' 
#' @slot nds Normalyzer dataset representing run data.
#' @slot ner Normalyzer evaluation results.
#' @slot methodnames Short names for included normalization methods.
#' @slot normfinderMaxThreshold Max threshold for running Normfinder normalization.
#' @slot furtherNormalizationMinThreshold Min threshold for running extended normalizations.
#' @slot houseKeepingVars Vector of house keeping variables found by Normfinder.
#' @slot houseKeepingVarsVals ?
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
#' @slot data2normfinder CTR-log normalized data 
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
                                  normfinderMaxThreshold="numeric",
                                  furtherNormalizationMinThreshold="numeric",
                                  
                                  houseKeepingVars="matrix",
                                  houseKeepingVarsVals="matrix",
                                  
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
                                  data2normfinder = "matrix",
                                  data2quantile = "matrix",
                                  
                                  data2rtMed = "matrix",
                                  data2rtMean = "matrix",
                                  data2rtLoess = "matrix",
                                  
                                  data2ctr = "matrix",
                                  data2mad = "matrix"
                              ),
                              prototype=prototype(nds=NULL, normfinderMaxThreshold=1000, furtherNormalizationMinThreshold=100))

#' Initialize Normalyzer results object
#'
#' @param nr Normalyzer results object.
#' @return None
#' @rdname initializeResultsObject
setGeneric(name="initializeResultsObject", 
           function(nr) standardGeneric("initializeResultsObject"))

#' Main function for executing normalizations
#'
#' @param nr Normalyzer results object.
#' @param forceAll Ignore dataset size limits and run all normalizations
#'  (only meant for testing purposes)
#' @param rtNorm Perform retention time based normalizations
#' @param rtWindow Retention time window size.
#' @param runNormfinder Perform Normfinder normalization.
#' @return None
#' @rdname performNormalizations
setGeneric(name="performNormalizations", 
           function(nr, forceAll, rtNorm, rtWindow, runNormfinder) standardGeneric("performNormalizations"))

#' Generate basic metrics normalizations
#'
#' @param nr Normalyzer results object.
#' @return None
#' @rdname basicMetricNormalizations
setGeneric(name="basicMetricNormalizations", 
           function(nr) standardGeneric("basicMetricNormalizations"))

#' Calculate and assign housekeeping variables
#'
#' @param nr Normalyzer results object.
#' @return None
#' @rdname calculateHKdataForNormObj
setGeneric(name="calculateHKdataForNormObj", 
           function(nr) standardGeneric("calculateHKdataForNormObj"))

#' Generate replicate based normalizations
#'
#' @param nr Normalyzer results object.
#' @return None
#' @rdname performReplicateBasedNormalizations
setGeneric(name="performReplicateBasedNormalizations", 
           function(nr) standardGeneric("performReplicateBasedNormalizations"))

#' Perform retention time normalizations
#'
#' @param nr Results object.
#' @param stepSizeMinutes Size of normalization windows in minutes retention time.
#' @param overlapWindows Number of overlapping normalization windows.
#' @return None
#' @rdname performRTNormalizations
setGeneric(name="performRTNormalizations", 
           function(nr, stepSizeMinutes, overlapWindows) standardGeneric("performRTNormalizations"))

#' Get vector of labels for used methods
#'
#' @param nr Normalyzer results object.
#' @return None
#' @rdname getUsedMethodNames
setGeneric(name="getUsedMethodNames", 
           function(nr) standardGeneric("getUsedMethodNames"))

#' Get list of names for slots
#'
#' @param nr Normalyzer results object.
#' @return None
#' @rdname getSlotNameList
setGeneric(name="getSlotNameList", 
           function(nr) standardGeneric("getSlotNameList"))

#' Get list with normalization matrices
#'
#' @param nr Normalyzer results object.
#' @return None
#' @rdname getNormalizationMatrices
setGeneric(name="getNormalizationMatrices", 
           function(nr) standardGeneric("getNormalizationMatrices"))


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

#' @rdname performNormalizations
setMethod("performNormalizations", "NormalyzerResults",
          function(nr, forceAll=FALSE, rtNorm=FALSE, rtWindow=0.1, runNormfinder=TRUE) {
              
              nds <- nr@nds
              nr <- basicMetricNormalizations(nr)
              rtColPresent <- length(nds@retentionTimes) > 0
              
              if (nrow(nds@normfinderFilterRawData) < nr@normfinderMaxThreshold && runNormfinder || forceAll) {

                  if (!nds@singleReplicateRun) {
                      nr <- normfinder(nr)
                      # nr@houseKeepingVars <- normfinder(nds)
                      nr <- calculateHKdataForNormObj(nr)
                  }
                  else {
                      print("Processing in single replicate mode, Normfinder is omitted")
                  }
              }
              
              if (nrow(nds@filterrawdata) > nr@furtherNormalizationMinThreshold || forceAll) {
                  
                  nr@data2vsn <- performVSNNormalization(nds@filterrawdata)
                  nr@data2quantile <- performQuantileNormalization(nds@filterrawdata)
                  nr@data2mad <- performSMADNormalization(nds@filterrawdata)
                  nr@data2loess <- performCyclicLoessNormalization(nds@filterrawdata)
                  nr@globalfittedRLR <- performGlobalRLRNormalization(nds@filterrawdata)
                  
                  if (!nds@singleReplicateRun && !rtColPresent) {
                      nr <- performReplicateBasedNormalizations(nr)
                  }
                  else if (nds@singleReplicateRun) {
                      print("Processing in single replicate mode, replicate based normalizations are omitted")
                  }
              }
              
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

#' @rdname calculateHKdataForNormObj
setMethod("calculateHKdataForNormObj", "NormalyzerResults",
        function(nr) {
              
            nds <- nr@nds
            filterrawdata <- nds@filterrawdata

            houseKeepVals <- nr@houseKeepingVarsVals
            # HKVarTemp <- as.matrix(nr@houseKeepingVars[, which(as.numeric(nds@sampleReplicateGroups) > 0)])
            # class(HKVarTemp) <- "numeric"
            colmedianctr <- apply(houseKeepVals, 2, FUN="mean")
              
            for(i in 1:nrow(filterrawdata)) {
                nr@data2ctr[i,] <- unlist(sapply(1:ncol(filterrawdata), 
                                                 function(zd) { (filterrawdata[i, zd] / colmedianctr[zd]) * mean(colmedianctr) }))
            }
              
            nr@data2normfinder <- log2(nr@data2ctr)
            colnames(nr@data2normfinder) <- colnames(nds@filterrawdata)
            
            # browser()
            
            nr
        })

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

#' @rdname getUsedMethodNames
setMethod("getUsedMethodNames", "NormalyzerResults",
        function(nr) {
              
            usedMethodNames <- c()
            
            for (i in 1:length(slotNames)) {
                slotName <- slotNames[i]
                fieldValue <- methods::slot(nr, slotName)
                  
                if (!all(is.na(fieldValue))) {
                    outputName <- outputNames[i]
                    usedMethodNames <- c(usedMethodNames, outputName)
                }
            }
              
            usedMethodNames
        })

#' @rdname getSlotNameList
setMethod("getSlotNameList", "NormalyzerResults",
        function(nr) {
            methodDataList <- c()
              
            for (i in 1:length(slotNames)) {
                slotName <- slotNames[i]
                fieldValue <- methods::slot(nr, slotName)
                  
                if (!all(is.na(fieldValue))) {
                    methodDataList <- c(methodDataList, slotName)
                }
            }
              
            methodDataList
        })

#' @rdname getNormalizationMatrices
setMethod("getNormalizationMatrices", "NormalyzerResults",
        function(nr) {
            
            methodDataList <- list()
            listCounter <- 1
            for (i in 1:length(slotNames)) {
                slotName <- slotNames[i]
                fieldValue <- methods::slot(nr, slotName)
                  
                if (!all(is.na(fieldValue))) {
                    methodDataList[[listCounter]] <- fieldValue
                    listCounter <- listCounter + 1
                }
            }
            
            methodDataList
        })



