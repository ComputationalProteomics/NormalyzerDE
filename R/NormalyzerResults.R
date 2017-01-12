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
               "data2ctrlog",
               "data2quantile",
               "data2rt")

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
                 "RT-test")

#' S4 class to represent dataset information
#' 
#' @slot nds Normalyzer dataset representing run data.
#' @slot ner Normalyzer evaluation results.
#' @slot methodnames Short names for included normalization methods.
#' @slot normfinderMaxThreshold Max threshold for running Normfinder normalization.
#' @slot furtherNormalizationMinThreshold Min threshold for running extended normalizations.
#' @slot houseKeepingVars Vector of house keeping variables found by Normfinder.
#' @slot data2log2 Log2 of filtered raw data.
#' @slot data2limloess Lim-Loess normalized raw data.
#' @slot fittedLR Fitted Loess regression (TODO: Investigate this further).
#' @slot data2vsnrep Replicate based VSN normalization.
#' @slot data2loess Loess normalization.
#' @slot globalfittedRLR Global fitted RLR normalization (TODO: Investigate further).
#' @slot data2vsn Global VSN normalized data.
#' @slot data2GI GI-normalized data (TODO: Check out).
#' @slot data2med Median normalized data.
#' @slot data2mean Mean normalized data.
#' @slot data2ctrlog CTR-log normalized data (TODO: Check out).
#' @slot data2quantile Quantile normalized data.
#' @slot data2ctr CTR normalized data (TODO: Check out).
#' @slot data2mad MAD normalized data (TODO: Check out).
#' @export
NormalyzerResults <- setClass("NormalyzerResults",
                              slots=c(
                                  nds = "NormalyzerDataset",
                                  ner = "NormalizationEvaluationResults",
                                  
                                  methodnames = "character",
                                  normfinderMaxThreshold="numeric",
                                  furtherNormalizationMinThreshold="numeric",
                                  
                                  houseKeepingVars="matrix",
                                  
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
                                  data2ctrlog = "matrix",
                                  data2quantile = "matrix",
                                  
                                  data2rt = "matrix",
                                  
                                  data2ctr = "matrix",
                                  data2mad = "matrix"
                              ),
                              prototype=prototype(nds=NULL, normfinderMaxThreshold=1000, furtherNormalizationMinThreshold=100))

#' Initialize Normalyzer results object
#'
#' @param nr Normalyzer results object.
#' @return None
#' @rdname initializeResultsObject
#' @export
setGeneric(name="initializeResultsObject", 
           function(nr) standardGeneric("initializeResultsObject"))

#' Main function for executing normalizations
#'
#' @param nr Normalyzer results object.
#' @param forceAll Ignore dataset size limits and run all normalizations
#'  (only meant for testing purposes)
#' @param rtNorm Perform retention time based normalizations
#' @return None
#' @rdname performNormalizations
#' @export
setGeneric(name="performNormalizations", 
           function(nr, forceAll, rtNorm, rtWindow) standardGeneric("performNormalizations"))

#' Generate basic metrics normalizations
#'
#' @param nr Normalyzer results object.
#' @return None
#' @rdname basicMetricNormalizations
#' @export
setGeneric(name="basicMetricNormalizations", 
           function(nr) standardGeneric("basicMetricNormalizations"))

#' Calculate and assign housekeeping variables
#'
#' @param nr Normalyzer results object.
#' @return None
#' @rdname calculateHKdataForNormObj
#' @export
setGeneric(name="calculateHKdataForNormObj", 
           function(nr) standardGeneric("calculateHKdataForNormObj"))

#' Generate replicate based normalizations
#'
#' @param nr Normalyzer results object.
#' @return None
#' @rdname performReplicateBasedNormalizations
#' @export
setGeneric(name="performReplicateBasedNormalizations", 
           function(nr) standardGeneric("performReplicateBasedNormalizations"))

#' Perform retention time normalizations
#'
#' @param nr Normalyzer results object.
#' @param stepSize Normalyzer results object.
#' @return None
#' @rdname performRTNormalizations
#' @export
setGeneric(name="performRTNormalizations", 
           function(nr, stepSizeMinutes) standardGeneric("performRTNormalizations"))

#' Get vector of labels for used methods
#'
#' @param nr Normalyzer results object.
#' @return None
#' @rdname getUsedMethodNames
#' @export
setGeneric(name="getUsedMethodNames", 
           function(nr) standardGeneric("getUsedMethodNames"))

#' Get list of names for slots
#'
#' @param nr Normalyzer results object.
#' @return None
#' @rdname getSlotNameList
#' @export
setGeneric(name="getSlotNameList", 
           function(nr) standardGeneric("getSlotNameList"))

#' Get list with normalization matrices
#'
#' @param nr Normalyzer results object.
#' @return None
#' @rdname getNormalizationMatrices
#' @export
setGeneric(name="getNormalizationMatrices", 
           function(nr) standardGeneric("getNormalizationMatrices"))


#' @rdname initializeResultsObject
setMethod("initializeResultsObject", "NormalyzerResults",
          function(nr) {
              
              nds <- nr@nds
              
              nr@data2log2 <- log2(nds@filterrawdata)
              # nr@data2GI <- matrix(nrow=nrow(nds@filterrawdata),
              #                      ncol=ncol(nds@filterrawdata), byrow=TRUE)
              nr@data2ctr <- matrix(nrow=nrow(nds@filterrawdata),
                                    ncol=ncol(nds@filterrawdata), byrow=TRUE)
              # nr@data2med <- matrix(nrow=nrow(nds@filterrawdata),
              #                       ncol=ncol(nds@filterrawdata), byrow=TRUE)
              # nr@data2mean <- matrix(nrow=nrow(nds@filterrawdata),
              #                        ncol=ncol(nds@filterrawdata), byrow=TRUE)
              
              nr
          }
)

#' @rdname performNormalizations
setMethod("performNormalizations", "NormalyzerResults",
          function(nr, forceAll=FALSE, rtNorm=FALSE, rtWindow=0.1) {
              
              nds <- nr@nds
              nr <- basicMetricNormalizations(nr)
              doRt <- length(nds@retentionTimes) > 0
              
              if (nrow(nds@normfinderFilterRawData) < nr@normfinderMaxThreshold || forceAll) {
                  
                  if (!nds@singleReplicateRun) {
                      nr@houseKeepingVars <- normfinder(nds)
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
                  
                  if (!nds@singleReplicateRun && !doRt) {
                      nr <- performReplicateBasedNormalizations(nr)
                  }
                  else if (nds@singleReplicateRun) {
                      print("Processing in single replicate mode, replicate based normalizations are omitted")
                  }
              }
              
              if (doRt) {
                  nr <- performRTNormalizations(nr, rtWindow)
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

            HKVarTemp <- as.matrix(nr@houseKeepingVars[, which(as.numeric(nds@inputHeaderValues) > 0)])
            class(HKVarTemp) <- "numeric"
            colmedianctr <- apply(HKVarTemp, 2, FUN="mean")
              
            for(i in 1:nrow(filterrawdata)) {
                nr@data2ctr[i,] <- unlist(sapply(1:ncol(filterrawdata), 
                                                 function(zd) { (filterrawdata[i, zd] / colmedianctr[zd]) * mean(colmedianctr) }))
            }
              
            nr@data2ctrlog <- log2(nr@data2ctr)
            colnames(nr@data2ctrlog) <- colnames(nds@filterrawdata)
              
            nr
        })

#' @rdname performReplicateBasedNormalizations
setMethod("performReplicateBasedNormalizations", "NormalyzerResults",
        function(nr) {
            
            ## NORMALIZATION within REPLICATES RLR, VSN and Loess
            
            nds <- nr@nds
              
            sampleReplicateGroups <- nds@sampleReplicateGroups
            filterrawdata <- nds@filterrawdata

            firstIndices <- getFirstIndicesInVector(sampleReplicateGroups, reverse=FALSE)
            lastIndices <- getFirstIndicesInVector(sampleReplicateGroups, reverse=TRUE)
            stopifnot(length(firstIndices) == length(lastIndices))

            limloess <- matrix(, nrow=nrow(filterrawdata), ncol=0)
            vsn <- matrix(, nrow=nrow(filterrawdata), ncol=0)
            fittedLR <- matrix(, nrow=nrow(filterrawdata), ncol=0)
            
            for (sampleIndex in 1:length(firstIndices)) {

                startColIndex <- firstIndices[sampleIndex]
                endColIndex <- lastIndices[sampleIndex]
                
                rawDataWindow <- filterrawdata[, startColIndex:endColIndex]
                
                LRWindow <- performGlobalRLRNormalization(rawDataWindow)
                fittedLR <- cbind(fittedLR, LRWindow)
                
                loessDataWindow <- performCyclicLoessNormalization(rawDataWindow)
                limloess <- cbind(limloess, loessDataWindow)
                
                vsnDataWindow <- performVSNNormalization(rawDataWindow)
                vsn <- cbind(vsn, vsnDataWindow)
            }

            nr@fittedLR <- fittedLR
            nr@data2limloess <- limloess
            nr@data2vsnrep <- vsn
            
            nr
        })

#' @rdname performRTNormalizations
setMethod("performRTNormalizations", "NormalyzerResults",
          function(nr, stepSizeMinutes) {
              
              print(stepSizeMinutes)
              
              nds <- nr@nds
              
              startVal <- min(na.omit(nds@retentionTimes))
              endVal <- max(na.omit(nds@retentionTimes))
              
              rowNumbers <- c()
              isFirstSample <- TRUE
              
              for (windowStart in seq(startVal, endVal, stepSizeMinutes)) {
                  
                  windowEnd <- windowStart + stepSizeMinutes
                  sliceRows <- which(nds@retentionTimes >= windowStart & nds@retentionTimes < windowEnd)
                  rowNumbers <- c(rowNumbers, sliceRows)
                  
                  currentRows <- nds@filterrawdata[sliceRows,]
                  normalizedRow <- medianSortData(currentRows)

                  if (isFirstSample) {
                      dataRows <- normalizedRow
                      isFirstSample <- FALSE
                  }
                  else {
                      dataRows <- rbind(dataRows, normalizedRow)
                  }
              }
              
              print(dataRows)

              # Keep track of original ordering
              orderedDataRows <- dataRows[order(rowNumbers),]
              nr@data2rt <- orderedDataRows
              nr
          })

medianSortData <- function(rawDf) {
    
    medOfData <- apply(rawDf, 2, FUN="median", na.rm=TRUE)
    medSortDf <- matrix(ncol=ncol(rawDf), nrow=nrow(rawDf))
    
    for (i in 1:nrow(rawDf)) {
        
        medSortDf[i,] <- unlist(sapply(1:ncol(rawDf), 
                                       function(zd) { (rawDf[i, zd] / medOfData[zd]) * mean(na.omit(medOfData)) }))
            
    }
    
    medSortDf
}

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



