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
               "data2quantile")

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
                 "Quantile")

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
#' @return None
#' @rdname performNormalizations
#' @export
setGeneric(name="performNormalizations", 
           function(nr) standardGeneric("performNormalizations"))

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

#' Generate VSN normalization
#'
#' @param nr Normalyzer results object.
#' @return None
#' @rdname performVSNNormalization
#' @export
setGeneric(name="performVSNNormalization", 
           function(nr) standardGeneric("performVSNNormalization"))

#' Generate quantile normalization
#'
#' @param nr Normalyzer results object.
#' @return None
#' @rdname performQuantileNormalization
#' @export
setGeneric(name="performQuantileNormalization", 
           function(nr) standardGeneric("performQuantileNormalization"))

#' Generate SMAD normalization
#'
#' @param nr Normalyzer results object.
#' @return None
#' @rdname performSMADNormalization
#' @export
setGeneric(name="performSMADNormalization", 
           function(nr) standardGeneric("performSMADNormalization"))

#' Generate cyclic Loess normalization
#'
#' @param nr Normalyzer results object.
#' @return None
#' @rdname performCyclicLoessNormalization
#' @export
setGeneric(name="performCyclicLoessNormalization", 
           function(nr) standardGeneric("performCyclicLoessNormalization"))

#' Generate global RLR normalization
#'
#' @param nr Normalyzer results object.
#' @return None
#' @rdname performGlobalRLRNormalization
#' @export
setGeneric(name="performGlobalRLRNormalization", 
           function(nr) standardGeneric("performGlobalRLRNormalization"))

#' Generate replicate based normalizations
#'
#' @param nr Normalyzer results object.
#' @return None
#' @rdname performReplicateBasedNormalizations
#' @export
setGeneric(name="performReplicateBasedNormalizations", 
           function(nr) standardGeneric("performReplicateBasedNormalizations"))

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
              nr@data2GI <- matrix(nrow=nrow(nds@filterrawdata), 
                                   ncol=ncol(nds@filterrawdata), byrow=TRUE)
              nr@data2ctr <- matrix(nrow=nrow(nds@filterrawdata), 
                                    ncol=ncol(nds@filterrawdata), byrow=TRUE)
              nr@data2med <- matrix(nrow=nrow(nds@filterrawdata), 
                                    ncol=ncol(nds@filterrawdata), byrow=TRUE) 
              nr@data2mean <- matrix(nrow=nrow(nds@filterrawdata), 
                                     ncol=ncol(nds@filterrawdata), byrow=TRUE)
              
              nr
          }
)

#' @rdname performNormalizations
setMethod("performNormalizations", "NormalyzerResults",
          function(nr) {
              
              nds <- nr@nds
              nr <- basicMetricNormalizations(nr)
              
              if (nrow(nds@normfinderFilterRawData) < nr@normfinderMaxThreshold) {
                  nr@houseKeepingVars <- normfinder(nds)
                  nr <- calculateHKdataForNormObj(nr)
              }
              
              if (nrow(nds@filterrawdata) > nr@furtherNormalizationMinThreshold) {
                  nr <- performVSNNormalization(nr)
                  nr <- performQuantileNormalization(nr)
                  nr <- performSMADNormalization(nr)
                  nr <- performCyclicLoessNormalization(nr)
                  nr <- performGlobalRLRNormalization(nr)
                  
                  ## NORMALIZATION within REPLICATES RLR, VSN and Loess
                  nr <- performReplicateBasedNormalizations(nr)
              }
              
              nr
          }
)

#' @rdname basicMetricNormalizations
setMethod("basicMetricNormalizations", "NormalyzerResults",
          function(nr) {
            ## Normalizing using average sample sum, sample mean and sample median
              
            nds <- nr@nds

            avgcolsum <- stats::median(nds@colsum)
            for (i in 1:nrow(nds@filterrawdata)) {
                nr@data2GI[i, ] <- unlist(sapply(1:ncol(nds@filterrawdata), 
                                                 function(zd) { (nds@filterrawdata[i, zd] / nds@colsum[zd]) * (avgcolsum) }))
                nr@data2med[i, ] <- unlist(sapply(1:ncol(nds@filterrawdata), 
                                                  function(zd) { (nds@filterrawdata[i, zd] / nds@medofdata[zd]) * mean(nds@medofdata) }))
                nr@data2mean[i, ] <- unlist(sapply(1:ncol(nds@filterrawdata), 
                                                   function(zd) { (nds@filterrawdata[i, zd] / nds@meanofdata[zd]) * mean(nds@meanofdata) }))
            } 
              
            nr@data2GI <- log2(nr@data2GI)
            nr@data2med <- log2(nr@data2med)
            nr@data2mean <- log2(nr@data2mean)
            colnames(nr@data2GI) <- colnames(nr@data2log2)
            colnames(nr@data2med) <- colnames(nr@data2log2)
            colnames(nr@data2mean) <- colnames(nr@data2log2)
              
            nr
        }
)

#' @rdname calculateHKdataForNormObj
setMethod("calculateHKdataForNormObj", "NormalyzerResults",
        function(nr) {
              
            nds <- nr@nds
            filterrawdata <- nds@filterrawdata
              
            Hkvartemp <- as.matrix(nr@houseKeepingVars[, -(1:(length(nds@inputHeaderValues) - length(nds@sampleReplicateGroups)))])
            class(Hkvartemp) <- "numeric"
            colmedianctr <- apply(Hkvartemp, 2, FUN="mean")
              
            for(i in 1:nrow(filterrawdata)) {
                nr@data2ctr[i,] <- unlist(sapply(1:ncol(filterrawdata), 
                                                 function(zd) { (filterrawdata[i, zd] / colmedianctr[zd]) * mean(colmedianctr) }))
            }
              
            nr@data2ctrlog <- log2(nr@data2ctr)
            colnames(nr@data2ctrlog) <- colnames(nds@filterrawdata)
              
            nr
        })

#' @rdname performVSNNormalization
setMethod("performVSNNormalization", "NormalyzerResults",
        function(nr) {
            nds <- nr@nds
            filterrawdata <- nds@filterrawdata
            nr@data2vsn <- vsn::justvsn(filterrawdata)
            nr
        })

#' @rdname performQuantileNormalization
setMethod("performQuantileNormalization", "NormalyzerResults",
        function(nr) {
            stopifnot(!is.null(nr@data2log2))
              
            nr@data2quantile <- preprocessCore::normalize.quantiles(nr@data2log2, copy=TRUE)
            nr
        })

#' @rdname performSMADNormalization
setMethod("performSMADNormalization", "NormalyzerResults",
        function(nr) {
            mediandata <- apply(nr@data2log2, 2, "median", na.rm=TRUE)
            maddata <- apply(nr@data2log2, 2, function(x) stats::mad(x, na.rm=TRUE))
            nr@data2mad <- t(apply(nr@data2log2, 1, function(x) ((x - mediandata) / maddata)))
            nr@data2mad <- nr@data2mad + mean(mediandata)
            nr
        })

#' @rdname performCyclicLoessNormalization
setMethod("performCyclicLoessNormalization", "NormalyzerResults",
        function(nr) {
            nr@data2loess <- limma::normalizeCyclicLoess(nr@data2log2, method="fast")
            nr
        })

#' @rdname performGlobalRLRNormalization
setMethod("performGlobalRLRNormalization", "NormalyzerResults",
        function(nr) {
            mediandata <- apply(nr@data2log2, 1, "median", na.rm=TRUE)
            isFirstSample <- TRUE
              
            for (j in 1:ncol(nr@data2log2)) {
                  
                # print(sprintf("Index: %i", j))
                  
                lrFit <- MASS::rlm(as.matrix(nr@data2log2[, j])~mediandata, na.action=stats::na.exclude)
                coeffs <- lrFit$coefficients
                coeffs2 <- coeffs[2]
                coeffs1 <- coeffs[1]
                  
                if (isFirstSample) {
                    value_to_assign <- (nr@data2log2[, j] - coeffs1) / coeffs2
                    globalfittedRLR <- (nr@data2log2[, j] - coeffs1) / coeffs2
                    isFirstSample <- FALSE
                }
                else {
                    globalfittedRLR <- cbind(globalfittedRLR, (nr@data2log2[, j] - coeffs1) / coeffs2)
                }
            }
              
            nr@globalfittedRLR <- globalfittedRLR
            colnames(nr@globalfittedRLR) <- colnames(nr@data2log2)
              
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
              
            nr@fittedLR <- matrix(, nrow=nrow(filterrawdata), ncol=0)
            nr@data2limloess <- matrix(, nrow=nrow(filterrawdata), ncol=0)
            nr@data2vsnrep <- matrix(, nrow=nrow(filterrawdata), ncol=0)
              
            for (sampleIndex in 1:length(firstIndices)) {
                  
                startColIndex <- firstIndices[sampleIndex]
                endColIndex <- lastIndices[sampleIndex]
                  
                # Median based LR normalization
                mediandata <- apply(nr@data2log2[, startColIndex:endColIndex], 1, "median", na.rm=TRUE)
                  
                for (currentCol in startColIndex:endColIndex) {
                      
                    LRfit <- MASS::rlm(as.matrix(nr@data2log2[, currentCol])~mediandata, na.action=stats::na.exclude)
                    LRcoeffs <- LRfit$coefficients
                    nr@fittedLR <- cbind(nr@fittedLR, (nr@data2log2[, currentCol] - LRcoeffs[1]) / LRcoeffs[2])
                }
                  
                nr@data2limloess <- cbind(nr@data2limloess, limma::normalizeCyclicLoess(nr@data2log2[, startColIndex:endColIndex], method="fast"))
                nr@data2vsnrep <- cbind(nr@data2vsnrep, vsn::justvsn(as.matrix(filterrawdata[, startColIndex:endColIndex])))
            }
              
            colnames(nr@fittedLR) <- colnames(nr@data2log2)
            colnames(nr@data2limloess) <- colnames(nr@data2log2)
            colnames(nr@data2vsnrep) <- colnames(nr@data2log2)
              
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

