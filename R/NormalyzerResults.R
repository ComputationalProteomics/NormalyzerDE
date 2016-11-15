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

NormalyzerResults <- setClass("NormalyzerResults",
                              slots=c(
                                  nds = "NormalyzerDataset",
                                  
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
                              prototype=prototype(nds=NULL, normfinderMaxThreshold=1000, furtherNormalizationMinThreshold=50))

# setGeneric(name="", function(nr) standardGeneric(""))
# setMethod("", "NormalyzerResults",
#           function(nr) {
#               nds <- nr@nds
# 
#               nr
#           })

setGeneric(name="initializeResultsObject", function(nr) standardGeneric("initializeResultsObject"))
setGeneric(name="performNormalizations", function(nr) standardGeneric("performNormalizations"))
setGeneric(name="basicMetricNormalizations", function(nr) standardGeneric("basicMetricNormalizations"))
setGeneric(name="calculateHKdataForNormObj", function(nr) standardGeneric("calculateHKdataForNormObj"))

setGeneric(name="performVSNNormalization", function(nr) standardGeneric("performVSNNormalization"))
setGeneric(name="performQuantileNormalization", function(nr) standardGeneric("performQuantileNormalization"))
setGeneric(name="performSMADNormalization", function(nr) standardGeneric("performSMADNormalization"))
setGeneric(name="performCyclicLoessNormalization", function(nr) standardGeneric("performCyclicLoessNormalization"))
setGeneric(name="performGlobalRLRNormalization", function(nr) standardGeneric("performGlobalRLRNormalization"))
setGeneric(name="performReplicateBasedNormalizations", function(nr) standardGeneric("performReplicateBasedNormalizations"))

setGeneric(name="getMethodNames", function(nr) standardGeneric("getMethodNames"))
setGeneric(name="getSlotNameList", function(nr) standardGeneric("getSlotNameList"))



setMethod("initializeResultsObject", "NormalyzerResults",
          function(nr) {
              
              nds <- nr@nds

              nr@data2log2 <- log2(nds@filterrawdata)
              nr@data2GI <- matrix(nrow=nrow(nds@filterrawdata), ncol=ncol(nds@filterrawdata), byrow=T)
              nr@data2ctr <- matrix(nrow=nrow(nds@filterrawdata), ncol=ncol(nds@filterrawdata), byrow=T)
              nr@data2med <- matrix(nrow=nrow(nds@filterrawdata), ncol=ncol(nds@filterrawdata), byrow=T) 
              nr@data2mean <- matrix(nrow=nrow(nds@filterrawdata), ncol=ncol(nds@filterrawdata), byrow=T)
              
              nr
          })

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
          })

setMethod("basicMetricNormalizations", "NormalyzerResults",
              function(nr) {
                  ## Normalizing using average sample sum, sample mean and sample median
                  
                  nds <- nr@nds
                  
                  avgcolsum <- median(nds@colsum)
                  for (i in 1:nrow(nds@filterrawdata)) {
                      nr@data2GI[i, ] <- unlist(sapply(1:ncol(nds@filterrawdata), function(zd) { (nds@filterrawdata[i, zd] / nds@colsum[zd]) * (avgcolsum) }))
                      nr@data2med[i, ] <- unlist(sapply(1:ncol(nds@filterrawdata), function(zd) { (nds@filterrawdata[i, zd] / nds@medofdata[zd]) * mean(nds@medofdata) }))
                      nr@data2mean[i, ] <- unlist(sapply(1:ncol(nds@filterrawdata), function(zd) { (nds@filterrawdata[i, zd] / nds@meanofdata[zd]) * mean(nds@meanofdata) }))
                  } 
                  
                  nr@data2GI <- log2(nr@data2GI)
                  nr@data2med <- log2(nr@data2med)
                  nr@data2mean <- log2(nr@data2mean)
                  colnames(nr@data2GI) <- colnames(nr@data2log2)
                  colnames(nr@data2med) <- colnames(nr@data2log2)
                  colnames(nr@data2mean) <- colnames(nr@data2log2)
                  
                  nr
            })

setMethod("calculateHKdataForNormObj", "NormalyzerResults",
          function(nr) {
              nds <- nr@nds

              filterrawdata <- nds@filterrawdata
              
              Hkvartemp <- as.matrix(nr@houseKeepingVars[, -(1:(length(nds@inputHeaderValues) - length(nds@sampleReplicateGroups)))])
              class(Hkvartemp) <- "numeric"
              colmedianctr <- apply(Hkvartemp, 2, FUN="mean")
              
              for(i in 1:nrow(filterrawdata)) {
                  nr@data2ctr[i,] <- unlist(sapply(1:ncol(filterrawdata), 
                                                   function(zd) { (filterrawdata[i, zd] / colmedianctr[zd]) * (mean(colmedianctr)) }))
              }
              
              nr@data2ctrlog <- log2(nr@data2ctr)
              colnames(nr@data2ctrlog) <- colnames(nds@filterrawdata)
                            
              nr
          })

setMethod("performVSNNormalization", "NormalyzerResults",
          function(nr) {
              nds <- nr@nds
              filterrawdata <- nds@filterrawdata
              nr@data2vsn <- justvsn(filterrawdata)
              nr
          })

setMethod("performQuantileNormalization", "NormalyzerResults",
          function(nr) {
              stopifnot(!is.null(nr@data2log2))
              
              nr@data2quantile <- normalize.quantiles(nr@data2log2, copy=T)
              nr
          })

setMethod("performSMADNormalization", "NormalyzerResults",
          function(nr) {
              mediandata <- apply(nr@data2log2, 2, "median", na.rm=T)
              maddata <- apply(nr@data2log2, 2, function(x) mad(x, na.rm=T))
              nr@data2mad <- t(apply(nr@data2log2, 1, function(x) ((x - mediandata) / maddata)))
              nr@data2mad <- nr@data2mad + mean(mediandata)
              nr
          })

setMethod("performCyclicLoessNormalization", "NormalyzerResults",
          function(nr) {
              nr@data2loess <- normalizeCyclicLoess(nr@data2log2, method="fast")
              nr
          })

setMethod("performGlobalRLRNormalization", "NormalyzerResults",
          function(nr) {
              mediandata <- apply(nr@data2log2, 1, "median", na.rm=T)
              isFirstSample <- TRUE
              
              for (j in 1:ncol(nr@data2log2)) {
                  
                  # print(sprintf("Index: %i", j))
                  
                  lrFit <- rlm(as.matrix(nr@data2log2[, j])~mediandata, na.action=na.exclude)
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

setMethod("performReplicateBasedNormalizations", "NormalyzerResults",
          function(nr) {
              ## NORMALIZATION within REPLICATES RLR, VSN and Loess
              
              nds <- nr@nds
              
              sampleReplicateGroups <- nds@sampleReplicateGroups
              filterrawdata <- nds@filterrawdata

              firstIndices <- getFirstIndicesInVector(sampleReplicateGroups, reverse=F)
              lastIndices <- getFirstIndicesInVector(sampleReplicateGroups, reverse=T)
              
              stopifnot(length(firstIndices) == length(lastIndices))
              
              nr@fittedLR <- matrix(, nrow=nrow(filterrawdata), ncol=0)
              nr@data2limloess <- matrix(, nrow=nrow(filterrawdata), ncol=0)
              nr@data2vsnrep <- matrix(, nrow=nrow(filterrawdata), ncol=0)
              
              for (sampleIndex in 1:length(firstIndices)) {
                  
                  startColIndex <- firstIndices[sampleIndex]
                  endColIndex <- lastIndices[sampleIndex]
                  
                  # Median based LR normalization
                  mediandata <- apply(nr@data2log2[, startColIndex:endColIndex], 1, "median", na.rm=T)
                  
                  for (currentCol in startColIndex:endColIndex) {
                      
                      LRfit <- rlm(as.matrix(nr@data2log2[, currentCol])~mediandata, na.action=na.exclude)
                      LRcoeffs <- LRfit$coefficients
                      nr@fittedLR <- cbind(nr@fittedLR, (nr@data2log2[, currentCol] - LRcoeffs[1]) / LRcoeffs[2])
                  }
                  
                  nr@data2limloess <- cbind(nr@data2limloess, normalizeCyclicLoess(nr@data2log2[, startColIndex:endColIndex], method="fast"))
                  nr@data2vsnrep <- cbind(nr@data2vsnrep, justvsn(as.matrix(filterrawdata[, startColIndex:endColIndex])))
              }
              
              colnames(nr@fittedLR) <- colnames(nr@data2log2)
              colnames(nr@data2limloess) <- colnames(nr@data2log2)
              colnames(nr@data2vsnrep) <- colnames(nr@data2log2)
              
              nr
          })

setMethod("getMethodNames", "NormalyzerResults",
          function(nr) {
              
              usedMethodNames <- c()
              
              for (i in 1:length(slotNames)) {
                  slotName <- slotNames[i]
                  fieldValue <- slot(nr, slotName)
                  
                  if (!all(is.na(fieldValue))) {
                      outputName <- outputNames[i]
                      usedMethodNames <- c(usedMethodNames, outputName)
                  }
              }

              usedMethodNames
          })

setMethod("getSlotNameList", "NormalyzerResults",
          function(nr) {
              methodDataList <- c()
              
              for (i in 1:length(slotNames)) {
                  slotName <- slotNames[i]
                  fieldValue <- slot(nr, slotName)
                  
                  if (!all(is.na(fieldValue))) {
                      methodDataList <- c(methodDataList, slotName)
                  }
              }
              
              methodDataList
          })



