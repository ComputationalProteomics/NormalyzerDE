#' S4 class to represent normalization evaluations
#' 
#' @slot avgcvmem TODO: Look into
#' @slot avgmadmem TODO: Look into
#' @slot avgvarmem TODO: Look into
#' @slot avgcvmempdiff TODO: Look into
#' @slot avgmadmempdiff TODO: Look into
#' @slot avgvarmempdiff TODO: Look into
#' @slot nonsiganfdrlist TODO: Look into
#' @slot nonsiganfdrlistcvpdiff TODO: Look into
#' @slot anfdr TODO: Look into
#' @slot kwfdr TODO: Look into
#' @slot avgpercorsum TODO: Look into
#' @slot avgspecorsum TODO: Look into
#' @export

NormalizationEvaluationResults <- setClass("NormalizationEvaluationResults",
                                           slots=c(
                                               avgcvmem = "matrix",
                                               avgcvmempdiff = "numeric",
                                               
                                               avgmadmem = "matrix",
                                               avgmadmempdiff = "numeric",
                                               
                                               avgvarmem = "matrix",
                                               avgvarmempdiff = "numeric",
                                               
                                               nonsiganfdrlist = "numeric",
                                               nonsiganfdrlistcvpdiff = "numeric",
                                               
                                               anfdr = "matrix",
                                               kwfdr = "matrix",
                                               
                                               avgpercorsum = "list",
                                               avgspecorsum = "list"
                                           ),
                                           prototype=prototype())



#' Calculate coefficient of variation (relative standard deviation)
#' Stores percentage variation compared to 
#'
#' @param nr Normalyzer results object.
#' @return None
#' @rdname Calculate calculateCV
setGeneric(name="calculateCV", 
           function(ner, nr) standardGeneric("calculateCV"))

#' @rdname calculateCV
setMethod("calculateCV", "NormalizationEvaluationResults",
          function(ner, nr) {
              
              sampleReplicateGroups <- nr@nds@sampleReplicateGroups
              methodCount <- length(getUsedMethodNames(nr))
              methodList <- getNormalizationMatrices(nr)
              
              rowCount <- length(levels(as.factor(unlist(sampleReplicateGroups))))
              
              avgCVPerNormAndReplicates <- matrix(nrow=rowCount, ncol=methodCount, byrow=TRUE)
              for (methodIndex in 1:methodCount) {
                  
                  processedDataMatrix <- methodList[[methodIndex]]
                  tempCVMatrix <- matrix(nrow=nrow(processedDataMatrix), ncol=length(levels(as.factor(unlist(sampleReplicateGroups)))), byrow=TRUE)
                  
                  for (i in 1:nrow(processedDataMatrix)) {
                      tempCV <- RcmdrMisc::numSummary(processedDataMatrix[i, ], statistics=c("cv"), groups=unlist(sampleReplicateGroups))
                      tempCVMatrix[i, ] <- tempCV$table
                  }
                  
                  # print(head(tempCVMatrix))
                  
                  tempCVMatrixSum <- apply(tempCVMatrix, 2, mean, na.rm=TRUE)
                  
                  # print("tempCVMatrixSum")
                  # print(head(tempCVMatrixSum))
                  
                  avgCVPerNormAndReplicates[, methodIndex] <- tempCVMatrixSum * 100
                  
                  # print(head(avgCVPerNormAndReplicates))
              }
              
              # Requires log normalized data to be at first index
              cvPercentVarFromLog <- sapply(1:ncol(avgCVPerNormAndReplicates), 
                                            function (sampleIndex) (mean(avgCVPerNormAndReplicates[, sampleIndex]) * 100) / mean(avgCVPerNormAndReplicates[, 1]))
              
              ner@avgcvmem <- avgCVPerNormAndReplicates
              ner@avgcvmempdiff <- cvPercentVarFromLog
              
              ner
          }
)

#' 
#'
#' @param nr Normalyzer results object.
#' @return None
#' @rdname calculateMAD
setGeneric(name="calculateMAD", 
           function(ner, nr) standardGeneric("calculateMAD"))

#' @rdname calculateMAD
setMethod("calculateMAD", "NormalizationEvaluationResults",
          function(ner, nr) {
              
              # Setup
              sampleReplicateGroups <- nr@nds@sampleReplicateGroups
              methodCount <- length(getUsedMethodNames(nr))
              methodList <- getNormalizationMatrices(nr)
              rowCount <- length(levels(as.factor(unlist(sampleReplicateGroups))))
              
              # Start
              avgmadmem <- matrix(nrow=rowCount, ncol=methodCount, byrow=TRUE)
              firstIndices <- getFirstIndicesInVector(sampleReplicateGroups, reverse=FALSE)
              lastIndices <- getFirstIndicesInVector(sampleReplicateGroups, reverse=TRUE)
              
              for (methodIndex in 1:methodCount) {
                  
                  processedDataMatrix <- methodList[[methodIndex]]
                  
                  medianAbsDevMem <- matrix(nrow=nrow(processedDataMatrix), 
                                            ncol=length(levels(as.factor(unlist(sampleReplicateGroups)))), 
                                            byrow=TRUE)
                  
                  for (sampleIndex in 1:length(firstIndices)) {
                      startColIndex <- firstIndices[sampleIndex]
                      endColIndex <- lastIndices[sampleIndex]
                      
                      medianAbsDevMem[, sampleIndex] <- apply(processedDataMatrix[, startColIndex:endColIndex], 1, function(x) { stats::mad(x, na.rm=TRUE) })
                      rowNonNACount <- apply(processedDataMatrix[, startColIndex:endColIndex], 1, function(x) { sum(!is.na(x)) }) - 1
                  }
                  
                  temmadmatsum <- apply(medianAbsDevMem, 2, mean, na.rm=TRUE)
                  avgmadmem[, methodIndex] <- temmadmatsum
              }
              
              avgmadmempdiff <- sapply(1:ncol(avgmadmem), function (sampleIndex) (mean(avgmadmem[, sampleIndex]) * 100) / mean(avgmadmem[, 1]))
              
              print("Assignment")
              
              ner@avgmadmem <- avgmadmem
              ner@avgmadmempdiff <- avgmadmempdiff
              
              print(head(ner@avgmadmem))
              print(head(ner@avgmadmempdiff))
              
              print("Return")
              # stop("")
              
              ner
          }
)

#' 
#'
#' @param nr Normalyzer results object.
#' @return None
#' @rdname calculateAvgVar
setGeneric(name="calculateAvgVar", 
           function(ner, nr) standardGeneric("calculateAvgVar"))

#' @rdname calculateAvgVar
setMethod("calculateAvgVar", "NormalizationEvaluationResults",
          function(ner, nr) {
              
              # Setup
              sampleReplicateGroups <- nr@nds@sampleReplicateGroups
              methodCount <- length(getUsedMethodNames(nr))
              methodList <- getNormalizationMatrices(nr)
              rowCount <- length(levels(as.factor(unlist(sampleReplicateGroups))))
              
              # Start
              avgvarmem <- matrix(nrow=rowCount, ncol=methodCount, byrow=TRUE)
              firstIndices <- getFirstIndicesInVector(sampleReplicateGroups, reverse=FALSE)
              lastIndices <- getFirstIndicesInVector(sampleReplicateGroups, reverse=TRUE)
              
              for (methodIndex in 1:methodCount) {
                  
                  replicateGroupVariance <- vector()
                  rowVariances <- vector()
                  
                  processedDataMatrix <- methodList[[methodIndex]]
                  
                  for (sampleIndex in 1:length(firstIndices)) {
                      
                      startColIndex <- firstIndices[sampleIndex]
                      endColIndex <- lastIndices[sampleIndex]
                      
                      rowNonNACount <- apply(processedDataMatrix[, startColIndex:endColIndex], 1, function(x) { sum(!is.na(x)) }) - 1
                      rowVariances <- rowNonNACount * apply(processedDataMatrix[, startColIndex:endColIndex], 1, function(x) { stats::var(x, na.rm=TRUE) })
                      
                      replicateGroupVariance <- c(replicateGroupVariance, sum(rowVariances, na.rm=TRUE) / sum(rowNonNACount, na.rm=TRUE))
                  }
                  
                  avgvarmem[, methodIndex] <- replicateGroupVariance
                  
              }
              
              avgvarmempdiff <- sapply(1:ncol(avgvarmem), function (sampleIndex) (mean(avgvarmem[, sampleIndex]) * 100) / mean(avgvarmem[, 1]))
              
              ner@avgvarmem <- avgvarmem
              ner@avgvarmempdiff <- avgvarmempdiff
              
              ner
          }
)


#' 
#'
#' @param nr Normalyzer results object.
#' @return None
#' @rdname calculateSignificanceMeasures
setGeneric(name="calculateSignificanceMeasures", 
           function(ner, nr) standardGeneric("calculateSignificanceMeasures"))

#' @rdname calculateSignificanceMeasures
setMethod("calculateSignificanceMeasures", "NormalizationEvaluationResults",
          function(ner, nr) {
              
              # Setup
              sampleReplicateGroups <- nr@nds@sampleReplicateGroups
              methodCount <- length(getUsedMethodNames(nr))
              methodList <- getNormalizationMatrices(nr)
              rowCount <- length(levels(as.factor(unlist(sampleReplicateGroups))))
              
              anovaPVal <- vector()
              krusWalPVal <- vector()
              anovaFDR <- vector()
              krusValFDR <- vector()
              
              firstIndices <- getFirstIndicesInVector(sampleReplicateGroups, reverse=FALSE)
              lastIndices <- getFirstIndicesInVector(sampleReplicateGroups, reverse=TRUE)

              for (methodIndex in 1:methodCount) {
                  
                  processedDataMatrix <- methodList[[methodIndex]]
                  rowNonNACount <- vector()
                  
                  for (sampleIndex in 1:length(firstIndices)) {
                      
                      startColIndex <- firstIndices[sampleIndex]
                      endColIndex <- lastIndices[sampleIndex]
                      
                      rowNonNACount <- apply(processedDataMatrix[, startColIndex:endColIndex], 1, function(x) { sum(!is.na(x)) }) - 1
                  }
                  
                  # ANOVA
                  nbsNAperLine <- rowSums(is.na(processedDataMatrix))
                  
                  # Retrieving lines where at least half is non-NAs
                  datastoretmp <- processedDataMatrix[nbsNAperLine < ncol(processedDataMatrix) / 2, ]
                  dataStoreReplicateNAFiltered <- filterLinesWithEmptySamples(datastoretmp, sampleReplicateGroups)
                  
                  anovaPVal <- cbind(anovaPVal, apply(dataStoreReplicateNAFiltered, 1, function(sampleIndex) summary(stats::aov(unlist(sampleIndex)~sampleReplicateGroups))[[1]][[5]][1]))
                  anovaFDR <- cbind(anovaFDR, stats::p.adjust(anovaPVal[, methodIndex], method="BH"))

                  krusWalPVal <- cbind(krusWalPVal, apply(dataStoreReplicateNAFiltered, 1, function(sampleIndex) stats::kruskal.test(unlist(sampleIndex)~sampleReplicateGroups, na.action="na.exclude")[[3]][1]))
                  krusValFDR <- cbind(krusValFDR, stats::p.adjust(krusWalPVal[, methodIndex], method="BH"))
              }
              
              # Finds to 5% of least DE variables in log2 data based on ANOVA
              # Generates error if it doesn't find least DE peptides
              if (sum(anovaFDR[, 1] >= min(utils::head(rev(sort(anovaFDR[, 1])), n=(5 * nrow(anovaFDR) / 100)))) > 0) {
                  nonsiganfdrlist <- which(anovaFDR[, 1] >= min(utils::head(rev(sort(anovaFDR[, 1])), n=(5 * nrow(anovaFDR) / 100))))
              }
              
              nonsiganfdrlistcv <- vector()
              for (mlist in 1:methodCount) {
                  tmpdata <- methodList[[mlist]][nonsiganfdrlist, ]
                  nonsiganfdrlistcv[mlist] <- mean(apply(tmpdata, 1, function(sampleIndex) raster::cv(sampleIndex, na.rm=TRUE)), na.rm=TRUE)
              }
              
              nonsiganfdrlistcvpdiff <- sapply(1:length(nonsiganfdrlistcv), function(sampleIndex) (nonsiganfdrlistcv[sampleIndex] * 100) / nonsiganfdrlistcv[1])
              
              ner@nonsiganfdrlist <- nonsiganfdrlistcv
              ner@nonsiganfdrlistcvpdiff <- nonsiganfdrlistcvpdiff
              ner@anfdr <- anovaFDR
              ner@kwfdr <- krusValFDR
              
              ner
          }
)

#' 
#'
#' @param nr Normalyzer results object.
#' @return None
#' @rdname calculateAnfdr
setGeneric(name="calculateAnfdr", 
           function(nr) standardGeneric("calculateAnfdr"))

#' @rdname calculateAnfdr
setMethod("calculateAnfdr", "NormalizationEvaluationResults",
          function(nr) {
              nr
          }
)

#' 
#'
#' @param nr Normalyzer results object.
#' @return None
#' @rdname calculateKwfdr
setGeneric(name="calculateKwfdr", 
           function(nr) standardGeneric("calculateKwfdr"))

#' @rdname calculateKwfdr
setMethod("calculateKwfdr", "NormalizationEvaluationResults",
          function(nr) {
              nr
          }
)


#' 
#'
#' @param nr Normalyzer results object.
#' @return None
#' @rdname calculateAvgpercorsum
setGeneric(name="calculateAvgpercorsum", 
           function(nr) standardGeneric("calculateAvgpercorsum"))

#' @rdname calculateAvgpercorsum
setMethod("calculateAvgpercorsum", "NormalizationEvaluationResults",
          function(nr) {
              nr
          }
)

#' 
#'
#' @param nr Normalyzer results object.
#' @return None
#' @rdname calculateAvgspecorsum
setGeneric(name="calculateAvgspecorsum", 
           function(nr) standardGeneric("calculateAvgspecorsum"))

#' @rdname calculateAvgspecorsum
setMethod("calculateAvgspecorsum", "NormalizationEvaluationResults",
          function(nr) {
              nr
          }
)
























