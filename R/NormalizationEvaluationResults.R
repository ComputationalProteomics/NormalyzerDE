#' S4 class to represent normalization evaluations
#' 
#' @slot avgcvmem Average coefficient of variance per method
#' @slot avgmadmem Average median absolute deviation 
#' @slot avgvarmem Average variance per method
#' @slot avgcvmempdiff Percentage difference of mean coefficient of variance
#'  compared to log2-transformed data
#' @slot avgmadmempdiff Percentage difference of median absolute deviation
#'  compared to log2-transformed data
#' @slot avgvarmempdiff Percentage difference of mean variance compared
#'  to log2-transformed data
#' @slot nonsiganfdrlist List of 5% least variable entries based on ANOVA
#'  for log2-transformed data
#' @slot nonsiganfdrlistcvpdiff Coefficient of variance for least variable
#'  entries
#' @slot anfdr TODO: Look into (FDR for ANOVA)
#' @slot kwfdr TODO: Look into (FDR for Kruskal-Wallis)
#' @slot avgpercorsum TODO: Look into (Average Pearson correlation sum)
#' @slot avgspecorsum TODO: Look into (Average Spearman correlation sum)
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
                  
                  tempCVMatrixSum <- apply(tempCVMatrix, 2, mean, na.rm=TRUE)
                  avgCVPerNormAndReplicates[, methodIndex] <- tempCVMatrixSum * 100
              }
              
              # Requires log normalized data to be at first index
              cvPercentVarFromLog <- sapply(1:ncol(avgCVPerNormAndReplicates), 
                                            function (sampleIndex) (mean(avgCVPerNormAndReplicates[, sampleIndex]) * 100) / mean(avgCVPerNormAndReplicates[, 1]))
              
              ner@avgcvmem <- avgCVPerNormAndReplicates
              ner@avgcvmempdiff <- cvPercentVarFromLog
              
              ner
          }
)

#' Calculate median absolute deviation measures
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
              
              ner@avgmadmem <- avgmadmem
              ner@avgmadmempdiff <- avgmadmempdiff
              
              ner
          }
)

#' Calculate average variation measures
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


#' Calculate significance measures (p-value and FDR-value) for ANOVA and
#' Kruskal-Wallis
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



















