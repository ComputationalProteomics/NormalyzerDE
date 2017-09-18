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
                                               featureCVPerMethod = "matrix",
                                               
                                               avgmadmem = "matrix",
                                               avgmadmempdiff = "numeric",
                                               
                                               avgvarmem = "matrix",
                                               avgvarmempdiff = "numeric",
                                               
                                               nonsiganfdrlist = "numeric",
                                               nonsiganfdrlistcvpdiff = "numeric",
                                               
                                               anfdr = "list",
                                               kwfdr = "list",
                                               anovaFDRWithNA = "matrix",
                                               krusWalFDRWithNA = "matrix",
                                               anova_p = "matrix",
                                               
                                               pairwise_comps = "list",
                                               pairwise_comps_fdr = "list",
                                               
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
              numberFeatures <- nrow(methodList[[1]])
              rowCount <- length(levels(as.factor(unlist(sampleReplicateGroups))))
              
              matrix_cvs <- matrix(nrow=numberFeatures, ncol=methodCount)
              avgCVPerNormAndReplicates <- matrix(nrow=rowCount, ncol=methodCount, byrow=TRUE)
              
              for (methodIndex in 1:methodCount) {
                  
                  processedDataMatrix <- methodList[[methodIndex]]
                  tempCVMatrix <- matrix(nrow=nrow(processedDataMatrix), ncol=length(levels(as.factor(unlist(sampleReplicateGroups)))), byrow=TRUE)
                  
                  for (i in 1:nrow(processedDataMatrix)) {
                      # tempCV <- rcmdrNumSummary(processedDataMatrix[i, ], statistics=c("cv"), groups=unlist(sampleReplicateGroups))
                      tempCV <- RcmdrMisc::numSummary(processedDataMatrix[i, ], statistics=c("cv"), groups=unlist(sampleReplicateGroups))
                      tempCVMatrix[i, ] <- tempCV$table
                  }
                  
                  tempCVMatrixSum <- apply(tempCVMatrix, 2, mean, na.rm=TRUE)
                  avgCVPerNormAndReplicates[, methodIndex] <- tempCVMatrixSum * 100
                  
                  
                  # Calculate sample-wise CV
                  cv <- function(row) {
                      st_dev <- sd(row, na.rm=TRUE)
                      mean_val <- mean(unlist(row), na.rm=TRUE)
                      st_dev / mean(mean_val) * 100
                  }
                  method_cvs <- apply(processedDataMatrix, 1, cv)
                  matrix_cvs[,methodIndex] <- method_cvs
                  
              }
              
              # Requires log normalized data to be at first index
              cvPercentVarFromLog <- sapply(1:ncol(avgCVPerNormAndReplicates), 
                                            function (sampleIndex) (mean(avgCVPerNormAndReplicates[, sampleIndex]) * 100) / mean(avgCVPerNormAndReplicates[, 1]))
              
              ner@avgcvmem <- avgCVPerNormAndReplicates
              ner@avgcvmempdiff <- cvPercentVarFromLog
              ner@featureCVPerMethod <- matrix_cvs
              
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
           function(ner, nr, categorical_anova, var_filter_frac) standardGeneric("calculateSignificanceMeasures"))

#' @rdname calculateSignificanceMeasures
setMethod("calculateSignificanceMeasures", "NormalizationEvaluationResults",
          function(ner, nr, 
                   categorical_anova=FALSE, 
                   var_filter_frac=NULL) {
              
              # browser()
              
              # Setup
              sampleReplicateGroups <- nr@nds@sampleReplicateGroups
              methodCount <- length(getUsedMethodNames(nr))
              methodList <- getNormalizationMatrices(nr)

              anovaPValsWithNA <- matrix(NA, ncol=methodCount, nrow=nrow(methodList[[1]]))
              anovaFDRs <- list()
              anovaFDRsWithNA <- matrix(NA, ncol=methodCount, nrow=nrow(methodList[[1]]))
              krusValFDRs <- list()
              krusWalFDRsWithNA <- matrix(NA, ncol=methodCount, nrow=nrow(methodList[[1]]))
              
              firstIndices <- getFirstIndicesInVector(sampleReplicateGroups, reverse=FALSE)
              lastIndices <- getFirstIndicesInVector(sampleReplicateGroups, reverse=TRUE)

              for (methodIndex in 1:methodCount) {
                  
                  # print(paste("CV col length: ", length(ner@featureCVPerMethod[,methodIndex])))
                  
                  processedDataMatrix <- methodList[[methodIndex]]
                  naFilterContrast <- getRowNAFilterContrast(processedDataMatrix,
                                                             sampleReplicateGroups,
                                                             var_filter_frac=var_filter_frac)
                  
                  # print(paste("NA filter contrast length: ", length(naFilterContrast)))
                  
                  dataStoreReplicateNAFiltered <- processedDataMatrix[naFilterContrast,]

                  if (categorical_anova) {
                      testLevels <- factor(sampleReplicateGroups)
                  }
                  else {
                      testLevels <- sampleReplicateGroups
                  }
                  
                  # browser()
                  
                  
                  
                  anovaPValCol <- apply(dataStoreReplicateNAFiltered, 1, 
                                        function(sampleIndex) summary(stats::aov(unlist(sampleIndex)~testLevels))[[1]][[5]][1])
                  krusWalPValCol <- apply(dataStoreReplicateNAFiltered, 1, 
                                          function(sampleIndex) stats::kruskal.test(unlist(sampleIndex)~testLevels, na.action="na.exclude")[[3]][1])
                  
                  # print(paste("anovaPVals length: ", length(anovaPValCol)))
                  
                  anovaPValsWithNA[naFilterContrast, methodIndex] <- anovaPValCol
                  anovaFDRCol <- stats::p.adjust(anovaPValCol, method="BH")
                  anovaFDRs[[methodIndex]] <- anovaFDRCol
                  anovaFDRsWithNA[naFilterContrast, methodIndex] <- anovaFDRCol
                  
                  krusWalFDRCol <- stats::p.adjust(krusWalPValCol, method="BH")
                  krusValFDRs[[methodIndex]] <- krusWalFDRCol
                  krusWalFDRsWithNA[naFilterContrast, methodIndex] <- krusWalFDRCol
              }
              
              # Finds to 5% of least DE variables in log2 data based on ANOVA
              # Generates error if it doesn't find least DE peptides
              if (sum(anovaFDRs[[1]] >= min(utils::head(rev(sort(anovaFDRs[[1]])), n=(5 * length(anovaFDRs[[1]]) / 100)))) > 0) {
                  lowlyVariableFeatures <- which(anovaFDRs[[1]] >= min(utils::head(rev(sort(anovaFDRs[[1]])), n=(5 * length(anovaFDRs[[1]]) / 100))))
              }
              
              # if (sum(anovaFDRs[, 1] >= min(utils::head(rev(sort(anovaFDRs[, 1])), n=(5 * nrow(anovaFDRs) / 100)))) > 0) {
              #     lowlyVariableFeatures <- which(anovaFDRs[, 1] >= min(utils::head(rev(sort(anovaFDRs[, 1])), n=(5 * nrow(anovaFDRs) / 100))))
              # }

              nonsiganfdrlistcv <- vector()
              for (mlist in 1:methodCount) {
                  tmpdata <- methodList[[mlist]][lowlyVariableFeatures, ]
                  nonsiganfdrlistcv[mlist] <- mean(apply(tmpdata, 1, function(sampleIndex) raster::cv(sampleIndex, na.rm=TRUE)), na.rm=TRUE)
              }
              
              nonsiganfdrlistcvpdiff <- sapply(1:length(nonsiganfdrlistcv), function(sampleIndex) (nonsiganfdrlistcv[sampleIndex] * 100) / nonsiganfdrlistcv[1])
              
              ner@nonsiganfdrlist <- nonsiganfdrlistcv
              ner@nonsiganfdrlistcvpdiff <- nonsiganfdrlistcvpdiff
              
              ner@anova_p <- anovaPValsWithNA
              
              ner@anfdr <- anovaFDRs
              ner@kwfdr <- krusValFDRs
              ner@anovaFDRWithNA = anovaFDRsWithNA
              ner@krusWalFDRWithNA = krusWalFDRsWithNA

              ner
          }
)


#' Calculate Welch t-test comparisons between target samples
#'
#' @param nr Normalyzer results object.
#' @return None
#' @rdname calculatePairwiseComparisons
setGeneric(name="calculatePairwiseComparisons", 
           function(ner, nr, comparisons) standardGeneric("calculatePairwiseComparisons"))

#' @rdname calculatePairwiseComparisons
setMethod("calculatePairwiseComparisons", "NormalizationEvaluationResults",
          function(ner, nr, comparisons) {
              
              # Setup
              sampleReplicateGroups <- nr@nds@sampleReplicateGroups
              methodCount <- length(getUsedMethodNames(nr))
              methodList <- getNormalizationMatrices(nr)

              firstIndices <- getFirstIndicesInVector(sampleReplicateGroups, reverse=FALSE)
              lastIndices <- getFirstIndicesInVector(sampleReplicateGroups, reverse=TRUE)
              
              pairwisePList <- list()
              pairwiseFDRList <- list()
              for (comp in comparisons) {
                  pairwisePList[[comp]] <- matrix(NA, ncol=methodCount, nrow=nrow(methodList[[1]]))
                  pairwiseFDRList[[comp]] <- matrix(NA, ncol=methodCount, nrow=nrow(methodList[[1]]))
              }

              do_t_test <- function(c1_vals, c2_vals, default=NA) {
                  if (length(na.omit(c1_vals)) > 1 && length(na.omit(c2_vals)) > 1) {
                      tryCatch(stats::t.test(c1_vals, c2_vals)[[3]],
                               error=function(x) NA)
                  }
                  else {
                      NA
                  }
              }
                            
              for (methodIndex in 1:methodCount) {
                  
                  processedDataMatrix <- methodList[[methodIndex]]
                  naFilterContrast <- getRowNAFilterContrast(processedDataMatrix, sampleReplicateGroups)
                  dataStoreReplicateNAFiltered <- processedDataMatrix[naFilterContrast,]
                  
                  for (comp in comparisons) {
                      
                      sample1 <- strtoi(substr(comp, 1, 1))
                      sample2 <- strtoi(substr(comp, 2, 2))
                      
                      s1start <- firstIndices[sample1]
                      s1end <- lastIndices[sample1]
                      s2start <- firstIndices[sample2]
                      s2end <- lastIndices[sample2]
                      
                      
                      welchPValCol <- apply(dataStoreReplicateNAFiltered, 1, function(row) do_t_test(row[s1start:s1end], row[s2start:s2end], default=NA))
                      welchFDRCol <- stats::p.adjust(welchPValCol, method="BH")

                      pairwisePList[[comp]][naFilterContrast, methodIndex] <- welchPValCol
                      pairwiseFDRList[[comp]][naFilterContrast, methodIndex] <- welchFDRCol
                  }
              }
              
              ner@pairwise_comps <- pairwisePList
              ner@pairwise_comps_fdr <- pairwiseFDRList
              ner
          }
)







