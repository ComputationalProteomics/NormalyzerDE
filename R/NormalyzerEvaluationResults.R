#' S4 class to represent normalization evaluations
#' 
#' @slot avgcvmem Average coefficient of variance per method
#' @slot avgcvmempdiff Percentage difference of mean coefficient of variance
#'  compared to log2-transformed data
#' @slot featureCVPerMethod CV calculated per feature and normalization method.  
#' @slot avgmadmem Average median absolute deviation 
#' @slot avgmadmempdiff Percentage difference of median absolute deviation
#'  compared to log2-transformed data
#' @slot avgvarmem Average variance per method
#' @slot avgvarmempdiff Percentage difference of mean variance compared
#'  to log2-transformed data
#' @slot nonsiganfdrlist List of 5% least variable entries based on ANOVA
#'  for log2-transformed data
#' @slot nonsiganfdrlistcvpdiff Coefficient of variance for least variable
#'  entries
#' @slot anfdr ANOVA FDR values
#' @slot anovaFDRWithNA ANOVA FDR with NA when not applicable
#' @slot anovaP ANOVA calculated p-values
#' @slot avgpercorsum Within group Pearson correlations
#' @slot avgspecorsum Within group Spearman correlations
#' @export
NormalyzerEvaluationResults <- setClass("NormalyzerEvaluationResults",
                                           representation(
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
                                               anovaFDRWithNA = "matrix",
                                               anovaP = "matrix",
                                               
                                               avgpercorsum = "list",
                                               avgspecorsum = "list"
                                           ),
                                           prototype=prototype())


#' Calculate coefficient of variation (relative standard deviation)
#' Stores percentage variation compared to 
#'
#' @param nr NormalyzerDE results object.
#' @param ner NormalyzerDE evaluation object.
#' @return None
#' @rdname calculateCV
#' @keywords internal
setGeneric(name="calculateCV", 
           function(ner, nr) standardGeneric("calculateCV"))

#' @rdname calculateCV
setMethod("calculateCV", "NormalyzerEvaluationResults",
          function(ner, nr) {
              
              sampleReplicateGroups <- nr@nds@sampleReplicateGroups
              methodCount <- length(getUsedMethodNames(nr))
              methodList <- getNormalizationMatrices(nr)
              numberFeatures <- nrow(methodList[[1]])
              rowCount <- length(levels(as.factor(unlist(sampleReplicateGroups))))
              
              matrixCvs <- matrix(nrow=numberFeatures, ncol=methodCount)
              avgCVPerNormAndReplicates <- matrix(nrow=rowCount, ncol=methodCount, byrow=TRUE)
              
              for (methodIndex in seq_len(methodCount)) {
                  
                  processedDataMatrix <- methodList[[methodIndex]]
                  tempCVMatrix <- matrix(nrow=nrow(processedDataMatrix), ncol=length(levels(as.factor(unlist(sampleReplicateGroups)))), byrow=TRUE)
                  
                  for (i in seq_len(nrow(processedDataMatrix))) {
                      
                      tempCV <- RcmdrMisc::numSummary(processedDataMatrix[i, ], statistics=c("cv"), groups=unlist(sampleReplicateGroups))
                      tempCVMatrix[i, ] <- tempCV$table
                  }
                  
                  tempCVMatrixSum <- apply(tempCVMatrix, 2, mean, na.rm=TRUE)
                  avgCVPerNormAndReplicates[, methodIndex] <- tempCVMatrixSum * 100
                  
                  
                  # Calculate sample-wise CV
                  cv <- function(row) {
                      stDev <- stats::sd(row, na.rm=TRUE)
                      meanVal <- mean(unlist(row), na.rm=TRUE)
                      stDev / mean(meanVal) * 100
                  }
                  methodCvs <- apply(processedDataMatrix, 1, cv)
                  matrixCvs[, methodIndex] <- methodCvs
              }
              
              # Requires log normalized data to be at first index
              cvPercentVarFromLog <- vapply(
                  seq_len(ncol(avgCVPerNormAndReplicates)), 
                      function (sampleIndex) (mean(avgCVPerNormAndReplicates[, sampleIndex]) * 100) / mean(avgCVPerNormAndReplicates[, 1]),
                  0
              )
              
              ner@avgcvmem <- avgCVPerNormAndReplicates
              ner@avgcvmempdiff <- cvPercentVarFromLog
              ner@featureCVPerMethod <- matrixCvs
              
              ner
          }
)

#' Calculate median absolute deviation measures
#'
#' @param nr Normalyzer results object.
#' @param ner NormalyzerDE evaluation object.
#' @return None
#' @rdname calculateMAD
#' @keywords internal
setGeneric(name="calculateMAD", 
           function(ner, nr) standardGeneric("calculateMAD"))

#' @rdname calculateMAD
setMethod("calculateMAD", "NormalyzerEvaluationResults",
          function(ner, nr) {
              
              # Setup
              sampleReplicateGroups <- nr@nds@sampleReplicateGroups
              methodCount <- length(getUsedMethodNames(nr))
              methodList <- getNormalizationMatrices(nr)
              rowCount <- length(levels(as.factor(unlist(sampleReplicateGroups))))
              
              # Start
              avgmadmem <- matrix(nrow=rowCount, ncol=methodCount, byrow=TRUE)
              indexList <- getIndexList(sampleReplicateGroups)
              
              for (methodIndex in seq_len(methodCount)) {
                  
                  processedDataMatrix <- methodList[[methodIndex]]

                  medianAbsDevMem <- matrix(nrow=nrow(processedDataMatrix), 
                                            ncol=length(levels(as.factor(unlist(sampleReplicateGroups)))), 
                                            byrow=TRUE)
                  
                  for (sampleIndex in seq_len(length(names(indexList)))) {
                      repVal <- names(indexList)[sampleIndex]
                      cols <- indexList[[repVal]]
                      medianAbsDevMem[, sampleIndex] <- apply(processedDataMatrix[, cols], 1, function(x) { stats::mad(x, na.rm=TRUE) })
                      rowNonNACount <- apply(processedDataMatrix[, cols], 1, function(x) { sum(!is.na(x)) }) - 1
                  }
                  
                  temmadmatsum <- apply(medianAbsDevMem, 2, mean, na.rm=TRUE)
                  avgmadmem[, methodIndex] <- temmadmatsum
              }
              
              avgmadmempdiff <- vapply(seq_len(ncol(avgmadmem)), 
                                       function (sampleIndex) (mean(avgmadmem[, sampleIndex]) * 100) / mean(avgmadmem[, 1]),
                                       0)
              ner@avgmadmem <- avgmadmem
              ner@avgmadmempdiff <- avgmadmempdiff
              
              ner
          }
)

#' Calculate average variation measures
#'
#' @param ner NormalyzerDE evaluation object.
#' @param nr NormalyzerDE results object.
#' @return None
#' @rdname calculateAvgVar
#' @keywords internal
setGeneric(name="calculateAvgVar", 
           function(ner, nr) standardGeneric("calculateAvgVar"))

#' @rdname calculateAvgVar
setMethod("calculateAvgVar", "NormalyzerEvaluationResults",
          function(ner, nr) {
              
              # Setup
              sampleReplicateGroups <- nr@nds@sampleReplicateGroups
              methodCount <- length(getUsedMethodNames(nr))
              methodList <- getNormalizationMatrices(nr)
              rowCount <- length(levels(as.factor(unlist(sampleReplicateGroups))))
              
              # Start
              avgvarmem <- matrix(nrow=rowCount, ncol=methodCount, byrow=TRUE)
              indexList <- getIndexList(sampleReplicateGroups)
              
              for (methodIndex in seq_len(methodCount)) {
                  
                  replicateGroupVariance <- vector()
                  rowVariances <- vector()
                  
                  processedDataMatrix <- methodList[[methodIndex]]
                  
                  for (sampleIndex in seq_len(length(names(indexList)))) {
                      
                      repVal <- names(indexList)[sampleIndex]
                      cols <- indexList[[repVal]]
                      
                      rowNonNACount <- apply(processedDataMatrix[, cols], 1, function(x) { sum(!is.na(x)) }) - 1
                      rowVariances <- rowNonNACount * apply(processedDataMatrix[, cols], 1, function(x) { stats::var(x, na.rm=TRUE) })
                      replicateGroupVariance <- c(replicateGroupVariance, sum(rowVariances, na.rm=TRUE) / sum(rowNonNACount, na.rm=TRUE))
                  }
                  
                  avgvarmem[, methodIndex] <- replicateGroupVariance
                  
              }
              
              avgvarmempdiff <- vapply(
                  seq_len(ncol(avgvarmem)), 
                  function (sampleIndex) (mean(avgvarmem[, sampleIndex]) * 100) / mean(avgvarmem[, 1]),
                  0
              )
              
              ner@avgvarmem <- avgvarmem
              ner@avgvarmempdiff <- avgvarmempdiff
              
              ner
          }
)


#' Calculate significance measures (p-value and FDR-value) for ANOVA and
#' Kruskal-Wallis
#'
#' @param ner NormalyzerDE evaluation object.
#' @param nr Normalyzer results object.
#' @param categoricalAnova Categorical groupwise comparison instead of numeric
#' @return None
#' @rdname calculateSignificanceMeasures
#' @keywords internal
setGeneric(name="calculateSignificanceMeasures", 
           function(ner, nr, categoricalAnova) standardGeneric("calculateSignificanceMeasures"))

#' @rdname calculateSignificanceMeasures
setMethod("calculateSignificanceMeasures", "NormalyzerEvaluationResults",
          function(ner, nr, categoricalAnova=FALSE) {
              
              sampleReplicateGroups <- nr@nds@sampleReplicateGroups
              methodCount <- length(getUsedMethodNames(nr))
              methodList <- getNormalizationMatrices(nr)

              anovaPValsWithNA <- matrix(NA, ncol=methodCount, nrow=nrow(methodList[[1]]))
              anovaFDRs <- list()
              anovaFDRsWithNA <- matrix(NA, ncol=methodCount, nrow=nrow(methodList[[1]]))

              for (methodIndex in seq_len(methodCount)) {
                  
                  processedDataMatrix <- methodList[[methodIndex]]
                  naFilterContrast <- getRowNAFilterContrast(
                      processedDataMatrix,
                      sampleReplicateGroups,
                      2)
                  
                  dataStoreReplicateNAFiltered <- processedDataMatrix[naFilterContrast,]

                  if (categoricalAnova) {
                      testLevels <- factor(sampleReplicateGroups)
                  }
                  else {
                      testLevels <- sampleReplicateGroups
                  }
                  
                  anovaPValCol <- apply(
                      dataStoreReplicateNAFiltered, 
                      1, 
                      function(sampleIndex) 
                          summary(stats::aov(unlist(sampleIndex)~testLevels))[[1]][[5]][1])
                  
                  anovaPValsWithNA[naFilterContrast, methodIndex] <- anovaPValCol
                  anovaFDRCol <- stats::p.adjust(anovaPValCol, method="BH")
                  anovaFDRs[[methodIndex]] <- anovaFDRCol
                  anovaFDRsWithNA[naFilterContrast, methodIndex] <- anovaFDRCol
              }
              
              # Finds to 5% of least DE variables in log2 data based on ANOVA
              minVal <- min(utils::head(rev(sort(anovaFDRs[[1]])), n=(5 * length(anovaFDRs[[1]]) / 100)))
              nbrAboveThres <- sum(anovaFDRs[[1]] >= minVal)
              if (!is.infinite(minVal) && nbrAboveThres > 0) {
                  lowlyVariableFeatures <- which(anovaFDRs[[1]] >= min(utils::head(rev(sort(anovaFDRs[[1]])), n=(5 * length(anovaFDRs[[1]]) / 100))))
                  
                  nonsiganfdrlistcv <- vector()
                  for (mlist in seq_len(methodCount)) {
                      tmpdata <- methodList[[mlist]][lowlyVariableFeatures, ]
                      nonsiganfdrlistcv[mlist] <- mean(apply(tmpdata, 1, function(sampleIndex) raster::cv(sampleIndex, na.rm=TRUE)), na.rm=TRUE)
                  }
                  
                  nonsiganfdrlistcvpdiff <- vapply(
                      seq_len(length(nonsiganfdrlistcv)), 
                      function(sampleIndex) (nonsiganfdrlistcv[sampleIndex] * 100) / nonsiganfdrlistcv[1],
                      0
                  )
                  ner@nonsiganfdrlist <- nonsiganfdrlistcv
                  ner@nonsiganfdrlistcvpdiff <- nonsiganfdrlistcvpdiff
              }
              else {
                  warning("Too few successful ANOVA calculations to generate lowly variable features")
              }
              
              ner@anovaP <- anovaPValsWithNA
              
              ner@anfdr <- anovaFDRs
              ner@anovaFDRWithNA = anovaFDRsWithNA

              ner
          }
)


