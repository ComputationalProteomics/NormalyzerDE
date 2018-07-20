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
#' @slot repCorPear Within group Pearson correlations
#' @slot repCorSpear Within group Spearman correlations
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
                                               
                                               repCorPear = "list",
                                               repCorSpear = "list"
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
              methodList <- getNormalizationMatrices(nr)

              avgCVPerNormAndReplicates <- calculateReplicateCV(methodList, sampleReplicateGroups)
              ner@avgcvmem <- avgCVPerNormAndReplicates
              ner@featureCVPerMethod <- calculateFeatureCV(methodList)
              ner@avgcvmempdiff <- calculatePercentageAvgDiffInMat(avgCVPerNormAndReplicates)
              
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
              
              sampleReplicateGroups <- nr@nds@sampleReplicateGroups
              methodList <- getNormalizationMatrices(nr)
              avgmadmem <- calculateAvgMadMem(methodList, sampleReplicateGroups)

              avgMadPdiff<- calculatePercentageAvgDiffInMat(avgmadmem)

              ner@avgmadmem <- avgmadmem
              ner@avgmadmempdiff <- avgMadPdiff
              
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
              
              sampleReplicateGroups <- nr@nds@sampleReplicateGroups
              methodList <- getNormalizationMatrices(nr)

              avgVarianceMat <- calculateAvgReplicateVariation(methodList, sampleReplicateGroups)
              ner@avgvarmem <- avgVarianceMat
              ner@avgvarmempdiff <- calculatePercentageAvgDiffInMat(avgVarianceMat)
              
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

#' Pearson and Spearman correlation calculations for methods and samples
#' Calculates internal correlation per condition
#' 
#' @param nr Normalyzer results object with calculated results.
#' @param ner Normalyzer evaluation object.
#' @return ner Normalyzer evaluation object with attached evaluation results.
#' @keywords internal
setGeneric(name="calculateCorrelations", 
           function(ner, nr) standardGeneric("calculateCorrelations"))

#' @rdname calculateCorrelations
setMethod("calculateCorrelations", "NormalyzerEvaluationResults",
          function(ner, nr) {

            methodlist <- getNormalizationMatrices(nr)
            allReplicateGroups <- nr@nds@sampleReplicateGroups
            sampleGroupsWithReplicates <- nr@nds@samplesGroupsWithReplicates
            
            ner@repCorPear <- calculateSummarizedCorrelationVector(
                methodlist, allReplicateGroups, sampleGroupsWithReplicates, "pearson")
            ner@repCorSpear <- calculateSummarizedCorrelationVector(
                methodlist, allReplicateGroups, sampleGroupsWithReplicates, "spearman")
            ner
          }
)


