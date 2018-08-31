#' Representation of evaluation results by calculating performance measures
#' for an an NormalyzerResults instance
#' 
#' Contains the resulting information from the processing which subsequently
#' can be used to generate the quality assessment report.
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
#' @slot lowVarFeaturesCVs List of 5% least variable entries based on ANOVA
#'  for log2-transformed data
#' @slot lowVarFeaturesCVsPercDiff Coefficient of variance for least variable
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
                                               
                                               lowVarFeaturesCVs = "numeric",
                                               lowVarFeaturesCVsPercDiff = "numeric",
                                               
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
              
              nds <- nds(nr)
              sampleReplicateGroups <- sampleReplicateGroups(nds)
              sampleReplicateGroups <- sampleReplicateGroups(nds)
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
              
              nds <- nds(nr)
              sampleReplicateGroups <- sampleReplicateGroups(nds)
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
              
              nds <- nds(nr)
              sampleReplicateGroups <- sampleReplicateGroups(nds)
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
          function(ner, nr, categoricalAnova=TRUE) {
              
              nds <- nds(nr)
              sampleReplicateGroups <- sampleReplicateGroups(nds)
              methodList <- getNormalizationMatrices(nr)

              anovaPValsWithNA <- calculateANOVAPValues(methodList, sampleReplicateGroups, categoricalAnova)
              log2AnovaFDR <- stats::p.adjust(anovaPValsWithNA[, 1][!is.na(anovaPValsWithNA[, 1])], method="BH")
              
              lowVarFeaturesCVs <- findLowlyVariableFeatures(log2AnovaFDR, methodList)
              lowVarFeaturesCVsPercDiff <- vapply(
                      seq_len(length(lowVarFeaturesCVs)),
                      function(sampleIndex) {
                          (lowVarFeaturesCVs[sampleIndex] * 100) / lowVarFeaturesCVs[1]
                      },
                      0
                  )
              
              ner@anovaP <- anovaPValsWithNA
              ner@lowVarFeaturesCVs <- lowVarFeaturesCVs
              ner@lowVarFeaturesCVsPercDiff <- lowVarFeaturesCVsPercDiff
              
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

            nds <- nds(nr)
            methodlist <- getNormalizationMatrices(nr)
            allReplicateGroups <- sampleReplicateGroups(nds)
            sampleGroupsWithReplicates <- samplesGroupsWithReplicates(nds)
            
            ner@repCorPear <- calculateSummarizedCorrelationVector(
                methodlist, 
                allReplicateGroups, 
                sampleGroupsWithReplicates, "pearson"
            )
            ner@repCorSpear <- calculateSummarizedCorrelationVector(
                methodlist, 
                allReplicateGroups, 
                sampleGroupsWithReplicates, 
                "spearman"
            )
            
            ner
          }
)


