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


setGeneric("avgcvmem", function(object) { standardGeneric("avgcvmem") })
setMethod("avgcvmem", signature(object="NormalyzerEvaluationResults"), 
          function(object) { slot(object, "avgcvmem") })
setGeneric("avgcvmem<-", function(object, value) { standardGeneric("avgcvmem<-") })
setReplaceMethod("avgcvmem", signature(object="NormalyzerEvaluationResults"), 
                 function(object, value) { 
                     slot(object, "avgcvmem") <- value
                     validObject(object)
                     object
                 })

setGeneric("avgcvmempdiff", function(object) { standardGeneric("avgcvmempdiff") })
setMethod("avgcvmempdiff", signature(object="NormalyzerEvaluationResults"), 
          function(object) { slot(object, "avgcvmempdiff") })
setGeneric("avgcvmempdiff<-", function(object, value) { standardGeneric("avgcvmempdiff<-") })
setReplaceMethod("avgcvmempdiff", signature(object="NormalyzerEvaluationResults"), 
                 function(object, value) { 
                     slot(object, "avgcvmempdiff") <- value
                     validObject(object)
                     object
                 })

setGeneric("featureCVPerMethod", function(object) { standardGeneric("featureCVPerMethod") })
setMethod("featureCVPerMethod", signature(object="NormalyzerEvaluationResults"), 
          function(object) { slot(object, "featureCVPerMethod") })
setGeneric("featureCVPerMethod<-", function(object, value) { standardGeneric("featureCVPerMethod<-") })
setReplaceMethod("featureCVPerMethod", signature(object="NormalyzerEvaluationResults"), 
                 function(object, value) { 
                     slot(object, "featureCVPerMethod") <- value
                     validObject(object)
                     object
                 })

setGeneric("avgmadmem", function(object) { standardGeneric("avgmadmem") })
setMethod("avgmadmem", signature(object="NormalyzerEvaluationResults"), 
          function(object) { slot(object, "avgmadmem") })
setGeneric("avgmadmem<-", function(object, value) { standardGeneric("avgmadmem<-") })
setReplaceMethod("avgmadmem", signature(object="NormalyzerEvaluationResults"), 
                 function(object, value) { 
                     slot(object, "avgmadmem") <- value
                     validObject(object)
                     object
                 })

setGeneric("avgmadmempdiff", function(object) { standardGeneric("avgmadmempdiff") })
setMethod("avgmadmempdiff", signature(object="NormalyzerEvaluationResults"), 
          function(object) { slot(object, "avgmadmempdiff") })
setGeneric("avgmadmempdiff<-", function(object, value) { standardGeneric("avgmadmempdiff<-") })
setReplaceMethod("avgmadmempdiff", signature(object="NormalyzerEvaluationResults"), 
                 function(object, value) { 
                     slot(object, "avgmadmempdiff") <- value
                     validObject(object)
                     object
                 })

setGeneric("avgvarmem", function(object) { standardGeneric("avgvarmem") })
setMethod("avgvarmem", signature(object="NormalyzerEvaluationResults"), 
          function(object) { slot(object, "avgvarmem") })
setGeneric("avgvarmem<-", function(object, value) { standardGeneric("avgvarmem<-") })
setReplaceMethod("avgvarmem", signature(object="NormalyzerEvaluationResults"), 
                 function(object, value) { 
                     slot(object, "avgvarmem") <- value
                     validObject(object)
                     object
                 })

setGeneric("avgvarmempdiff", function(object) { standardGeneric("avgvarmempdiff") })
setMethod("avgvarmempdiff", signature(object="NormalyzerEvaluationResults"), 
          function(object) { slot(object, "avgvarmempdiff") })
setGeneric("avgvarmempdiff<-", function(object, value) { standardGeneric("avgvarmempdiff<-") })
setReplaceMethod("avgvarmempdiff", signature(object="NormalyzerEvaluationResults"), 
                 function(object, value) { 
                     slot(object, "avgvarmempdiff") <- value
                     validObject(object)
                     object
                 })

setGeneric("lowVarFeaturesCVs", function(object) { standardGeneric("lowVarFeaturesCVs") })
setMethod("lowVarFeaturesCVs", signature(object="NormalyzerEvaluationResults"), 
          function(object) { slot(object, "lowVarFeaturesCVs") })
setGeneric("lowVarFeaturesCVs<-", function(object, value) { standardGeneric("lowVarFeaturesCVs<-") })
setReplaceMethod("lowVarFeaturesCVs", signature(object="NormalyzerEvaluationResults"), 
                 function(object, value) { 
                     slot(object, "lowVarFeaturesCVs") <- value
                     validObject(object)
                     object
                 })

setGeneric("lowVarFeaturesCVsPercDiff", function(object) { standardGeneric("lowVarFeaturesCVsPercDiff") })
setMethod("lowVarFeaturesCVsPercDiff", signature(object="NormalyzerEvaluationResults"), 
          function(object) { slot(object, "lowVarFeaturesCVsPercDiff") })
setGeneric("lowVarFeaturesCVsPercDiff<-", function(object, value) { standardGeneric("lowVarFeaturesCVsPercDiff<-") })
setReplaceMethod("lowVarFeaturesCVsPercDiff", signature(object="NormalyzerEvaluationResults"), 
                 function(object, value) { 
                     slot(object, "lowVarFeaturesCVsPercDiff") <- value
                     validObject(object)
                     object
                 })

setGeneric("anfdr", function(object) { standardGeneric("anfdr") })
setMethod("anfdr", signature(object="NormalyzerEvaluationResults"), 
          function(object) { slot(object, "anfdr") })
setGeneric("anfdr<-", function(object, value) { standardGeneric("anfdr<-") })
setReplaceMethod("anfdr", signature(object="NormalyzerEvaluationResults"), 
                 function(object, value) { 
                     slot(object, "anfdr") <- value
                     validObject(object)
                     object
                 })

setGeneric("anovaFDRWithNA", function(object) { standardGeneric("anovaFDRWithNA") })
setMethod("anovaFDRWithNA", signature(object="NormalyzerEvaluationResults"), 
          function(object) { slot(object, "anovaFDRWithNA") })
setGeneric("anovaFDRWithNA<-", function(object, value) { standardGeneric("anovaFDRWithNA<-") })
setReplaceMethod("anovaFDRWithNA", signature(object="NormalyzerEvaluationResults"), 
                 function(object, value) { 
                     slot(object, "anovaFDRWithNA") <- value
                     validObject(object)
                     object
                 })

setGeneric("anovaP", function(object) { standardGeneric("anovaP") })
setMethod("anovaP", signature(object="NormalyzerEvaluationResults"), 
          function(object) { slot(object, "anovaP") })
setGeneric("anovaP<-", function(object, value) { standardGeneric("anovaP<-") })
setReplaceMethod("anovaP", signature(object="NormalyzerEvaluationResults"), 
                 function(object, value) { 
                     slot(object, "anovaP") <- value
                     validObject(object)
                     object
                 })

setGeneric("repCorPear", function(object) { standardGeneric("repCorPear") })
setMethod("repCorPear", signature(object="NormalyzerEvaluationResults"), 
          function(object) { slot(object, "repCorPear") })
setGeneric("repCorPear<-", function(object, value) { standardGeneric("repCorPear<-") })
setReplaceMethod("repCorPear", signature(object="NormalyzerEvaluationResults"), 
                 function(object, value) { 
                     slot(object, "repCorPear") <- value
                     validObject(object)
                     object
                 })

setGeneric("repCorSpear", function(object) { standardGeneric("repCorSpear") })
setMethod("repCorSpear", signature(object="NormalyzerEvaluationResults"), 
          function(object) { slot(object, "repCorSpear") })
setGeneric("repCorSpear<-", function(object, value) { standardGeneric("repCorSpear<-") })
setReplaceMethod("repCorSpear", signature(object="NormalyzerEvaluationResults"), 
                 function(object, value) { 
                     slot(object, "repCorSpear") <- value
                     validObject(object)
                     object
                 })

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
              avgcvmem(ner) <- avgCVPerNormAndReplicates
              featureCVPerMethod(ner) <- calculateFeatureCV(methodList)
              avgcvmempdiff(ner) <- calculatePercentageAvgDiffInMat(avgCVPerNormAndReplicates)
              
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

              avgmadmem(ner) <- avgmadmem
              avgmadmempdiff(ner) <- avgMadPdiff
              
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
              avgvarmem(ner) <- avgVarianceMat
              avgvarmempdiff(ner) <- calculatePercentageAvgDiffInMat(avgVarianceMat)
              
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
              
              anovaP(ner) <- anovaPValsWithNA
              lowVarFeaturesCVs(ner) <- lowVarFeaturesCVs
              lowVarFeaturesCVsPercDiff(ner) <- lowVarFeaturesCVsPercDiff
              
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
            
            repCorPear(ner) <- calculateSummarizedCorrelationVector(
                methodlist, 
                allReplicateGroups, 
                sampleGroupsWithReplicates, "pearson"
            )
            repCorSpear(ner) <- calculateSummarizedCorrelationVector(
                methodlist, 
                allReplicateGroups, 
                sampleGroupsWithReplicates, 
                "spearman"
            )
            
            ner
          }
)


