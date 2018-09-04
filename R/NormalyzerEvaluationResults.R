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
#' @slot anovaP ANOVA calculated p-values
#' @slot repCorPear Within group Pearson correlations
#' @slot repCorSpear Within group Spearman correlations
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
                                               
                                               anovaP = "matrix",
                                               
                                               repCorPear = "matrix",
                                               repCorSpear = "matrix"
                                           ))

#' Constructor for NormalyzerEvaluationResults
#' 
#' @param nr NormalyzerResults object
#' @return nds Generated NormalyzerEvaluationResults instance
#' @export
NormalyzerEvaluationResults <- function (nr) {

               nds <- nds(nr)
               sampleReplicateGroups <- sampleReplicateGroups(nds) 
               methodList <- normalizations(nr)
               sampleGroupsWithReplicates <- samplesGroupsWithReplicates(nds)
               singleReplicateRun <- singleReplicateRun(nds)
               
               # Calculate CV related measures               
               avgCVPerNormAndReplicates <- calculateReplicateCV(methodList, sampleReplicateGroups)
               avgcvmem <- avgCVPerNormAndReplicates
               featureCVPerMethod <- calculateFeatureCV(methodList)
               avgcvmempdiff <- calculatePercentageAvgDiffInMat(avgCVPerNormAndReplicates)
               
               if (!singleReplicateRun) {

                   # MAD
                   avgmadmem <- calculateAvgMadMem(methodList, sampleReplicateGroups)
                   avgmadmempdiff<- calculatePercentageAvgDiffInMat(avgmadmem)
                   
                   # Variance measures
                   avgVarianceMat <- calculateAvgReplicateVariation(methodList, sampleReplicateGroups)
                   avgvarmem <- avgVarianceMat
                   avgvarmempdiff <- calculatePercentageAvgDiffInMat(avgVarianceMat)
                   
                   # Significance measures
                   anovaPValsWithNAMat <- calculateANOVAPValues(methodList, sampleReplicateGroups, categoricalANOVA=TRUE)
                   validPValuesContrast <- !is.na(anovaPValsWithNAMat[, 1])
                   
                   log2AnovaFDR <- rep(NA, length(anovaPValsWithNAMat[, 1]))
                   log2AnovaFDR[validPValuesContrast] <- stats::p.adjust(
                       anovaPValsWithNAMat[, 1][validPValuesContrast],
                       method="BH")
                   
                   lowVarFeaturesCVs <- findLowlyVariableFeaturesCVs(log2AnovaFDR, methodList)
                   lowVarFeaturesCVsPercDiff <- vapply(
                       seq_along(lowVarFeaturesCVs),
                       function(sampleIndex) {
                           (lowVarFeaturesCVs[sampleIndex] * 100) / lowVarFeaturesCVs[1]
                       },
                       0
                   )
                   
                   anovaP <- anovaPValsWithNAMat
                   lowVarFeaturesCVs <- lowVarFeaturesCVs
                   lowVarFeaturesCVsPercDiff <- lowVarFeaturesCVsPercDiff
               }
               
               # Correlation measures
               repCorPear <- calculateSummarizedCorrelationVector(
                   methodList,
                   sampleReplicateGroups,
                   sampleGroupsWithReplicates, 
                   "pearson"
               )
               
               repCorSpear <- calculateSummarizedCorrelationVector(
                   methodList,
                   sampleReplicateGroups,
                   sampleGroupsWithReplicates,
                   "spearman"
               )
               
               object <- new(
                   "NormalyzerEvaluationResults",
                   avgcvmem=avgcvmem,
                   avgcvmempdiff=avgcvmempdiff,
                   featureCVPerMethod=featureCVPerMethod,
                   avgmadmem=avgmadmem,
                   avgmadmempdiff=avgmadmempdiff,
                   avgvarmem=avgvarmem,
                   avgvarmempdiff=avgvarmempdiff,
                   lowVarFeaturesCVs=lowVarFeaturesCVs,
                   lowVarFeaturesCVsPercDiff=lowVarFeaturesCVsPercDiff,
                   anovaP=anovaP,
                   repCorPear=repCorPear,
                   repCorSpear=repCorSpear
               )
               
               return (object)
           }

setGeneric("avgcvmem", function(object) { standardGeneric("avgcvmem") })
setMethod("avgcvmem", signature(object="NormalyzerEvaluationResults"),
          function(object) { slot(object, "avgcvmem") })

setGeneric("avgcvmempdiff", function(object) { standardGeneric("avgcvmempdiff") })
setMethod("avgcvmempdiff", signature(object="NormalyzerEvaluationResults"),
          function(object) { slot(object, "avgcvmempdiff") })

setGeneric("featureCVPerMethod", function(object) { standardGeneric("featureCVPerMethod") })
setMethod("featureCVPerMethod", signature(object="NormalyzerEvaluationResults"),
          function(object) { slot(object, "featureCVPerMethod") })

setGeneric("avgmadmem", function(object) { standardGeneric("avgmadmem") })
setMethod("avgmadmem", signature(object="NormalyzerEvaluationResults"),
          function(object) { slot(object, "avgmadmem") })

setGeneric("avgmadmempdiff", function(object) { standardGeneric("avgmadmempdiff") })
setMethod("avgmadmempdiff", signature(object="NormalyzerEvaluationResults"),
          function(object) { slot(object, "avgmadmempdiff") })

setGeneric("avgvarmem", function(object) { standardGeneric("avgvarmem") })
setMethod("avgvarmem", signature(object="NormalyzerEvaluationResults"),
          function(object) { slot(object, "avgvarmem") })

setGeneric("avgvarmempdiff", function(object) { standardGeneric("avgvarmempdiff") })
setMethod("avgvarmempdiff", signature(object="NormalyzerEvaluationResults"),
          function(object) { slot(object, "avgvarmempdiff") })

setGeneric("lowVarFeaturesCVs", function(object) { standardGeneric("lowVarFeaturesCVs") })
setMethod("lowVarFeaturesCVs", signature(object="NormalyzerEvaluationResults"),
          function(object) { slot(object, "lowVarFeaturesCVs") })

setGeneric("lowVarFeaturesCVsPercDiff", function(object) { standardGeneric("lowVarFeaturesCVsPercDiff") })
setMethod("lowVarFeaturesCVsPercDiff", signature(object="NormalyzerEvaluationResults"),
          function(object) { slot(object, "lowVarFeaturesCVsPercDiff") })

setGeneric("anovaP", function(object) { standardGeneric("anovaP") })
setMethod("anovaP", signature(object="NormalyzerEvaluationResults"),
          function(object) { slot(object, "anovaP") })

setGeneric("repCorPear", function(object) { standardGeneric("repCorPear") })
setMethod("repCorPear", signature(object="NormalyzerEvaluationResults"),
          function(object) { slot(object, "repCorPear") })

setGeneric("repCorSpear", function(object) { standardGeneric("repCorSpear") })
setMethod("repCorSpear", signature(object="NormalyzerEvaluationResults"),
          function(object) { slot(object, "repCorSpear") })
