#' S4 class to represent dataset information
#' 
#' @export
NormalyzerStatistics <- setClass("NormalyzerResults",
                                 representation(),
                                 prototype=prototype())

#' Calculate Welch t-test comparisons between target samples
#'
#' @param ner Results evaluation object.
#' @param nr Results object.
#' @param comparisons Vector containing statistical contrasts.
#' @return None
#' @rdname calculatePairwiseComparisons
setGeneric(name="calculatePairwiseComparisons", 
           function(ner, nr, comparisons) standardGeneric("calculatePairwiseComparisons"))

#' @rdname calculatePairwiseComparisons
setMethod(f="calculatePairwiseComparisons", 
          signature=c("NormalizationEvaluationResults", "NormalyzerResults", "character"),
          function(ner, nr, comparisons) {
              
              # Setup
              sampleReplicateGroups <- nr@nds@sampleReplicateGroups
              methodCount <- length(getUsedMethodNames(nr))
              methodList <- getNormalizationMatrices(nr)
              
              indexList <- getIndexList(sampleReplicateGroups)
              
              pairwisePList <- list()
              pairwiseFDRList <- list()
              for (comp in comparisons) {
                  pairwisePList[[comp]] <- matrix(NA, ncol=methodCount, nrow=nrow(methodList[[1]]))
                  pairwiseFDRList[[comp]] <- matrix(NA, ncol=methodCount, nrow=nrow(methodList[[1]]))
              }
              
              do_t_test <- function(c1_vals, c2_vals, default=NA) {
                  if (length(stats::na.omit(c1_vals)) > 1 && length(stats::na.omit(c2_vals)) > 1) {
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
                      
                      sample1 <- substr(comp, 1, 1)
                      sample2 <- substr(comp, 2, 2)
                      
                      s1cols <- indexList[[sample1]]
                      s2cols <- indexList[[sample2]]
                      
                      welchPValCol <- apply(dataStoreReplicateNAFiltered, 1, function(row) do_t_test(row[s1cols], row[s2cols], default=NA))
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