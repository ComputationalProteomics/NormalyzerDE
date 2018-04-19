#' S4 class to represent dataset information
#' 
#' @export
NormalyzerStatistics <- setClass("NormalyzerStatistics",
                                 representation(
                                     annotMat = "matrix",
                                     dataMat = "matrix",
                                     designDf = "data.frame",
                                     
                                     # sampleCol = "character",
                                     # conditionCol = "character",
                                     # batchCol = "character",
                                     
                                     pairwiseComps = "list",
                                     pairwiseCompsFdr = "list"
                                 ),
                                 prototype=prototype(
                                     annotMat=NULL,
                                     dataMat=NULL,
                                     designDf=NULL
                                 ))

#' Calculate Welch t-test comparisons between target samples
#'
#' @param ner Results evaluation object.
#' @param nr Results object.
#' @param comparisons Vector containing statistical contrasts.
#' @return None
#' @rdname calculatePairwiseComparisons
setGeneric(name="calculatePairwiseComparisons", 
           function(nst, comparisons, condCol, splitter="-") standardGeneric("calculatePairwiseComparisons"))

#' @rdname calculatePairwiseComparisons
setMethod(f="calculatePairwiseComparisons", 
          signature=c("NormalyzerStatistics", "character"),
          function(nst, comparisons, condCol, splitter="-") {
              
              sampleReplicateGroups <- nst@designDf[, condCol]
              sampleReplicateGroupsStrings <- as.character(nst@designDf[, condCol])
              
              pairwisePList <- list()
              pairwiseFDRList <- list()
              for (comp in comparisons) {
                  pairwisePList[[comp]] <- c()
                  pairwiseFDRList[[comp]] <- c()
              }
              
              do_t_test <- function(c1_vals, c2_vals, default=NA) {
                  if (length(stats::na.omit(c1_vals)) > 1 && length(stats::na.omit(c2_vals)) > 1) {
                      tryCatch(stats::t.test(c1_vals, c2_vals)[[3]], error=function(x) NA)
                  }
                  else {
                      NA
                  }
              }
              
              processedDataMatrix <- nst@dataMat
              naFilterContrast <- getRowNAFilterContrast(processedDataMatrix, sampleReplicateGroups)
              dataStoreReplicateNAFiltered <- processedDataMatrix[naFilterContrast,]
              
              for (comp in comparisons) {
                  
                  compSplit <- unlist(strsplit(comp, splitter))
                  
                  if (length(compSplit) != 2) {
                      stop(paste("Comparison should be in format cond1-cond2, here the split product was:", paste(compSplit, collapse=" ")))
                  }
                  
                  level1 <- compSplit[1]
                  level2 <- compSplit[2]
                  
                  if (length(sampleReplicateGroupsStrings %in% level1) == 0) {
                      stop(paste("No samples matching condition", level1, "found in conditions:", paste(sampleReplicateGroupsStrings, collapse=" ")))
                  }
                  
                  if (length(sampleReplicateGroupsStrings %in% level2) == 0) {
                      stop(paste("No samples matching condition", level2, "found in conditions:", paste(sampleReplicateGroupsStrings, collapse=" ")))
                  }
                  
                  s1cols <- which(sampleReplicateGroupsStrings %in% level1)
                  s2cols <- which(sampleReplicateGroupsStrings %in% level2)
                  
                  welchPValCol <- apply(dataStoreReplicateNAFiltered, 1, 
                                        function(row) do_t_test(row[s1cols], row[s2cols], default=NA))
                  welchFDRCol <- stats::p.adjust(welchPValCol, method="BH")
                  
                  pairwisePList[[comp]][naFilterContrast] <- welchPValCol
                  pairwiseFDRList[[comp]][naFilterContrast] <- welchFDRCol
              }

              nst@pairwiseComps <- pairwisePList
              nst@pairwiseCompsFdr <- pairwiseFDRList
              nst
          }
)

#' Calculate Limma comparisons between target samples
#'
#' @param ner Results evaluation object.
#' @param nr Results object.
#' @param comparisons Vector containing statistical contrasts.
#' @return None
#' @rdname calculatePairwiseComparisonsLimma
setGeneric(name="calculatePairwiseComparisonsLimma", 
           function(nst, comparisons, condCol, splitter="-", robustLimma) standardGeneric("calculatePairwiseComparisonsLimma"))

#' @rdname calculatePairwiseComparisonsLimma
setMethod(f="calculatePairwiseComparisonsLimma", 
          signature=c("NormalyzerStatistics", "character", "character"),
          function(nst, comparisons, condCol, splitter="-", robustLimma=FALSE) {
              
              Variable <- as.factor(nst@designDf[, condCol])
              model <- ~0+Variable
              
              sampleReplicateGroups <- nst@designDf[, condCol]
              sampleReplicateGroupsStrings <- as.character(nst@designDf[, condCol])
              
              pairwisePList <- list()
              pairwiseFDRList <- list()
              for (comp in comparisons) {
                  pairwisePList[[comp]] <- c()
                  pairwiseFDRList[[comp]] <- c()
              }
              
              dataMat <- nst@dataMat
              rownames(dataMat) <- 1:nrow(dataMat)
              naFilterContrast <- getRowNAFilterContrast(dataMat, sampleReplicateGroups)
              dataMatNAFiltered <- dataMat[naFilterContrast,]
              
              model <- ~0+Variable
              limmaDesign <- stats::model.matrix(model)
              limmaFit <- limma::lmFit(dataMatNAFiltered, limmaDesign)
              
              for (comp in comparisons) {
                  
                  compSplit <- unlist(strsplit(comp, splitter))
                  
                  if (length(compSplit) != 2) {
                      stop(paste("Comparison should be in format cond1-cond2, here the split product was:", paste(compSplit, collapse=" ")))
                  }
                  
                  level1 <- compSplit[1]
                  level2 <- compSplit[2]
                  
                  if (length(sampleReplicateGroupsStrings %in% level1) == 0) {
                      stop(paste("No samples matching condition", level1, "found in conditions:", paste(sampleReplicateGroupsStrings, collapse=" ")))
                  }
                  
                  if (length(sampleReplicateGroupsStrings %in% level2) == 0) {
                      stop(paste("No samples matching condition", level2, "found in conditions:", paste(sampleReplicateGroupsStrings, collapse=" ")))
                  }

                  # Limma contrast
                  myContrast <- paste0("Variable", level1, "-", "Variable", level2)
                  contrastMatrix <- limma::makeContrasts(contrasts=c(myContrast), levels=limmaDesign)
                  fitContrasts <- limma::contrasts.fit(limmaFit, contrastMatrix)
                  fitBayes <- limma::eBayes(fitContrasts, robust=robustLimma)
                  limmaTable <- limma::topTable(fitBayes, coef=1, number=Inf)
                  limmaTable <- limmaTable[rownames(dataMatNAFiltered), ]
                  

                  pairwisePList[[comp]][naFilterContrast] <- limmaTable$P.Value
                  pairwiseFDRList[[comp]][naFilterContrast] <- limmaTable$adj.P.Val
              }
              
              nst@pairwiseComps <- pairwisePList
              nst@pairwiseCompsFdr <- pairwiseFDRList
              nst
          }
)




