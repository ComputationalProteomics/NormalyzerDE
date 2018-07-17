#' Statistical dataset representation
#' 
#' @slot annotMat Matrix containing annotation information
#' @slot dataMat Matrix containing (normalized) expression data
#' @slot filteredDataMat Filtered matrix with low-count rows removed
#' @slot designDf Data frame containing design conditions
#' @slot filteringContrast Vector showing which entries are filtered (due to low count)
#' @slot pairwiseCompsP List with P-values for pairwise comparisons
#' @slot pairwiseCompsFdr List with FDR-values for pairwise comparisons
#' @slot pairwiseCompsAve List with average expression values
#' @slot pairwiseCompsFold List with log2 fold-change values for pairwise comparisons
#' @export
NormalyzerStatistics <- setClass("NormalyzerStatistics",
                                 representation(
                                     annotMat = "matrix",
                                     dataMat = "matrix",
                                     filteredDataMat = "matrix",
                                     designDf = "data.frame",
                                     filteringContrast = "vector",
                                     
                                     pairwiseCompsP = "list",
                                     pairwiseCompsFdr = "list",
                                     pairwiseCompsAve = "list",
                                     pairwiseCompsFold = "list"
                                 ),
                                 prototype=prototype(
                                     annotMat=NULL,
                                     dataMat=NULL,
                                     designDf=NULL
                                 ))

#' Main method for contrast calculation
#'
#' @param nst Results evaluation object.
#' @param comparisons String with comparisons for contrasts.
#' @param condCol Column name in design matrix containing condition information.
#' @param batchCol Column name in design matrix containing batch information.
#' @param splitter Character dividing contrast conditions.
#' @param type Type of statistical test (Limma or ANOVA).
#' @return nst Statistics object with statistical measures calculated
#' @rdname calculateContrasts
#' @export
setGeneric(name="calculateContrasts", 
           function(nst, comparisons, condCol, batchCol=NULL, splitter="-", type="limma") standardGeneric("calculateContrasts"))

#' @rdname calculateContrasts
setMethod(f="calculateContrasts", 
          signature=c("NormalyzerStatistics"),
          function(nst, comparisons, condCol, batchCol=NULL, splitter="-", type="limma") {
              
              sampleReplicateGroups <- nst@designDf[, condCol]
              sampleReplicateGroupsStrings <- as.character(nst@designDf[, condCol])
              statMeasures <- c("P", "FDR", "Ave", "Fold")
              compLists <- setupStatMeasureLists(statMeasures, comparisons)
              
              naFilterContrast <- nst@filteringContrast
              dataMatNAFiltered <- nst@filteredDataMat
              
              if (is.null(batchCol)) {
                  Variable <- as.factor(nst@designDf[, condCol])
                  model <- ~0+Variable
              }
              else {
                  if (type != "limma") {
                      stop(paste("Batch compensation only compatible with Limma, got:", type))
                  }
                  Variable <- as.factor(nst@designDf[, condCol])
                  Batch <- as.factor(nst@designDf[, batchCol])
                  model <- ~0+Variable+Batch
              }
              
              if (type == "limma") {
                  limmaDesign <- stats::model.matrix(model)
                  limmaFit <- limma::lmFit(dataMatNAFiltered, limmaDesign)
              }
              
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
                  
                  if (type == "welch") {
                      compLists <- calculateWelch(compLists, dataMatNAFiltered, naFilterContrast, s1cols, s2cols, comp)
                  }
                  else if (type == "limma") {
                      compLists <- calculateLimmaContrast(compLists, dataMatNAFiltered, naFilterContrast, limmaDesign, limmaFit, level1, level2, comp)
                  }
                  else {
                      stop(paste("Unknown statistics type:", type))
                  }
              }
              
              nst@pairwiseCompsP <- compLists[["P"]]
              nst@pairwiseCompsFdr <- compLists[["FDR"]]
              nst@pairwiseCompsAve <- compLists[["Ave"]]
              nst@pairwiseCompsFold <- compLists[["Fold"]]
              nst
          })

setupStatMeasureLists <- function(measures, comparisons) {
    
    compLists <- list()
    for (measure in measures) {
        compLists[[measure]] <- list()
    }
    
    for (measure in measures) {
        for (comp in comparisons) {
            compLists[[measure]][[comp]] <- c()
        }
    }
    compLists
}

calculateWelch <- function(compLists, dataMat, naFilterContrast, s1cols, s2cols, comp) {
    
    doTTest <- function(c1Vals, c2Vals, default=NA) {
        if (length(stats::na.omit(c1Vals)) > 1 && length(stats::na.omit(c2Vals)) > 1) {
            tryCatch(stats::t.test(c1Vals, c2Vals)[[3]], error=function(x) NA)
        }
        else {
            NA
        }
    }
    
    welchPValCol <- apply(dataMat, 1, 
                          function(row) doTTest(row[s1cols], row[s2cols], default=NA))
    welchFDRCol <- stats::p.adjust(welchPValCol, method="BH")

    compLists[["P"]][[comp]][naFilterContrast] <- welchPValCol
    compLists[["FDR"]][[comp]][naFilterContrast] <- welchFDRCol
    compLists[["Ave"]][[comp]][naFilterContrast] <- apply(dataMat, 1, mean)
    compLists[["Fold"]][[comp]][naFilterContrast] <- apply(dataMat, 1,
                                                         function(row) {
                                                             mean(row[s1cols]) - mean(row[s2cols])
                                                         })
    compLists
}

calculateANOVAContrast <- function(compLists, dataMat, naFilterContrast, s1cols, s2cols, comp) {
    
    doTTest <- function(c1Vals, c2Vals, default=NA) {
        if (length(stats::na.omit(c1Vals)) > 1 && length(stats::na.omit(c2Vals)) > 1) {
            tryCatch(stats::t.test(c1Vals, c2Vals)[[3]], error=function(x) NA)
        }
        else {
            NA
        }
    }
    
    welchPValCol <- apply(dataMat, 1, 
                          function(row) doTTest(row[s1cols], row[s2cols], default=NA))
    welchFDRCol <- stats::p.adjust(welchPValCol, method="BH")
    
    compLists[["P"]][[comp]][naFilterContrast] <- welchPValCol
    compLists[["FDR"]][[comp]][naFilterContrast] <- welchFDRCol
    compLists[["Ave"]][[comp]][naFilterContrast] <- apply(dataMat, 1, mean)
    compLists[["Fold"]][[comp]][naFilterContrast] <- apply(dataMat, 1,
                                                           function(row) {
                                                               mean(row[s1cols]) - mean(row[s2cols])
                                                           })
    compLists
}

calculateLimmaContrast <- function(compLists, dataMatNAFiltered, naFilterContrast, limmaDesign, limmaFit, level1, level2, comp) {

    myContrast <- paste0("Variable", level1, "-", "Variable", level2)
    contrastMatrix <- limma::makeContrasts(contrasts=c(myContrast), levels=limmaDesign)
    fitContrasts <- limma::contrasts.fit(limmaFit, contrastMatrix)
    fitBayes <- limma::eBayes(fitContrasts)
    limmaTable <- limma::topTable(fitBayes, coef=1, number=Inf)
    limmaTable <- limmaTable[rownames(dataMatNAFiltered), ]
    
    compLists[["P"]][[comp]][naFilterContrast] <- limmaTable$P.Value
    compLists[["FDR"]][[comp]][naFilterContrast] <- limmaTable$adj.P.Val
    compLists[["Ave"]][[comp]][naFilterContrast] <- limmaTable$AveExpr
    compLists[["Fold"]][[comp]][naFilterContrast] <- limmaTable$logFC
    compLists
}

#' Calculate Limma comparisons between target samples
#'
#' @param nst NormalyzerDE statistics object
#' @param comparisons Vector containing statistical contrasts.
#' @param condCol Column name in design matrix containing condition values
#' @param splitter String splitting the comparisons into condition levels
#' @return None
#' @rdname calculatePairwiseComparisonsLimma
#' @keywords internal
setGeneric(name="calculatePairwiseComparisonsLimma", 
           function(nst, comparisons, condCol, splitter="-") standardGeneric("calculatePairwiseComparisonsLimma"))

#' @rdname calculatePairwiseComparisonsLimma
setMethod(f="calculatePairwiseComparisonsLimma", 
          signature=c("NormalyzerStatistics", "character", "character"),
          function(nst, comparisons, condCol, splitter="-") {
              
              sampleReplicateGroups <- nst@designDf[, condCol]
              sampleReplicateGroupsStrings <- as.character(nst@designDf[, condCol])
              
              pairwisePList <- list()
              pairwiseFDRList <- list()
              pairwiseAveList <- list()
              pairwiseFoldList <- list()
              
              for (comp in comparisons) {
                  pairwisePList[[comp]] <- c()
                  pairwiseFDRList[[comp]] <- c()
                  pairwiseAveList[[comp]] <- c()
                  pairwiseFoldList[[comp]] <- c()
              }
              
              dataMat <- nst@dataMat
              rownames(dataMat) <- seq_len(nrow(dataMat))
              naFilterContrast <- getRowNAFilterContrast(dataMat, sampleReplicateGroups)
              dataMatNAFiltered <- dataMat[naFilterContrast,]
              
              # Limma specific
              Variable <- as.factor(nst@designDf[, condCol])
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
                  fitBayes <- limma::eBayes(fitContrasts)
                  limmaTable <- limma::topTable(fitBayes, coef=1, number=Inf)
                  limmaTable <- limmaTable[rownames(dataMatNAFiltered), ]

                  pairwisePList[[comp]][naFilterContrast] <- limmaTable$P.Value
                  pairwiseFDRList[[comp]][naFilterContrast] <- limmaTable$adj.P.Val
                  pairwiseAveList[[comp]][naFilterContrast] <- limmaTable$AveExpr
                  pairwiseFoldList[[comp]][naFilterContrast] <- limmaTable$logFC
              }
              
              nst@pairwiseCompsP <- pairwisePList
              nst@pairwiseCompsFdr <- pairwiseFDRList
              nst@pairwiseCompsAve <- pairwiseAveList
              nst@pairwiseCompsFold <- pairwiseFoldList
              nst
          }
)




