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

#' Calculate contrast
#'
#' @param ner Results evaluation object.
#' @param nr Results object.
#' @param comparisons Vector containing statistical contrasts.
#' @return None
#' @rdname calculateContrasts
setGeneric(name="calculateContrasts", 
           function(nst, comparisons, condCol, splitter="-", type="anova") standardGeneric("calculateContrasts"))

#' @rdname calculateContrasts
setMethod(f="calculateContrasts", 
          signature=c("NormalyzerStatistics", "character"),
          function(nst, comparisons, condCol, splitter="-", type="anova") {
              
              sampleReplicateGroups <- nst@designDf[, condCol]
              sampleReplicateGroupsStrings <- as.character(nst@designDf[, condCol])
              statMeasures <- c("P", "FDR", "Ave", "Fold")
              compLists <- setupStatMeasureLists(statMeasures, comparisons)
              
              processedDataMatrix <- nst@dataMat
              naFilterContrast <- getRowNAFilterContrast(processedDataMatrix, sampleReplicateGroups)
              dataMatNAFiltered <- processedDataMatrix[naFilterContrast, ]
              
              if (type == "limma") {
                  
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
                      compLists <- calculateWelch(dataMatNAFiltered, naFilterContrast, s1cols, s2cols, comp, compLists)
                  }
                  else if (type == "limma") {
                      compLists <- calculateLimmaContrast
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

calculateWelch <- function(dataMat, naFilterContrast, s1cols, s2cols, comp, compLists) {
    
    do_t_test <- function(c1_vals, c2_vals, default=NA) {
        if (length(stats::na.omit(c1_vals)) > 1 && length(stats::na.omit(c2_vals)) > 1) {
            tryCatch(stats::t.test(c1_vals, c2_vals)[[3]], error=function(x) NA)
        }
        else {
            NA
        }
    }
    
    welchPValCol <- apply(dataMat, 1, 
                          function(row) do_t_test(row[s1cols], row[s2cols], default=NA))
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

calculateLimmaContrast <- function(level1, level2, limmaDesign, limmaFit, dataMatNAFiltered) {
    
    myContrast <- paste0("Variable", level1, "-", "Variable", level2)
    contrastMatrix <- limma::makeContrasts(contrasts=c(myContrast), levels=limmaDesign)
    fitContrasts <- limma::contrasts.fit(limmaFit, contrastMatrix)
    fitBayes <- limma::eBayes(fitContrasts, robust=robustLimma)
    limmaTable <- limma::topTable(fitBayes, coef=1, number=Inf)
    limmaTable <- limmaTable[rownames(dataMatNAFiltered), ]
    
    pairwisePList[["P"]][[comp]][naFilterContrast] <- limmaTable$P.Value
    pairwiseFDRList[["FDR"]][[comp]][naFilterContrast] <- limmaTable$adj.P.Val
    pairwiseAveList[["Ave"]][[comp]][naFilterContrast] <- limmaTable$AveExpr
    pairwiseFoldList[["Fold"]][[comp]][naFilterContrast] <- limmaTable$logFC
    
}

#' #' Calculate Welch t-test comparisons between target samples
#' #'
#' #' @param ner Results evaluation object.
#' #' @param nr Results object.
#' #' @param comparisons Vector containing statistical contrasts.
#' #' @return None
#' #' @rdname calculatePairwiseComparisonsWelchT
#' setGeneric(name="calculatePairwiseComparisonsWelchT", 
#'            function(nst, comparisons, condCol, splitter="-") standardGeneric("calculatePairwiseComparisonsWelchT"))
#' 
#' #' @rdname calculatePairwiseComparisons
#' setMethod(f="calculatePairwiseComparisonsWelchT", 
#'           signature=c("NormalyzerStatistics", "character"),
#'           function(nst, comparisons, condCol, splitter="-") {
#'               
#'               # sampleReplicateGroups <- nst@designDf[, condCol]
#'               # sampleReplicateGroupsStrings <- as.character(nst@designDf[, condCol])
#'               # 
#'               # pairwisePList <- list()
#'               # pairwiseFDRList <- list()
#'               # pairwiseAveList <- list()
#'               # pairwiseFoldList <- list()
#'               # 
#'               # for (comp in comparisons) {
#'               #     pairwisePList[[comp]] <- c()
#'               #     pairwiseFDRList[[comp]] <- c()
#'               #     pairwiseAveList[[comp]] <- c()
#'               #     pairwiseFoldList[[comp]] <- c()
#'               # }
#'               
#'               do_t_test <- function(c1_vals, c2_vals, default=NA) {
#'                   if (length(stats::na.omit(c1_vals)) > 1 && length(stats::na.omit(c2_vals)) > 1) {
#'                       tryCatch(stats::t.test(c1_vals, c2_vals)[[3]], error=function(x) NA)
#'                   }
#'                   else {
#'                       NA
#'                   }
#'               }
#'               
#'               # processedDataMatrix <- nst@dataMat
#'               # naFilterContrast <- getRowNAFilterContrast(processedDataMatrix, sampleReplicateGroups)
#'               # dataMatNAFiltered <- processedDataMatrix[naFilterContrast, ]
#'               
#'               for (comp in comparisons) {
#'                   
#'                   # compSplit <- unlist(strsplit(comp, splitter))
#'                   # 
#'                   # if (length(compSplit) != 2) {
#'                   #     stop(paste("Comparison should be in format cond1-cond2, here the split product was:", paste(compSplit, collapse=" ")))
#'                   # }
#'                   # 
#'                   # level1 <- compSplit[1]
#'                   # level2 <- compSplit[2]
#'                   # 
#'                   # if (length(sampleReplicateGroupsStrings %in% level1) == 0) {
#'                   #     stop(paste("No samples matching condition", level1, "found in conditions:", paste(sampleReplicateGroupsStrings, collapse=" ")))
#'                   # }
#'                   # 
#'                   # if (length(sampleReplicateGroupsStrings %in% level2) == 0) {
#'                   #     stop(paste("No samples matching condition", level2, "found in conditions:", paste(sampleReplicateGroupsStrings, collapse=" ")))
#'                   # }
#'                   # 
#'                   # s1cols <- which(sampleReplicateGroupsStrings %in% level1)
#'                   # s2cols <- which(sampleReplicateGroupsStrings %in% level2)
#'                   
#'                   welchPValCol <- apply(dataMatNAFiltered, 1, 
#'                                         function(row) do_t_test(row[s1cols], row[s2cols], default=NA))
#'                   welchFDRCol <- stats::p.adjust(welchPValCol, method="BH")
#'                   
#'                   pairwisePList[[comp]][naFilterContrast] <- welchPValCol
#'                   pairwiseFDRList[[comp]][naFilterContrast] <- welchFDRCol
#'                   pairwiseAveList[[comp]][naFilterContrast] <- apply(dataMatNAFiltered, 1, mean)
#'                   pairwiseFoldList[[comp]][naFilterContrast] <- apply(dataMatNAFiltered, 1, 
#'                                                                       function(row) {
#'                                                                           mean(row[s1cols]) - mean(row[s2cols])
#'                                                                       })
#'               }
#' 
#'               # nst@pairwiseCompsP <- pairwisePList
#'               # nst@pairwiseCompsFdr <- pairwiseFDRList
#'               # nst@pairwiseCompsAve <- pairwiseAveList
#'               # nst@pairwiseCompsFold <- pairwiseFoldList
#'               # nst
#'           }
#' )

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
              rownames(dataMat) <- 1:nrow(dataMat)
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
                  fitBayes <- limma::eBayes(fitContrasts, robust=robustLimma)
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




