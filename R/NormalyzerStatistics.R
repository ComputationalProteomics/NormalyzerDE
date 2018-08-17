#' Class representing a dataset for statistical processing in NormalyzerDE
#' 
#' Is initialized with an annotation matrix, a data matrix and a design
#' data frame. This object can subsequently be processed to generate statistical
#' values and in turn used to write a full matrix with additional statistical
#' information as well as a graphical report of the comparisons.
#' 
#' @slot annotMat Matrix containing annotation information
#' @slot dataMat Matrix containing (normalized) expression data
#' @slot filteredDataMat Filtered matrix with low-count rows removed
#' @slot designDf Data frame containing design conditions
#' @slot filteringContrast Vector showing which entries are filtered 
#'   (due to low count)
#' @slot pairwiseCompsP List with P-values for pairwise comparisons
#' @slot pairwiseCompsFdr List with FDR-values for pairwise comparisons
#' @slot pairwiseCompsAve List with average expression values
#' @slot pairwiseCompsFold List with log2 fold-change values for pairwise 
#'   comparisons
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

#' Performs statistical comparisons between the supplied conditions.
#' It uses the design matrix and data matrix in the supplied 
#' NormalyzerStatistics object. A column is supplied specifying which of the
#' columns in the design matrix that is used for deciding the sample groups.
#' The comparisons vector specifies which pairwise comparisons between
#' condition levels that are to be calculated.
#' 
#' Optionally, a batch column can be specified allowing compensation for
#' covariate variation in the statistical model. This is only compatible
#' with a Limma-based statistical analysis.
#'
#' @param nst Results evaluation object.
#' @param comparisons String with comparisons for contrasts.
#' @param condCol Column name in design matrix containing condition information.
#' @param batchCol Column name in design matrix containing batch information.
#' @param splitter Character dividing contrast conditions.
#' @param type Type of statistical test (Limma or welch).
#' @return nst Statistics object with statistical measures calculated
#' @rdname calculateContrasts
#' @export
#' @examples
#' data("example_stat_data")
#' data("example_design")
#' nst <- setupStatisticsObject(example_design, example_stat_data)
#' results <- calculateContrasts(nst, c("1-2", "2-3"), "group")
#' resultsBatch <- calculateContrasts(nst, c("1-2", "2-3"), "group", "batch")
setGeneric(name="calculateContrasts", 
           function(nst, comparisons, condCol, batchCol=NULL, splitter="-", 
                    type="limma") standardGeneric("calculateContrasts"))

#' @rdname calculateContrasts
setMethod(f="calculateContrasts", 
          signature=c("NormalyzerStatistics"),
          function(nst, comparisons, condCol, batchCol=NULL, splitter="-", 
                   type="limma") {
              
              sampleReplicateGroups <- nst@designDf[, condCol]
              sampleReplicateGroupsStrings <- as.character(nst@designDf[, condCol])
              statMeasures <- c("P", "FDR", "Ave", "Fold")

              naFilterContrast <- nst@filteringContrast
              dataMatNAFiltered <- nst@filteredDataMat
              
              if (is.null(batchCol)) {
                  Variable <- as.factor(nst@designDf[, condCol])
                  model <- ~0+Variable
              }
              else {
                  if (type != "limma") {
                      stop(
                          "Batch compensation only compatible with Limma, got: ", 
                          type
                      )
                  }
                  Variable <- as.factor(nst@designDf[, condCol])
                  Batch <- as.factor(nst@designDf[, batchCol])
                  model <- ~0+Variable+Batch
              }
              
              if (type == "limma") {
                  limmaDesign <- stats::model.matrix(model)
                  limmaFit <- limma::lmFit(dataMatNAFiltered, limmaDesign)
              }
              
              compLists <- list()
              for (statMeasure in statMeasures) {
                  compLists[[statMeasure]] <- list()
              }

              for (comp in comparisons) {
                  
                  compSplit <- unlist(strsplit(comp, splitter))
                  
                  if (length(compSplit) != 2) {
                      stop("Comparison should be in format cond1-cond2 ", 
                           "here the split product was: ", 
                           paste(compSplit, collapse=" "))
                  }
                  
                  level1 <- compSplit[1]
                  level2 <- compSplit[2]
                  
                  if (length(sampleReplicateGroupsStrings %in% level1) == 0) {
                      stop("No samples matching condition ", 
                           level1, 
                           " found in conditions: ", 
                           paste(sampleReplicateGroupsStrings, collapse=" "))
                  }
                  
                  if (length(sampleReplicateGroupsStrings %in% level2) == 0) {
                      stop("No samples matching condition ", 
                           level2, " found in conditions: ", 
                           paste(sampleReplicateGroupsStrings, collapse=" "))
                  }
                  
                  if (type == "welch") {
                      statResults <- calculateWelch(
                          dataMatNAFiltered, 
                          sampleReplicateGroupsStrings, 
                          c(level1, level2))
                  }
                  else if (type == "limma") {
                      statResults <- calculateLimmaContrast(
                          dataMatNAFiltered, 
                          limmaDesign, 
                          limmaFit, 
                          c(level1, level2))
                  }
                  else {
                      stop("Unknown statistics type: ", type)
                  }
                  
                  for (statMeasure in statMeasures) {
                      compLists[[statMeasure]][[comp]] <- c()
                      compLists[[statMeasure]][[comp]][naFilterContrast] <- statResults[[statMeasure]]
                  }
              }
              
              nst@pairwiseCompsP <- compLists[["P"]]
              nst@pairwiseCompsFdr <- compLists[["FDR"]]
              nst@pairwiseCompsAve <- compLists[["Ave"]]
              nst@pairwiseCompsFold <- compLists[["Fold"]]
              nst
          })

calculateWelch <- function(dataMat, groupHeader, levels) {
    
    s1cols <- which(groupHeader %in% levels[1])
    s2cols <- which(groupHeader %in% levels[2])
    
    doTTest <- function(c1Vals, c2Vals, default=NA) {
        if (length(stats::na.omit(c1Vals)) > 1 && 
            length(stats::na.omit(c2Vals)) > 1) {
            tryCatch(stats::t.test(c1Vals, c2Vals)[[3]], error=function(x) NA)
        }
        else {
            NA
        }
    }
    
    welchPValCol <- apply(
        dataMat, 1, 
        function(row) doTTest(row[s1cols], row[s2cols], default=NA))
    welchFDRCol <- stats::p.adjust(welchPValCol, method="BH")

    statResults <- list()

    statResults[["P"]] <- welchPValCol
    statResults[["FDR"]] <- welchFDRCol
    statResults[["Ave"]] <- apply(dataMat, 1, mean, na.rm=TRUE)
    statResults[["Fold"]] <- apply(
        dataMat, 
        1,
        function(row) { mean(row[s1cols], na.rm=TRUE) - mean(row[s2cols], na.rm=TRUE) })

    statResults
}

calculateLimmaContrast <- function(dataMat, limmaDesign, limmaFit, levels) {

    myContrast <- paste0("Variable", levels[1], "-", "Variable", levels[2])
    contrastMatrix <- limma::makeContrasts(
        contrasts=c(myContrast), 
        levels=limmaDesign)
    fitContrasts <- limma::contrasts.fit(limmaFit, contrastMatrix)
    fitBayes <- limma::eBayes(fitContrasts)
    limmaTable <- limma::topTable(fitBayes, coef=1, number=Inf, sort.by="none")

    statResults <- list()
    statResults[["P"]] <- limmaTable$P.Value
    statResults[["FDR"]] <- limmaTable$adj.P.Val
    statResults[["Ave"]] <- limmaTable$AveExpr
    statResults[["Fold"]] <- limmaTable$logFC
    statResults
}

