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
#' @slot contrasts Spot for saving vector of last used contrasts
#' @slot condCol Column containing last used conditions
#' @slot batchCol Column containing last used batch conditions
NormalyzerStatistics <- setClass("NormalyzerStatistics",
                                 slots = c(
                                     annotMat = "matrix",
                                     dataMat = "matrix",
                                     designDf = "data.frame",

                                     pairwiseCompsP = "list",
                                     pairwiseCompsFdr = "list",
                                     pairwiseCompsAve = "list",
                                     pairwiseCompsFold = "list",
                                     pairwiseCompsSig = "list",
                                     
                                     comparisons = "character",
                                     condCol = "character",
                                     batchCol = "numeric"
                                 ))

#' Constructor for NormalyzerStatistics
#' 
#' @param experimentObj Instance of SummarizedExperiment containing matrix
#'   and design information as column data
#' @param logTrans Whether the input data should be log transformed
#' @return nds Generated NormalyzerStatistics instance
#' @export
#' @examples
#' data(example_stat_summarized_experiment)
#' nst <- NormalyzerStatistics(example_stat_summarized_experiment)
NormalyzerStatistics <- function(experimentObj, logTrans=FALSE) { 

              dataMat <- SummarizedExperiment::assay(experimentObj)
              if (logTrans) {
                  dataMat <- log2(dataMat)
              }
              
              annotMat <- SummarizedExperiment::rowData(experimentObj)
              designDf <- SummarizedExperiment::colData(experimentObj)

              nst <- new("NormalyzerStatistics",
                         annotMat=as.matrix(annotMat), 
                         dataMat=as.matrix(dataMat), 
                         designDf=as.data.frame(designDf)
              )

              nst
          }

setGeneric("condCol", function(object) { standardGeneric("condCol") })
setMethod("condCol", signature(object="NormalyzerStatistics"), 
          function(object) { slot(object, "condCol") })
setGeneric("condCol<-", function(object, value) { standardGeneric("condCol<-") })
setReplaceMethod("condCol", signature(object="NormalyzerStatistics"), 
                 function(object, value) { 
                     slot(object, "condCol") <- value
                     validObject(object)
                     object
                 })

setGeneric("batchCol", function(object) { standardGeneric("batchCol") })
setMethod("batchCol", signature(object="NormalyzerStatistics"), 
          function(object) { slot(object, "batchCol") })
setGeneric("batchCol<-", function(object, value) { standardGeneric("batchCol<-") })
setReplaceMethod("batchCol", signature(object="NormalyzerStatistics"), 
                 function(object, value) { 
                     slot(object, "batchCol") <- value
                     validObject(object)
                     object
                 })

setGeneric("comparisons", function(object) { standardGeneric("comparisons") })
setMethod("comparisons", signature(object="NormalyzerStatistics"), 
          function(object) { slot(object, "comparisons") })
setGeneric("comparisons<-", function(object, value) { standardGeneric("comparisons<-") })
setReplaceMethod("comparisons", signature(object="NormalyzerStatistics"), 
                 function(object, value) { 
                     slot(object, "comparisons") <- value
                     validObject(object)
                     object
                 })

setGeneric("annotMat", function(object) { standardGeneric("annotMat") })
setMethod("annotMat", signature(object="NormalyzerStatistics"), 
          function(object) { slot(object, "annotMat") })

setGeneric("dataMat", function(object) { standardGeneric("dataMat") })
setMethod("dataMat", signature(object="NormalyzerStatistics"), 
          function(object) { slot(object, "dataMat") })
setGeneric("dataMat<-", function(object, value) { standardGeneric("dataMat<-") })
setReplaceMethod("dataMat", signature(object="NormalyzerStatistics"), 
                 function(object, value) { 
                     
                     slot(object, "dataMat") <- value
                     validObject(object)
                     object
                 })

setGeneric("filteredDataMat", function(object) { standardGeneric("filteredDataMat") })
setMethod("filteredDataMat", signature(object="NormalyzerStatistics"), 
          function(object) { slot(object, "filteredDataMat") })

setGeneric("designDf", function(object) { standardGeneric("designDf") })
setMethod("designDf", signature(object="NormalyzerStatistics"), 
          function(object) { slot(object, "designDf") })
setGeneric("designDf<-", function(object, value) { standardGeneric("designDf<-") })
setReplaceMethod("designDf", signature(object="NormalyzerStatistics"), 
                 function(object, value) { 
                     slot(object, "designDf") <- value
                     validObject(object)
                     object
                 })

setGeneric("filteringContrast", function(object) { standardGeneric("filteringContrast") })
setMethod("filteringContrast", signature(object="NormalyzerStatistics"), 
          function(object) { slot(object, "filteringContrast") })

setGeneric("pairwiseCompsP", function(object) { standardGeneric("pairwiseCompsP") })
setMethod("pairwiseCompsP", signature(object="NormalyzerStatistics"), 
          function(object) { slot(object, "pairwiseCompsP") })
setGeneric("pairwiseCompsP<-", function(object, value) { standardGeneric("pairwiseCompsP<-") })
setReplaceMethod("pairwiseCompsP", signature(object="NormalyzerStatistics"), 
                 function(object, value) { 
                     slot(object, "pairwiseCompsP") <- value
                     validObject(object)
                     object
                 })

setGeneric("pairwiseCompsFdr", function(object) { standardGeneric("pairwiseCompsFdr") })
setMethod("pairwiseCompsFdr", signature(object="NormalyzerStatistics"), 
          function(object) { slot(object, "pairwiseCompsFdr") })
setGeneric("pairwiseCompsFdr<-", function(object, value) { standardGeneric("pairwiseCompsFdr<-") })
setReplaceMethod("pairwiseCompsFdr", signature(object="NormalyzerStatistics"), 
                 function(object, value) { 
                     slot(object, "pairwiseCompsFdr") <- value
                     validObject(object)
                     object
                 })

setGeneric("pairwiseCompsAve", function(object) { standardGeneric("pairwiseCompsAve") })
setMethod("pairwiseCompsAve", signature(object="NormalyzerStatistics"), 
          function(object) { slot(object, "pairwiseCompsAve") })
setGeneric("pairwiseCompsAve<-", function(object, value) { standardGeneric("pairwiseCompsAve<-") })
setReplaceMethod("pairwiseCompsAve", signature(object="NormalyzerStatistics"), 
                 function(object, value) { 
                     slot(object, "pairwiseCompsAve") <- value
                     validObject(object)
                     object
                 })

setGeneric("pairwiseCompsFold", function(object) { standardGeneric("pairwiseCompsFold") })
setMethod("pairwiseCompsFold", signature(object="NormalyzerStatistics"), 
          function(object) { slot(object, "pairwiseCompsFold") })
setGeneric("pairwiseCompsFold<-", function(object, value) { standardGeneric("pairwiseCompsFold<-") })
setReplaceMethod("pairwiseCompsFold", signature(object="NormalyzerStatistics"), 
                 function(object, value) { 
                     slot(object, "pairwiseCompsFold") <- value
                     validObject(object)
                     object
                 })


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
#' @param comparisons Character vector with pairwise comparisons for contrasts.
#'   Ignored if \code{oneVsRest=TRUE}.
#' @param condCol Column name in design matrix containing condition information.
#' @param batchCol Column name in design matrix containing batch information.
#' @param splitter Character dividing contrast conditions.
#' @param type Type of statistical test (Limma or welch).
#' @param leastRepCount Least replicates in each group to be retained for 
#'   contrast calculations
#' @param impute Whether to impute values
#' @param imputeMinFraction Minimum fraction non-NA values for an analyte in any group to impute in other groups
#' @param subsetByComparison If TRUE, subset data and design to each comparison
#'   before NA-filtering, imputation and model fitting.
#' @param oneVsRest If TRUE, compute one-vs-rest contrasts for each group in
#'   \code{condCol} (or the subset in \code{oneVsRestGroups}).
#' @param oneVsRestGroups Optional character vector specifying which groups in
#'   \code{condCol} to compare against all other samples.
#' @return nst Statistics object with statistical measures calculated
#' @rdname calculateContrasts 
#' @export
#' @examples
#' data(example_stat_summarized_experiment)
#' nst <- NormalyzerStatistics(example_stat_summarized_experiment)
#' results <- calculateContrasts(nst, c("1-2", "2-3"), "group")
#' resultsBatch <- calculateContrasts(nst, c("1-2", "2-3"), "group", batchCol="batch")
#' resultsOneVsRest <- calculateContrasts(nst, condCol="group", oneVsRest=TRUE)
setGeneric(name="calculateContrasts", 
           function(nst, comparisons=NULL, condCol, batchCol=NULL, splitter="-", 
                    type="limma", leastRepCount=1, impute = FALSE, imputeMinFraction=1,
                    subsetByComparison = FALSE, oneVsRest = FALSE, oneVsRestGroups = NULL) standardGeneric("calculateContrasts"))

#' @rdname calculateContrasts
setMethod(f="calculateContrasts", 
          signature=c("NormalyzerStatistics"),
          function(nst, comparisons=NULL, condCol, batchCol=NULL, splitter="-", 
                   type="limma", leastRepCount=1, impute = FALSE, imputeMinFraction=1,
                   subsetByComparison = FALSE, oneVsRest = FALSE, oneVsRestGroups = NULL) {
              
              dataMat <- dataMat(nst)
              designDf <- designDf(nst)

              condCol(nst) <- as.character(designDf[, condCol])

              if (!is.null(batchCol)) {
                  conditionCombs <- paste(designDf[, condCol], designDf[, batchCol], sep="_")
                  batchCol(nst) <- as.factor(designDf[, batchCol])
              }
              else {
                  conditionCombs <- designDf[, condCol]
                  batchCol(nst) <- numeric()
              }
              
              rownames(dataMat) <- seq_len(nrow(dataMat))

              sampleReplicateGroupsStrings <- as.character(designDf[, condCol])
              statMeasures <- c("P", "FDR", "Ave", "Fold")

              setupModelFromDesign <- function(designDf, condCol, batchCol=NULL, type="limma") {
                  if (is.null(batchCol)) {
                      Variable <- as.factor(designDf[, condCol])
                      model <- ~0+Variable
                  }
                  else {
                      if (!(type %in% c("limma", "limma_intensity"))) {
                          stop(
                              "Batch compensation only compatible with Limma, got: ", 
                              type
                          )
                      }
                      Variable <- as.factor(designDf[, condCol])
                      Batch <- as.factor(designDf[, batchCol])
                      model <- ~0+Variable+Batch
                  }
                  model
              }

              compLists <- list()
              for (statMeasure in statMeasures) {
                  compLists[[statMeasure]] <- list()
              }

	              if (oneVsRest) {

	                  restLabel <- chooseOneVsRestLabel(sampleReplicateGroupsStrings)

	                  targetGroups <- if (is.null(oneVsRestGroups)) {
	                      unique(sampleReplicateGroupsStrings)
	                  } else {
                      unique(as.character(oneVsRestGroups))
                  }

                  if (length(targetGroups) == 0) {
                      stop("No groups specified for one-vs-rest comparisons.")
                  }

                  missingGroups <- setdiff(targetGroups, unique(sampleReplicateGroupsStrings))
                  if (length(missingGroups) > 0) {
                      stop(
                          "Some groups in oneVsRestGroups were not found in condCol '",
                          condCol,
                          "': ",
                          paste(missingGroups, collapse=", ")
                      )
                  }

                  comparisonsGenerated <- paste(targetGroups, restLabel, sep=splitter)
                  comparisons(nst) <- comparisonsGenerated

                  for (groupLabel in targetGroups) {

                      compName <- paste(groupLabel, restLabel, sep=splitter)
                      groupHeader <- ifelse(sampleReplicateGroupsStrings %in% groupLabel, groupLabel, restLabel)

                      designDfOVR <- designDf
                      designDfOVR[, condCol] <- groupHeader

                      conditionCombsOVR <- if (!is.null(batchCol)) {
                          paste(groupHeader, designDfOVR[, batchCol], sep="_")
                      } else {
                          groupHeader
                      }

                      dataMatNAFiltered <- filterLowRep(
                          dataMat, 
                          conditionCombsOVR, 
                          leastRep=leastRepCount
                      )

                      if (nrow(dataMatNAFiltered) == 0) {
                          stop(
                              "No rows remained after NA-filtering for one-vs-rest comparison: '",
                              compName,
                              "' (condCol: '", condCol,
                              "', batchCol: '", batchCol,
                              "' (if empty then batchCol is not specified))\n",
                              "Consider whether you can reduce the 'leastRepCount' setting which sets the lower limit ",
                              "of number of NA values in each condition-level combination"
                          )
                      }

                      naFilterContrast <- rownames(dataMat) %in% rownames(dataMatNAFiltered)

                      if (leastRepCount == 0 && impute) { 
                          dataMatNAFiltered <- imputeGroupValues(
                              dataMatNAFiltered, 
                              conditionCombsOVR, 
                              minFraction=imputeMinFraction
                          )
                      }

                      if (type == "welch") {
                          statResults <- calculateWelch(
                              dataMatNAFiltered, 
                              groupHeader, 
                              c(groupLabel, restLabel)
                          )
                      }
	                      else if (type %in% c("limma", "limma_intensity")) {
	                          model <- setupModelFromDesign(designDfOVR, condCol, batchCol=batchCol, type=type)
	                          limmaDesignRaw <- stats::model.matrix(model)
	                          limmaPrepared <- sanitizeLimmaDesign(limmaDesignRaw)
	                          limmaDesign <- limmaPrepared$design
	                          limmaCoefMap <- limmaPrepared$coefMap
	                          limmaFit <- limma::lmFit(dataMatNAFiltered, limmaDesign)

	                          statResults <- calculateLimmaContrast(
	                              dataMatNAFiltered, 
	                              limmaDesign, 
	                              limmaFit, 
	                              c(groupLabel, restLabel), 
	                              useIntensityTrend = type == "limma_intensity",
	                              coefMap = limmaCoefMap
	                          )
	                      }
                      else {
                          stop("Unknown statistics type: ", type)
                      }

                      for (statMeasure in statMeasures) {
                          compLists[[statMeasure]][[compName]] <- c()
                          compLists[[statMeasure]][[compName]][naFilterContrast] <- statResults[[statMeasure]]
                      }
                  }
              }
              else {
                  if (is.null(comparisons)) {
                      stop("Argument 'comparisons' must be provided unless oneVsRest=TRUE.")
                  }

	                  comparisons <- as.character(comparisons)
	                  comparisons(nst) <- comparisons
	                  verifyContrasts(sampleReplicateGroupsStrings, comparisons, splitter=splitter)

                  if (!subsetByComparison) {

                  dataMatNAFiltered <- filterLowRep(
                      dataMat, 
                      conditionCombs, 
                      leastRep=leastRepCount
                  )

                  if (nrow(dataMatNAFiltered) == 0) {
                      stop("No rows remained after NA-filtering for condition: '", 
                           condCol, 
                           "' and batchCol: '", batchCol, "' (if empty then batchCol is not specified)\n",
                           "Consider whether you can reduce the 'leastRepCount' setting which sets the lower limit ",
                           "of number of NA values in each condition-level combination ",
                           "You could also try running without batchCol and see if there is enough data per condition then"
                           )
                  }

                  naFilterContrast <- rownames(dataMat) %in% rownames(dataMatNAFiltered)

                  if (leastRepCount == 0 && impute) { 
                      dataMatNAFiltered <- imputeGroupValues(
                          dataMatNAFiltered, 
                          conditionCombs, 
                          minFraction=imputeMinFraction
                      )
                  }

	                  model <- setupModel(nst, condCol, batchCol=batchCol, type=type)

	                  if (type %in% c("limma", "limma_intensity")) {
	                      limmaDesignRaw <- stats::model.matrix(model)
	                      limmaPrepared <- sanitizeLimmaDesign(limmaDesignRaw)
	                      limmaDesign <- limmaPrepared$design
	                      limmaCoefMap <- limmaPrepared$coefMap
	                      limmaFit <- limma::lmFit(dataMatNAFiltered, limmaDesign)
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
	                          c(level1, level2), 
	                          useIntensityTrend = FALSE,
	                          coefMap = limmaCoefMap)
	                  }
	                  else if (type == "limma_intensity") {
	                      
	                      statResults <- calculateLimmaContrast(
	                          dataMatNAFiltered, 
	                          limmaDesign, 
	                          limmaFit, 
	                          c(level1, level2), 
	                          useIntensityTrend = TRUE,
	                          coefMap = limmaCoefMap)
	                  }
                  else {
                      stop("Unknown statistics type: ", type)
                  }
                  
                  for (statMeasure in statMeasures) {
                      compLists[[statMeasure]][[comp]] <- c()
                      compLists[[statMeasure]][[comp]][naFilterContrast] <- statResults[[statMeasure]]
                  }

                  }
              }
              else {
                  for (comp in comparisons) {

                      compSplit <- unlist(strsplit(comp, splitter))

                      if (length(compSplit) != 2) {
                          stop("Comparison should be in format cond1-cond2 ", 
                               "here the split product was: ", 
                               paste(compSplit, collapse=" "))
                      }

                      level1 <- compSplit[1]
                      level2 <- compSplit[2]

                      groupMatch <- sampleReplicateGroupsStrings %in% c(level1, level2)
                      designDfComp <- designDf[groupMatch, , drop=FALSE]
                      dataMatComp <- dataMat[, groupMatch, drop=FALSE]

                      if (!is.null(batchCol)) {
                          conditionCombsComp <- paste(designDfComp[, condCol], designDfComp[, batchCol], sep="_")
                      }
                      else {
                          conditionCombsComp <- designDfComp[, condCol]
                      }

                      dataMatNAFiltered <- filterLowRep(
                          dataMatComp, 
                          conditionCombsComp, 
                          leastRep=leastRepCount
                      )

                      if (nrow(dataMatNAFiltered) == 0) {
                          stop(
                              "No rows remained after NA-filtering for comparison: '", 
                              comp,
                              "' (condCol: '", condCol,
                              "', batchCol: '", batchCol,
                              "' (if empty then batchCol is not specified))\n",
                              "Consider whether you can reduce the 'leastRepCount' setting which sets the lower limit ",
                              "of number of NA values in each condition-level combination"
                          )
                      }

                      naFilterContrast <- rownames(dataMat) %in% rownames(dataMatNAFiltered)

                      if (leastRepCount == 0 && impute) { 
                          dataMatNAFiltered <- imputeGroupValues(
                              dataMatNAFiltered, 
                              conditionCombsComp, 
                              minFraction=imputeMinFraction
                          )
                      }

                      if (type == "welch") {
                          statResults <- calculateWelch(
                              dataMatNAFiltered, 
                              as.character(designDfComp[, condCol]), 
                              c(level1, level2)
                          )
                      }
	                      else if (type %in% c("limma", "limma_intensity")) {
	                          model <- setupModelFromDesign(designDfComp, condCol, batchCol=batchCol, type=type)
	                          limmaDesignRaw <- stats::model.matrix(model)
	                          limmaPrepared <- sanitizeLimmaDesign(limmaDesignRaw)
	                          limmaDesign <- limmaPrepared$design
	                          limmaCoefMap <- limmaPrepared$coefMap
	                          limmaFit <- limma::lmFit(dataMatNAFiltered, limmaDesign)

	                          statResults <- calculateLimmaContrast(
	                              dataMatNAFiltered, 
	                              limmaDesign, 
	                              limmaFit, 
	                              c(level1, level2), 
	                              useIntensityTrend = type == "limma_intensity",
	                              coefMap = limmaCoefMap
	                          )
	                      }
                      else {
                          stop("Unknown statistics type: ", type)
                      }

                      for (statMeasure in statMeasures) {
                          compLists[[statMeasure]][[comp]] <- c()
                          compLists[[statMeasure]][[comp]][naFilterContrast] <- statResults[[statMeasure]]
                      }
                  }
              }
              }
              
              pairwiseCompsP(nst) <- compLists[["P"]]
              pairwiseCompsFdr(nst) <- compLists[["FDR"]]
              pairwiseCompsAve(nst) <- compLists[["Ave"]]
              pairwiseCompsFold(nst) <- compLists[["Fold"]]
              
              nst
          })


#' Check that a given contrast string is valid given a particular design
#' matrix. Each level tested for in the contrast should be present in the
#' condition column for the design matrix.
#' 
#' Mainly meant to verify strings received during server usage.
#'
#' @param designLevels Vector containing condition levels present in design
#' @param contrasts A string containing one or several (comma delimited)
#'   strings for which contrasts should be performed
#' @return None
#' @keywords internal
verifyContrasts <- function(designLevels, contrasts, splitter="-") {
    
    for (contrast in contrasts) {
        parts <- unlist(strsplit(contrast, splitter))
        
        if (length(parts) != 2) {
            stop("A contrast string delimited by one splitter (", splitter, ") was expected. Instead following was found: ", contrast)
        }
        
        if (!all(parts %in% designLevels)) {
            stop("There were issues in your contrast. \n", 
                 "All contrasts: ", paste(contrasts, collapse=", "), "\n",
                 "Part with issue: ", contrast, "\n", 
                 "Not all parts was found in the design column levels. Levels present in design: \n",
                 paste(unique(designLevels), collapse=", ")
                 )
        }
    }
}

setupModel <- function(nst, condCol, batchCol=NULL, type="limma") {
    
    if (is.null(batchCol)) {
        Variable <- as.factor(designDf(nst)[, condCol])
        model <- ~0+Variable
    }
    else {
        if (!(type %in% c("limma", "limma_intensity"))) {
            stop(
                "Batch compensation only compatible with Limma, got: ", 
                type
            )
        }
        Variable <- as.factor(designDf(nst)[, condCol])
        Batch <- as.factor(designDf(nst)[, batchCol])
        model <- ~0+Variable+Batch
    }
    model
}

sanitizeLimmaDesign <- function(limmaDesign) {

    rawNames <- colnames(limmaDesign)
    safeNames <- base::make.names(rawNames, unique=TRUE)
    colnames(limmaDesign) <- safeNames
    list(design=limmaDesign, coefMap=stats::setNames(safeNames, rawNames))
}

chooseOneVsRestLabel <- function(existingLabels,
                                 candidates=c("rest", "others", "all_other", "all_others")) {

    if (length(candidates) == 0) {
        stop("Expected at least one candidate label")
    }

    existing <- unique(as.character(existingLabels))

    for (candidate in candidates) {
        if (!(candidate %in% existing)) {
            return(candidate)
        }
    }

    labelBase <- candidates[1]
    suffix <- 1
    label <- paste0(labelBase, suffix)
    while (label %in% existing) {
        suffix <- suffix + 1
        label <- paste0(labelBase, suffix)
    }

    label
}

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
    statResults[["Ave"]] <- rowMeans(dataMat, na.rm=TRUE)
    statResults[["Fold"]] <- rowMeans(dataMat[, s1cols, drop=FALSE], na.rm=TRUE) - rowMeans(dataMat[, s2cols, drop=FALSE], na.rm=TRUE)

    statResults
}

calculateLimmaContrast <- function(dataMat, limmaDesign, limmaFit, levels, useIntensityTrend, coefMap=NULL) {

    coefLevel1 <- paste0("Variable", levels[1])
    coefLevel2 <- paste0("Variable", levels[2])

    if (!is.null(coefMap)) {
        if (!(coefLevel1 %in% names(coefMap))) {
            stop("Could not find limma coefficient name for level '", levels[1], "' (expected: '", coefLevel1, "')")
        }
        if (!(coefLevel2 %in% names(coefMap))) {
            stop("Could not find limma coefficient name for level '", levels[2], "' (expected: '", coefLevel2, "')")
        }

        coefLevel1 <- coefMap[[coefLevel1]]
        coefLevel2 <- coefMap[[coefLevel2]]
    }

    myContrast <- paste0(coefLevel1, "-", coefLevel2)
    contrastMatrix <- limma::makeContrasts(
        contrasts=c(myContrast), 
        levels=limmaDesign)
    fitContrasts <- limma::contrasts.fit(limmaFit, contrastMatrix)
    fitBayes <- limma::eBayes(fitContrasts, trend=useIntensityTrend)
    limmaTable <- limma::topTable(fitBayes, coef=1, number=Inf, sort.by="none")

    statResults <- list()
    statResults[["P"]] <- limmaTable$P.Value
    statResults[["FDR"]] <- limmaTable$adj.P.Val
    statResults[["Ave"]] <- limmaTable$AveExpr
    statResults[["Fold"]] <- limmaTable$logFC
    statResults
}
