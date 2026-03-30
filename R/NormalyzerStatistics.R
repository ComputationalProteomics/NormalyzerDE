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
#' @slot backendData Named list of backend-specific cached objects.
#' @slot filteringContrast Vector showing which entries are filtered
#'   (due to low count)
#' @slot pairwiseCompsP List with P-values for pairwise comparisons
#' @slot pairwiseCompsFdr List with FDR-values for pairwise comparisons
#' @slot pairwiseCompsAve List with average expression values
#' @slot pairwiseCompsFold List with log2 fold-change values for pairwise
#'   comparisons
#' @slot comparisons Character vector of the most recently used comparisons
#' @slot condCol Column containing last used conditions
#' @slot batchCol Column containing last used batch conditions
#' @slot splitter Delimiter separating contrast groups
NormalyzerStatistics <- setClass(
  "NormalyzerStatistics",
  slots = c(
    annotMat = "matrix",
    dataMat = "matrix",
    designDf = "data.frame",
    backendData = "list",

    pairwiseCompsP = "list",
    pairwiseCompsFdr = "list",
    pairwiseCompsAve = "list",
    pairwiseCompsFold = "list",
    pairwiseCompsSig = "list",

    comparisons = "character",
    condCol = "character",
    batchCol = "character",
    splitter = "character"
  )
)

#' Constructor for NormalyzerStatistics
#'
#' @param experimentObj Instance of SummarizedExperiment containing matrix
#'   and design information as column data
#' @param logTrans Whether the input data should be log2-transformed. When
#'   \code{TRUE}, non-finite values produced by the transform (e.g.
#'   \code{log2(0)} returning \code{-Inf}) are treated as missing (set to
#'   \code{NA}).
#' @return nds Generated NormalyzerStatistics instance
#' @export
#' @examples
#' data(example_stat_summarized_experiment)
#' nst <- NormalyzerStatistics(example_stat_summarized_experiment)
NormalyzerStatistics <- function(experimentObj, logTrans = FALSE) {
  dataMat <- as.matrix(SummarizedExperiment::assay(experimentObj))
  if (logTrans) {
    wasMissing <- is.na(dataMat)
    dataMat <- log2(dataMat)
    nonFinite <- !is.finite(dataMat) & !wasMissing
    if (any(nonFinite)) {
      cli::cli_warn(
        "Non-finite values produced by log2 transform (e.g. zeros or negative values) were treated as missing (set to NA).",
        class = "normalyzerde_warning",
        call = NULL
      )
      dataMat[nonFinite] <- NA_real_
    }
  }

  annotMat <- SummarizedExperiment::rowData(experimentObj)
  designDf <- SummarizedExperiment::colData(experimentObj)

  nst <- new(
    "NormalyzerStatistics",
    annotMat = as.matrix(annotMat),
    dataMat = dataMat,
    designDf = as.data.frame(designDf),
    backendData = list()
  )

  nst
}

setGeneric("condCol", function(object) {
  standardGeneric("condCol")
})
setMethod(
  "condCol",
  signature(object = "NormalyzerStatistics"),
  function(object) {
    slot(object, "condCol")
  }
)
setGeneric("condCol<-", function(object, value) {
  standardGeneric("condCol<-")
})
setReplaceMethod(
  "condCol",
  signature(object = "NormalyzerStatistics"),
  function(object, value) {
    slot(object, "condCol") <- value
    validObject(object)
    object
  }
)

setGeneric("batchCol", function(object) {
  standardGeneric("batchCol")
})
setMethod(
  "batchCol",
  signature(object = "NormalyzerStatistics"),
  function(object) {
    slot(object, "batchCol")
  }
)
setGeneric("batchCol<-", function(object, value) {
  standardGeneric("batchCol<-")
})
setReplaceMethod(
  "batchCol",
  signature(object = "NormalyzerStatistics"),
  function(object, value) {
    if (is.null(value) || length(value) == 0) {
      value <- character()
    } else {
      value <- as.character(value)
    }
    slot(object, "batchCol") <- value
    validObject(object)
    object
  }
)

setGeneric("comparisons", function(object) {
  standardGeneric("comparisons")
})
setMethod(
  "comparisons",
  signature(object = "NormalyzerStatistics"),
  function(object) {
    slot(object, "comparisons")
  }
)
setGeneric("comparisons<-", function(object, value) {
  standardGeneric("comparisons<-")
})
setReplaceMethod(
  "comparisons",
  signature(object = "NormalyzerStatistics"),
  function(object, value) {
    slot(object, "comparisons") <- value
    validObject(object)
    object
  }
)

setGeneric("contrastSplitter", function(object) {
  standardGeneric("contrastSplitter")
})
setMethod(
  "contrastSplitter",
  signature(object = "NormalyzerStatistics"),
  function(object) {
    slot(object, "splitter")
  }
)
setGeneric("contrastSplitter<-", function(object, value) {
  standardGeneric("contrastSplitter<-")
})
setReplaceMethod(
  "contrastSplitter",
  signature(object = "NormalyzerStatistics"),
  function(object, value) {
    slot(object, "splitter") <- as.character(value)
    validObject(object)
    object
  }
)

setGeneric("annotMat", function(object) {
  standardGeneric("annotMat")
})
setMethod(
  "annotMat",
  signature(object = "NormalyzerStatistics"),
  function(object) {
    slot(object, "annotMat")
  }
)

setGeneric("dataMat", function(object) {
  standardGeneric("dataMat")
})
setMethod(
  "dataMat",
  signature(object = "NormalyzerStatistics"),
  function(object) {
    slot(object, "dataMat")
  }
)
setGeneric("dataMat<-", function(object, value) {
  standardGeneric("dataMat<-")
})
setReplaceMethod(
  "dataMat",
  signature(object = "NormalyzerStatistics"),
  function(object, value) {
    slot(object, "dataMat") <- value
    validObject(object)
    object
  }
)

#' Backend-specific cached objects
#'
#' Some backends can optionally store intermediate objects for reuse or
#' debugging. For example, \code{calculateContrasts(..., type="limpa")} can keep
#' the completed expression matrix (as a limma \code{EList}), fitted
#' \code{MArrayLM} object(s), and estimated sample weights when available.
#'
#' @param object A \code{NormalyzerStatistics} object.
#' @return Named list of backend-specific objects.
#' @aliases backendData,NormalyzerStatistics-method
#' @aliases backendData<-,NormalyzerStatistics-method
#' @export
setGeneric("backendData", function(object) {
  standardGeneric("backendData")
})
setMethod(
  "backendData",
  signature(object = "NormalyzerStatistics"),
  function(object) {
    slot(object, "backendData")
  }
)

#' @rdname backendData
#' @param value Named list of backend-specific objects.
#' @export
setGeneric("backendData<-", function(object, value) {
  standardGeneric("backendData<-")
})
setReplaceMethod(
  "backendData",
  signature(object = "NormalyzerStatistics"),
  function(object, value) {
    slot(object, "backendData") <- value
    validObject(object)
    object
  }
)
setGeneric("filteredDataMat", function(object) {
  standardGeneric("filteredDataMat")
})
setMethod(
  "filteredDataMat",
  signature(object = "NormalyzerStatistics"),
  function(object) {
    slot(object, "filteredDataMat")
  }
)

setGeneric("designDf", function(object) {
  standardGeneric("designDf")
})
setMethod(
  "designDf",
  signature(object = "NormalyzerStatistics"),
  function(object) {
    slot(object, "designDf")
  }
)
setGeneric("designDf<-", function(object, value) {
  standardGeneric("designDf<-")
})
setReplaceMethod(
  "designDf",
  signature(object = "NormalyzerStatistics"),
  function(object, value) {
    slot(object, "designDf") <- value
    validObject(object)
    object
  }
)

setGeneric("filteringContrast", function(object) {
  standardGeneric("filteringContrast")
})
setMethod(
  "filteringContrast",
  signature(object = "NormalyzerStatistics"),
  function(object) {
    slot(object, "filteringContrast")
  }
)

setGeneric("pairwiseCompsP", function(object) {
  standardGeneric("pairwiseCompsP")
})
setMethod(
  "pairwiseCompsP",
  signature(object = "NormalyzerStatistics"),
  function(object) {
    slot(object, "pairwiseCompsP")
  }
)
setGeneric("pairwiseCompsP<-", function(object, value) {
  standardGeneric("pairwiseCompsP<-")
})
setReplaceMethod(
  "pairwiseCompsP",
  signature(object = "NormalyzerStatistics"),
  function(object, value) {
    slot(object, "pairwiseCompsP") <- value
    validObject(object)
    object
  }
)

setGeneric("pairwiseCompsFdr", function(object) {
  standardGeneric("pairwiseCompsFdr")
})
setMethod(
  "pairwiseCompsFdr",
  signature(object = "NormalyzerStatistics"),
  function(object) {
    slot(object, "pairwiseCompsFdr")
  }
)
setGeneric("pairwiseCompsFdr<-", function(object, value) {
  standardGeneric("pairwiseCompsFdr<-")
})
setReplaceMethod(
  "pairwiseCompsFdr",
  signature(object = "NormalyzerStatistics"),
  function(object, value) {
    slot(object, "pairwiseCompsFdr") <- value
    validObject(object)
    object
  }
)

setGeneric("pairwiseCompsAve", function(object) {
  standardGeneric("pairwiseCompsAve")
})
setMethod(
  "pairwiseCompsAve",
  signature(object = "NormalyzerStatistics"),
  function(object) {
    slot(object, "pairwiseCompsAve")
  }
)
setGeneric("pairwiseCompsAve<-", function(object, value) {
  standardGeneric("pairwiseCompsAve<-")
})
setReplaceMethod(
  "pairwiseCompsAve",
  signature(object = "NormalyzerStatistics"),
  function(object, value) {
    slot(object, "pairwiseCompsAve") <- value
    validObject(object)
    object
  }
)

setGeneric("pairwiseCompsFold", function(object) {
  standardGeneric("pairwiseCompsFold")
})
setMethod(
  "pairwiseCompsFold",
  signature(object = "NormalyzerStatistics"),
  function(object) {
    slot(object, "pairwiseCompsFold")
  }
)
setGeneric("pairwiseCompsFold<-", function(object, value) {
  standardGeneric("pairwiseCompsFold<-")
})
setReplaceMethod(
  "pairwiseCompsFold",
  signature(object = "NormalyzerStatistics"),
  function(object, value) {
    slot(object, "pairwiseCompsFold") <- value
    validObject(object)
    object
  }
)


#' Performs statistical comparisons between the supplied conditions.
#' It uses the design matrix and data matrix in the supplied
#' NormalyzerStatistics object. A column is supplied specifying which of the
#' columns in the design matrix that is used for deciding the sample groups.
#' The comparisons vector specifies which pairwise comparisons between
#' condition levels that are to be calculated.
#'
#' Optionally, a batch column can be specified allowing compensation for
#' covariate variation in the statistical model. This is only compatible
#' with a Limma- or limpa-based statistical analysis.
#'
#' @param nst Results evaluation object.
#' @param comparisons Character vector with pairwise comparisons for contrasts.
#'   Ignored if \code{oneVsRest=TRUE}.
#' @param condCol Column name in design matrix containing condition information.
#' @param batchCol Column name in design matrix containing batch information.
#' @param splitter Delimiter used to separate contrast groups.
#' @param type Type of statistical test ("limma", "limma_intensity", "welch" or
#'   "limpa"). "limpa" uses the Bioconductor package \pkg{limpa} to
#'   handle missing values via a detection probability curve (DPC) model.
#' @param leastRepCount Minimum number of replicates required in each group for
#'   contrast calculations. For \code{type="limpa"}, a feature is retained if at
#'   least one group has \code{leastRepCount} observed samples (and features
#'   entirely missing across all samples are removed).
#' @param impute Whether to impute values (ignored for \code{type="limpa"}).
#' @param imputeMinFraction Minimum fraction non-NA values for an analyte in any
#'   group to impute in other groups (ignored for \code{type="limpa"}).
#' @param subsetByComparison If TRUE, subset data and design to each comparison
#'   before NA-filtering, imputation and model fitting. Use this when filtering
#'   or imputation should depend only on the samples in each comparison, not the
#'   full dataset.
#' @param oneVsRest If TRUE, create one comparison per selected group against
#'   all remaining samples (for all groups in \code{condCol}, or the subset in
#'   \code{oneVsRestGroups}).
#' @param oneVsRestGroups Optional character vector specifying which groups in
#'   \code{condCol} to compare against all other samples.
#' @param limpaOptions Optional helper created by \code{\link{limpaOptions}}.
#'   For \code{type = "limpa"}, use this to configure protein-level
#'   summarization (\code{proteinIdCol}, \code{byRow}), DPC estimation
#'   (\code{dpc}, \code{dpcMethod}, \code{dpcArgs}), quantification
#'   (\code{quantArgs}, \code{quantifiedRds}), differential expression
#'   (\code{deArgs}), optional post-quantification normalization
#'   (\code{postQuantNorm}), and retained backend objects (\code{keep}).
#'   See \code{\link{limpaOptions}} for the available fields.
#' @return nst Statistics object with statistical measures calculated
#' @rdname calculateContrasts
#' @export
#' @examples
#' data(example_stat_summarized_experiment)
#' nst <- NormalyzerStatistics(example_stat_summarized_experiment)
#' results <- calculateContrasts(nst, c("1-2", "2-3"), "group")
#' resultsBatch <- calculateContrasts(nst, c("1-2", "2-3"), "group", batchCol="batch")
#' resultsOneVsRest <- calculateContrasts(nst, condCol="group", oneVsRest=TRUE)
setGeneric(
  name = "calculateContrasts",
  function(
    nst,
    comparisons = NULL,
    condCol,
    batchCol = NULL,
    splitter = "-",
    type = "limma",
    leastRepCount = 1,
    impute = FALSE,
    imputeMinFraction = 0.75,
    subsetByComparison = FALSE,
    oneVsRest = FALSE,
    oneVsRestGroups = NULL,
    limpaOptions = NULL
  ) {
    standardGeneric("calculateContrasts")
  }
)

#' @rdname calculateContrasts
setMethod(
  f = "calculateContrasts",
  signature = c("NormalyzerStatistics"),
  function(
    nst,
    comparisons = NULL,
    condCol,
    batchCol = NULL,
    splitter = "-",
    type = "limma",
    leastRepCount = 1,
    impute = FALSE,
    imputeMinFraction = 0.75,
    subsetByComparison = FALSE,
    oneVsRest = FALSE,
    oneVsRestGroups = NULL,
    limpaOptions = NULL
  ) {
    dataMat <- dataMat(nst)
    dataMatOriginalRowNames <- rownames(dataMat)
    designDf <- designDf(nst)

    contrastSplitter(nst) <- splitter
    condCol(nst) <- as.character(designDf[, condCol])

    if (
      !is.null(batchCol) && !(type %in% c("limma", "limma_intensity", "limpa"))
    ) {
      cli::cli_abort(
        "Batch compensation is only compatible with {.val limma} and {.val limpa}, got {.val {type}}.",
        class = "normalyzerde_error",
        call = NULL
      )
    }

    if (!is.null(batchCol)) {
      conditionCombs <- paste(
        designDf[, condCol],
        designDf[, batchCol],
        sep = "_"
      )
      batchCol(nst) <- designDf[, batchCol]
    } else {
      conditionCombs <- designDf[, condCol]
      batchCol(nst) <- character()
    }

    rownames(dataMat) <- seq_len(nrow(dataMat))

    sampleReplicateGroupsStrings <- as.character(designDf[, condCol])
    statMeasures <- c("P", "FDR", "Ave", "Fold")
    isLimpa <- identical(type, "limpa")

    warnIfLimpaInputLooksNotLog2 <- function(dataMat, threshold = 50) {
      finiteVals <- dataMat[is.finite(dataMat)]
      if (length(finiteVals) == 0) {
        return(invisible(NULL))
      }

      maxVal <- max(finiteVals)
      if (is.finite(maxVal) && maxVal > threshold) {
        maxValDisp <- signif(maxVal, 4)
        cli::cli_warn(
          c(
            "For {.arg type}={.val limpa}, the input data should be on the log2 scale (missing values as NA).",
            i = "These values are larger than typical log2 intensities (max finite value = {.val {maxValDisp}}).",
            i = "If your matrix is on the linear scale, set {.arg logTrans}={.val TRUE} in {.fn normalyzerDE} (or log2-transform upstream before calling {.fn calculateContrasts})."
          ),
          class = "normalyzerde_warning",
          call = NULL
        )
      }

      invisible(NULL)
    }

    limpaOptions <- if (isLimpa) {
      resolveLimpaOptions(limpaOptions)
    } else {
      limpaOptions
    }
    limpaProteinIdCol <- if (isLimpa) limpaOptions$proteinIdCol else NULL
    limpaByRow <- if (isLimpa) limpaOptions$byRow else FALSE
    limpaDpc <- if (isLimpa) limpaOptions$dpc else NULL
    limpaDpcMethod <- if (isLimpa) limpaOptions$dpcMethod else "none"
    limpaDpcArgs <- if (isLimpa) limpaOptions$dpcArgs else list()
    limpaQuantArgs <- if (isLimpa) limpaOptions$quantArgs else list()
    limpaQuantifiedRds <- if (isLimpa) limpaOptions$quantifiedRds else NULL
    limpaDEArgs <- if (isLimpa) limpaOptions$deArgs else list()
    limpaKeep <- if (isLimpa) limpaOptions$keep else "none"
    limpaPostQuantNorm <- if (isLimpa) limpaOptions$postQuantNorm else "none"

    limpaPostQuantNormUse <- "none"
    if (isLimpa) {
      requireLimpaPackageInternal("Statistics type 'limpa'")
      limpaQuantByRow <- resolveLimpaQuantByRowFn()
      warnIfLimpaInputLooksNotLog2(dataMat)

      limpaDpcMethod <- match.arg(
        limpaDpcMethod,
        c("none", "dpc", "dpcON", "dpcCN")
      )
      limpaKeep <- match.arg(limpaKeep, c("none", "elist", "fit", "all"))
      limpaPostQuantNormUse <- match.arg(
        limpaPostQuantNorm,
        c(
          "none",
          "GI",
          "median",
          "mean",
          "Quantile",
          "CycLoess",
          "RLR",
          "quantile"
        )
      )

      byRowConfig <- normalizeLimpaByRowConfig(
        limpaByRow = limpaByRow,
        limpaProteinIdCol = limpaProteinIdCol
      )
      limpaByRow <- byRowConfig$limpaByRow
      limpaProteinIdCol <- byRowConfig$limpaProteinIdCol

      if (!is.null(limpaQuantifiedRds)) {
        if (
          !is.character(limpaQuantifiedRds) || length(limpaQuantifiedRds) != 1
        ) {
          cli::cli_abort(
            "{.arg quantifiedRds} must be a length-1 character file path (or {.val NULL}).",
            class = "normalyzerde_error",
            call = NULL
          )
        }
        limpaQuantifiedRds <- as.character(limpaQuantifiedRds)[1]
        if (is.na(limpaQuantifiedRds) || !nzchar(limpaQuantifiedRds)) {
          cli::cli_abort(
            "{.arg quantifiedRds} must be a non-empty file path (or {.val NULL}).",
            class = "normalyzerde_error",
            call = NULL
          )
        }
        if (!file.exists(limpaQuantifiedRds)) {
          cli::cli_abort(
            "{.arg quantifiedRds} file does not exist: {.path {limpaQuantifiedRds}}.",
            class = "normalyzerde_error",
            call = NULL
          )
        }
      }

      validateLimpaDpc(limpaDpc)
    }

    if (!isLimpa && !is.null(limpaQuantifiedRds)) {
      cli::cli_abort(
        "{.arg quantifiedRds} is only supported for {.arg type}={.val limpa}.",
        class = "normalyzerde_error",
        call = NULL
      )
    }

    limpaDpcArgsUse <- if (isLimpa) {
      sanitizeLimpaDpcArgs(limpaDpcArgs)
    } else {
      list()
    }

    limpaQuantArgsUse <- if (isLimpa) {
      applyLimpaQuantDefaultsAndValidate(sanitizeLimpaQuantArgs(limpaQuantArgs))
    } else {
      list()
    }
    limpaDEArgsUse <- if (isLimpa) {
      applyLimpaDEDefaultsAndValidate(sanitizeLimpaDEArgs(limpaDEArgs))
    } else {
      list()
    }
    limpaStoreSampleWeights <- if (isLimpa) {
      isTRUE(limpaDEArgsUse[["sample.weights"]])
    } else {
      FALSE
    }

    limpaVerboseUse <- if (isLimpa) {
      isTRUE(limpaQuantArgsUse[["verbose"]])
    } else {
      FALSE
    }
    limpaDpcSlopeUse <- if (isLimpa) {
      as.numeric(limpaQuantArgsUse[["dpc.slope"]])
    } else {
      0.8
    }

    if (isLimpa && leastRepCount > 1 && isTRUE(limpaVerboseUse)) {
      cli::cli_inform(
        c(
          "!" = "For {.arg type}={.val limpa}, {.arg leastRepCount}>1 can drop informative sparse features.",
          i = "Consider {.arg leastRepCount}=1 unless you explicitly want to filter sparse rows."
        )
      )
    }

    limpaDpcUse <- if (isLimpa) {
      if (!is.null(limpaDpc)) {
        if (limpaDpcMethod != "none" && isTRUE(limpaVerboseUse)) {
          cli::cli_inform(
            c(
              i = "{.arg dpc} was supplied; ignoring {.arg dpcMethod}={.val {limpaDpcMethod}}."
            )
          )
        }
        limpaDpc
      } else {
        estimateLimpaDpcFromData(
          dataMat = dataMat,
          limpaDpcMethod = limpaDpcMethod,
          limpaDpcArgs = limpaDpcArgsUse,
          dpcSlope = limpaDpcSlopeUse,
          verbose = limpaVerboseUse
        )
      }
    } else {
      NULL
    }

    limpaBackend <- NULL
    recordLimpaBackend <- function(key, y = NULL, fit = NULL, design = NULL) {
      invisible(NULL)
    }

    if (isLimpa && (limpaKeep != "none" || limpaStoreSampleWeights)) {
      dpcVec <- if (is.null(limpaDpcUse)) {
        NULL
      } else if (is.list(limpaDpcUse)) {
        limpaDpcUse$dpc
      } else {
        as.numeric(limpaDpcUse)
      }

      dpcMethodUsed <- if (!is.null(limpaDpc)) {
        "supplied"
      } else {
        limpaDpcMethod
      }

      limpaBackend <- list(
        dpc = dpcVec,
        dpcMethod = dpcMethodUsed,
        dpcSlopeInput = as.numeric(limpaDpcSlopeUse),
        dpcSlopeUsed = if (is.null(dpcVec)) {
          as.numeric(limpaDpcSlopeUse)
        } else {
          as.numeric(dpcVec[[2]])
        },
        postQuantNorm = limpaPostQuantNormUse,
        fits = list(),
        designs = list(),
        elists = list(),
        sampleWeights = list()
      )

      if (length(limpaDpcArgsUse) > 0) {
        limpaBackend$dpcArgs <- limpaDpcArgsUse
      }
      if (length(limpaQuantArgsUse) > 0) {
        limpaBackend$quantArgs <- limpaQuantArgsUse
      }

      recordLimpaBackend <- function(key, y = NULL, fit = NULL, design = NULL) {
        key <- as.character(key)[1]
        if (is.na(key) || !nzchar(key)) {
          key <- ".global"
        }

        if (limpaKeep %in% c("fit", "all") && !is.null(fit)) {
          limpaBackend$fits[[key]] <<- fit
          limpaBackend$designs[[key]] <<- design
        }
        if (limpaKeep %in% c("elist", "all") && !is.null(y)) {
          limpaBackend$elists[[key]] <<- y
        }
        if (limpaStoreSampleWeights && !is.null(fit)) {
          sampleWeights <- extractLimpaSampleWeightsFromFit(fit)
          if (!is.null(sampleWeights)) {
            limpaBackend$sampleWeights[[key]] <<- sampleWeights
          }
        }

        invisible(NULL)
      }
    }

    prepareLimpaQuantified <- function(
      dataMat,
      annotationMat,
      proteinIdCol
    ) {
      if (is.null(proteinIdCol)) {
        return(NULL)
      }

      proteinId <- annotationMat[, proteinIdCol]
      proteinId <- as.character(proteinId)

      if (anyDuplicated(proteinId) == 0) {
        return(NULL)
      }

      keepRows <- rowSums(!is.na(dataMat)) > 0
      dataMat <- dataMat[keepRows, , drop = FALSE]
      proteinId <- proteinId[keepRows]

      genesInput <- as.data.frame(
        annotationMat[keepRows, , drop = FALSE],
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      genesInput[[proteinIdCol]] <- proteinId

      quantifyLimpaByProteinInternal(
        dataMat = dataMat,
        genesDf = genesInput,
        proteinIdCol = proteinIdCol,
        dpc = limpaDpcUse,
        quantArgs = limpaQuantArgsUse
      )
    }

    filterLowRepLimpa <- function(
      df,
      groups,
      leastRep = 1,
      nObservations = NULL
    ) {
      obsMat <- if (!is.null(nObservations)) {
        if (!is.matrix(nObservations)) {
          cli::cli_abort(
            "{.arg nObservations} must be a matrix when provided.",
            class = "normalyzerde_error",
            call = NULL
          )
        }
        if (!all(dim(nObservations) == dim(df))) {
          cli::cli_abort(
            "{.arg nObservations} must have the same dimensions as {.arg df}.",
            class = "normalyzerde_error",
            call = NULL
          )
        }
        nObservations > 0
      } else {
        !is.na(df)
      }

      hasAnyObs <- rowSums(obsMat) > 0
      df <- df[hasAnyObs, , drop = FALSE]
      obsMat <- obsMat[hasAnyObs, , drop = FALSE]
      if (nrow(df) == 0) {
        return(df)
      }

      leastRep <- as.integer(leastRep)
      if (leastRep <= 1) {
        return(df)
      }

      groups <- as.character(groups)
      groupLevels <- unique(groups)

      maxCount <- integer(nrow(df))
      for (groupLevel in groupLevels) {
        cols <- which(groups == groupLevel)
        if (length(cols) == 0) {
          next
        }
        counts <- rowSums(obsMat[, cols, drop = FALSE])
        maxCount <- pmax(maxCount, counts)
      }

      df[maxCount >= leastRep, , drop = FALSE]
    }

    warnIfLimpaProteinSummarizationWillBeSkipped <- function(
      annotationMat,
      proteinIdCol
    ) {
      if (is.null(proteinIdCol)) {
        return(invisible(NULL))
      }

      proteinId <- annotationMat[, proteinIdCol]
      proteinId <- as.character(proteinId)

      if (anyDuplicated(proteinId) == 0) {
        cli::cli_warn(
          c(
            "{.arg proteinIdCol}={.val {proteinIdCol}} contains no duplicated identifiers, so no peptide/precursor-to-protein summarization will be performed.",
            i = "The analysis will proceed with input rows as features (this may indicate the input is already protein-level)."
          ),
          class = "normalyzerde_warning",
          call = NULL
        )
      }

      invisible(NULL)
    }

    resolveLimpaCachedRowMap <- function(
      dataRowIds,
      originalRowIds,
      cachedRowIds
    ) {
      dataRowIds <- as.character(dataRowIds)
      cachedRowIds <- as.character(cachedRowIds)

      if (all(dataRowIds %in% cachedRowIds)) {
        return(stats::setNames(dataRowIds, dataRowIds))
      }

      if (!is.null(originalRowIds)) {
        originalRowIds <- as.character(originalRowIds)
        if (
          length(originalRowIds) == length(dataRowIds) &&
            anyDuplicated(originalRowIds) == 0 &&
            all(originalRowIds %in% cachedRowIds)
        ) {
          return(stats::setNames(originalRowIds, dataRowIds))
        }
      }

      cli::cli_abort(
        c(
          "{.arg quantifiedRds} row identifiers could not be matched to the data matrix.",
          i = "Row-order fallback is not allowed because it can silently misalign cached missingness metadata.",
          i = "Ensure the quantified {.code EList} was generated from the same matrix (rows), or provide matching row names."
        ),
        class = "normalyzerde_error",
        call = NULL
      )
    }

    limpaCachedRowMap <- NULL
    limpaQuantifiedCached <- NULL
    if (type == "limpa" && !is.null(limpaQuantifiedRds)) {
      limpaQuantifiedCached <- readRDS(limpaQuantifiedRds)
      if (!inherits(limpaQuantifiedCached, "EList")) {
        cli::cli_abort(
          c(
            "{.arg quantifiedRds} must contain a limma {.code EList} object.",
            i = "Got class(es): {.val {paste(class(limpaQuantifiedCached), collapse = \", \")}}."
          ),
          class = "normalyzerde_error",
          call = NULL
        )
      }
      if (
        is.null(limpaQuantifiedCached$E) || !is.matrix(limpaQuantifiedCached$E)
      ) {
        cli::cli_abort(
          "{.arg quantifiedRds} {.code EList} must contain a matrix element {.code E}.",
          class = "normalyzerde_error",
          call = NULL
        )
      }
      if (is.null(colnames(limpaQuantifiedCached$E))) {
        cli::cli_abort(
          "{.arg quantifiedRds} {.code EList$E} must contain sample column names.",
          class = "normalyzerde_error",
          call = NULL
        )
      }

      cachedRowIds <- rownames(limpaQuantifiedCached$E)
      if (is.null(cachedRowIds)) {
        cachedRowIds <- as.character(seq_len(nrow(limpaQuantifiedCached$E)))
      }
      cachedRowIds <- as.character(cachedRowIds)
      if (
        anyNA(cachedRowIds) ||
          any(!nzchar(cachedRowIds)) ||
          anyDuplicated(cachedRowIds) > 0
      ) {
        cli::cli_abort(
          "{.arg quantifiedRds} {.code EList$E} must have unique, non-empty row names.",
          class = "normalyzerde_error",
          call = NULL
        )
      }
      rownames(limpaQuantifiedCached$E) <- cachedRowIds

      if (
        is.null(limpaQuantifiedCached$other$n.observations) ||
          !is.matrix(limpaQuantifiedCached$other$n.observations)
      ) {
        cli::cli_abort(
          c(
            "{.arg quantifiedRds} {.code EList} must contain {.code other$n.observations} (a matrix).",
            i = "Use an {.code EList} returned by {.code limpa::dpcQuant()} or {.code limpa::dpcQuantByRow()}."
          ),
          class = "normalyzerde_error",
          call = NULL
        )
      }
      if (
        !all(
          dim(limpaQuantifiedCached$other$n.observations) ==
            dim(limpaQuantifiedCached$E)
        )
      ) {
        cli::cli_abort(
          "{.arg quantifiedRds} {.code other$n.observations} must have the same dimensions as {.code EList$E}.",
          class = "normalyzerde_error",
          call = NULL
        )
      }
      rownames(limpaQuantifiedCached$other$n.observations) <- cachedRowIds
      if (is.null(colnames(limpaQuantifiedCached$other$n.observations))) {
        colnames(limpaQuantifiedCached$other$n.observations) <- colnames(
          limpaQuantifiedCached$E
        )
      }

      if (
        is.null(limpaQuantifiedCached$other$standard.error) ||
          !is.matrix(limpaQuantifiedCached$other$standard.error)
      ) {
        cli::cli_abort(
          c(
            "{.arg quantifiedRds} {.code EList} must contain {.code other$standard.error} (a matrix).",
            i = "Use an {.code EList} returned by {.code limpa::dpcQuant()} or {.code limpa::dpcQuantByRow()}."
          ),
          class = "normalyzerde_error",
          call = NULL
        )
      }
      if (
        !all(
          dim(limpaQuantifiedCached$other$standard.error) ==
            dim(limpaQuantifiedCached$E)
        )
      ) {
        cli::cli_abort(
          "{.arg quantifiedRds} {.code other$standard.error} must have the same dimensions as {.code EList$E}.",
          class = "normalyzerde_error",
          call = NULL
        )
      }
      rownames(limpaQuantifiedCached$other$standard.error) <- cachedRowIds
      if (is.null(colnames(limpaQuantifiedCached$other$standard.error))) {
        colnames(limpaQuantifiedCached$other$standard.error) <- colnames(
          limpaQuantifiedCached$E
        )
      }

      if (anyNA(dataMat)) {
        cli::cli_abort(
          c(
            "{.arg quantifiedRds} requires a completed expression matrix without NA values.",
            i = "Use the post-quant matrix (for example the {.val log2} output from {.code normalyzer(preQuant = 'limpa')})."
          ),
          class = "normalyzerde_error",
          call = NULL
        )
      }
    }

    limpaProteinIdColUsed <- if (
      type == "limpa" && is.null(limpaQuantifiedCached)
    ) {
      inferLimpaProteinIdCol(annotMat(nst), limpaProteinIdCol)
    } else {
      NULL
    }

    if (type == "limpa" && is.null(limpaQuantifiedCached)) {
      warnIfLimpaProteinSummarizationWillBeSkipped(
        annotMat(nst),
        limpaProteinIdColUsed
      )
    }

    limpaQuantified <- if (type != "limpa") {
      NULL
    } else if (!is.null(limpaQuantifiedCached)) {
      limpaQuantifiedCached
    } else if (!is.null(limpaProteinIdColUsed)) {
      prepareLimpaQuantified(dataMat, annotMat(nst), limpaProteinIdColUsed)
    } else {
      NULL
    }

    if (type == "limpa" && !is.null(limpaQuantifiedCached)) {
      missingCols <- setdiff(colnames(dataMat), colnames(limpaQuantified$E))
      if (length(missingCols) > 0) {
        cli::cli_abort(
          "{.arg quantifiedRds} is missing sample columns required by the data matrix: {.val {paste(missingCols, collapse = \", \")}}.",
          class = "normalyzerde_error",
          call = NULL
        )
      }

      limpaCachedRowMap <- resolveLimpaCachedRowMap(
        dataRowIds = rownames(dataMat),
        originalRowIds = dataMatOriginalRowNames,
        cachedRowIds = rownames(limpaQuantified$E)
      )
    }

    if (!is.null(limpaBackend)) {
      limpaBackend$proteinIdCol <- limpaProteinIdColUsed
      limpaBackend$usedDpcQuant <- !is.null(limpaQuantified) &&
        is.null(limpaQuantifiedCached)
      if (!is.null(limpaQuantifiedCached)) {
        limpaBackend$quantifiedRds <- limpaQuantifiedRds
      }
      if (!is.null(limpaQuantified) && limpaKeep %in% c("elist", "all")) {
        limpaBackend$quantifiedEList <- limpaQuantified
      }
    }

    if (!is.null(limpaQuantified) && is.null(limpaQuantifiedCached)) {
      dataMat <- as.matrix(limpaQuantified$E)
      slot(nst, "dataMat") <- dataMat
      slot(nst, "annotMat") <- as.matrix(limpaQuantified$genes)
    }

    applyLimpaPostQuantNorm <- function(y) {
      if (identical(limpaPostQuantNormUse, "none")) {
        return(y)
      }

      if (!is.null(y$E) && anyNA(y$E)) {
        cli::cli_abort(
          "{.arg postQuantNorm}={.val {limpaPostQuantNormUse}} requires a completed expression matrix without NA values.",
          class = "normalyzerde_error",
          call = NULL
        )
      }

      if (identical(limpaPostQuantNormUse, "GI")) {
        y$E <- globalIntensityNormalization(y$E, noLogTransform = TRUE)
      } else if (identical(limpaPostQuantNormUse, "median")) {
        y$E <- medianNormalization(y$E, noLogTransform = TRUE)
      } else if (identical(limpaPostQuantNormUse, "mean")) {
        y$E <- meanNormalization(y$E, noLogTransform = TRUE)
      } else if (limpaPostQuantNormUse %in% c("Quantile", "quantile")) {
        y$E <- performQuantileNormalization(y$E, noLogTransform = TRUE)
      } else if (identical(limpaPostQuantNormUse, "CycLoess")) {
        y$E <- performCyclicLoessNormalization(y$E, noLogTransform = TRUE)
      } else if (identical(limpaPostQuantNormUse, "RLR")) {
        y$E <- performGlobalRLRNormalization(y$E, noLogTransform = TRUE)
      } else {
        cli::cli_abort(
          "Unknown {.arg postQuantNorm} value: {.val {limpaPostQuantNormUse}}.",
          class = "normalyzerde_error",
          call = NULL
        )
      }

      y
    }

    calculateLimpaFit <- function(
      dataMat,
      limmaDesign,
      backendKey = ".global"
    ) {
      if (is.null(limpaQuantified)) {
        yImputed <- do.call(
          limpaQuantByRow,
          c(
            list(
              y = dataMat,
              dpc = limpaDpcUse
            ),
            limpaQuantArgsUse
          )
        )
        yImputed <- applyLimpaPostQuantNorm(yImputed)

        fit <- do.call(
          limpa::dpcDE,
          c(
            list(y = yImputed, design = limmaDesign, plot = FALSE),
            limpaDEArgsUse
          )
        )
        recordLimpaBackend(
          backendKey,
          y = yImputed,
          fit = fit,
          design = limmaDesign
        )
        return(fit)
      }

      yUse <- limpaQuantified

      cols <- colnames(dataMat)
      yUse$E <- yUse$E[, cols, drop = FALSE]

      if (!is.null(yUse$other$n.observations)) {
        yUse$other$n.observations <- yUse$other$n.observations[,
          cols,
          drop = FALSE
        ]
      }
      if (!is.null(yUse$other$standard.error)) {
        yUse$other$standard.error <- yUse$other$standard.error[,
          cols,
          drop = FALSE
        ]
      }

      rows <- rownames(dataMat)
      rowsInQuant <- if (!is.null(limpaCachedRowMap)) {
        unname(limpaCachedRowMap[rows])
      } else {
        rows
      }
      if (anyNA(rowsInQuant)) {
        cli::cli_abort(
          c(
            "Failed to map data matrix rows to {.arg quantifiedRds} rows.",
            i = "Please ensure both inputs come from the same quantified matrix."
          ),
          class = "normalyzerde_error",
          call = NULL
        )
      }

      yUse$E <- yUse$E[rowsInQuant, , drop = FALSE]
      rownames(yUse$E) <- rows

      if (!is.null(yUse$other$n.observations)) {
        yUse$other$n.observations <- yUse$other$n.observations[
          rowsInQuant,
          ,
          drop = FALSE
        ]
        rownames(yUse$other$n.observations) <- rows
      }
      if (!is.null(yUse$other$standard.error)) {
        yUse$other$standard.error <- yUse$other$standard.error[
          rowsInQuant,
          ,
          drop = FALSE
        ]
        rownames(yUse$other$standard.error) <- rows
      }
      if (!is.null(yUse$genes)) {
        genesDf <- as.data.frame(yUse$genes, check.names = FALSE)
        geneRows <- rownames(genesDf)
        if (!is.null(geneRows) && all(rowsInQuant %in% geneRows)) {
          genesDf <- genesDf[rowsInQuant, , drop = FALSE]
        } else if (nrow(genesDf) == nrow(limpaQuantified$E)) {
          genesDf <- genesDf[
            match(rowsInQuant, rownames(limpaQuantified$E)),
            ,
            drop = FALSE
          ]
        } else {
          cli::cli_abort(
            "{.arg quantifiedRds} genes rows could not be aligned to {.code EList} rows.",
            class = "normalyzerde_error",
            call = NULL
          )
        }
        rownames(genesDf) <- rows
        yUse$genes <- genesDf
      }

      if (!is.null(limpaQuantifiedCached)) {
        if (anyNA(dataMat)) {
          cli::cli_abort(
            c(
              "{.arg quantifiedRds} requires a completed expression matrix without NA values.",
              i = "Use the post-quant matrix (for example the {.val log2} output from {.code normalyzer(preQuant = 'limpa')})."
            ),
            class = "normalyzerde_error",
            call = NULL
          )
        }
        yUse$E <- dataMat
      }

      yUse <- applyLimpaPostQuantNorm(yUse)

      fit <- do.call(
        limpa::dpcDE,
        c(list(y = yUse, design = limmaDesign, plot = FALSE), limpaDEArgsUse)
      )
      recordLimpaBackend(backendKey, y = yUse, fit = fit, design = limmaDesign)
      fit
    }

    subsetLimpaNObservations <- function(dataMatSubset) {
      if (
        is.null(limpaQuantified) ||
          is.null(limpaQuantified$other$n.observations)
      ) {
        return(NULL)
      }

      rowIds <- rownames(dataMatSubset)
      if (!is.null(limpaCachedRowMap)) {
        rowIds <- unname(limpaCachedRowMap[rowIds])
      }
      if (anyNA(rowIds)) {
        cli::cli_abort(
          "Failed to map filtered data rows to {.arg quantifiedRds} {.code n.observations} rows.",
          class = "normalyzerde_error",
          call = NULL
        )
      }

      obs <- limpaQuantified$other$n.observations[
        rowIds,
        colnames(dataMatSubset),
        drop = FALSE
      ]
      rownames(obs) <- rownames(dataMatSubset)
      obs
    }

    filterContrastData <- function(
      dataMatCurrent,
      conditionHeaders,
      nObservations = NULL
    ) {
      dataMatNAFiltered <- if (type == "limpa") {
        filterLowRepLimpa(
          dataMatCurrent,
          conditionHeaders,
          leastRep = leastRepCount,
          nObservations = nObservations
        )
      } else {
        filterLowRep(
          dataMatCurrent,
          conditionHeaders,
          leastRep = leastRepCount
        )
      }

      if (type != "limpa" && leastRepCount == 0 && impute) {
        dataMatNAFiltered <- imputeGroupValues(
          dataMatNAFiltered,
          conditionHeaders,
          minFraction = imputeMinFraction
        )
      }

      dataMatNAFiltered
    }

    prepareLimmaContrastState <- function(
      dataMatFiltered,
      designDfCurrent,
      backendKey = ".global"
    ) {
      model <- setupModelFromDesign(
        designDfCurrent,
        condCol,
        batchCol = batchCol,
        type = type
      )
      limmaDesignRaw <- stats::model.matrix(model)
      limmaPrepared <- sanitizeLimmaDesign(limmaDesignRaw)
      limmaFit <- if (type == "limpa") {
        calculateLimpaFit(
          dataMatFiltered,
          limmaPrepared$design,
          backendKey = backendKey
        )
      } else {
        limma::lmFit(dataMatFiltered, limmaPrepared$design)
      }

      list(
        design = limmaPrepared$design,
        fit = limmaFit,
        coefMap = limmaPrepared$coefMap
      )
    }

    calculateContrastStatistics <- function(
      dataMatFiltered,
      groupHeader,
      levels,
      designDfCurrent = NULL,
      backendKey = ".global",
      limmaState = NULL
    ) {
      if (type == "welch") {
        return(calculateWelch(dataMatFiltered, groupHeader, levels))
      }

      if (type %in% c("limma", "limma_intensity", "limpa")) {
        if (is.null(limmaState)) {
          limmaState <- prepareLimmaContrastState(
            dataMatFiltered,
            designDfCurrent,
            backendKey = backendKey
          )
        }

        return(calculateLimmaContrast(
          dataMatFiltered,
          limmaState$design,
          limmaState$fit,
          levels,
          useIntensityTrend = type == "limma_intensity",
          coefMap = limmaState$coefMap
        ))
      }

      cli::cli_abort(
        "Unknown statistics {.arg type}: {.val {type}}.",
        class = "normalyzerde_error",
        call = NULL
      )
    }

    compLists <- initializeContrastResultLists(statMeasures)

    if (oneVsRest) {
      restLabel <- chooseOneVsRestLabel(sampleReplicateGroupsStrings)

      targetGroups <- if (is.null(oneVsRestGroups)) {
        unique(sampleReplicateGroupsStrings)
      } else {
        unique(as.character(oneVsRestGroups))
      }

      if (length(targetGroups) == 0) {
        cli::cli_abort(
          "No groups specified for one-vs-rest comparisons.",
          class = "normalyzerde_error",
          call = NULL
        )
      }

      missingGroups <- setdiff(
        targetGroups,
        unique(sampleReplicateGroupsStrings)
      )
      if (length(missingGroups) > 0) {
        cli::cli_abort(
          c(
            "Some groups in {.arg oneVsRestGroups} were not found in {.arg condCol}={.val {condCol}}.",
            i = "Missing groups: {.val {paste(missingGroups, collapse = \", \")}}."
          ),
          class = "normalyzerde_error",
          call = NULL
        )
      }

      comparisonsGenerated <- paste(targetGroups, restLabel, sep = splitter)
      comparisons(nst) <- comparisonsGenerated

      for (groupLabel in targetGroups) {
        compName <- paste(groupLabel, restLabel, sep = splitter)
        groupHeader <- ifelse(
          sampleReplicateGroupsStrings %in% groupLabel,
          groupLabel,
          restLabel
        )

        designDfOVR <- designDf
        designDfOVR[, condCol] <- groupHeader

        conditionCombsOVR <- if (!is.null(batchCol)) {
          paste(groupHeader, designDfOVR[, batchCol], sep = "_")
        } else {
          groupHeader
        }

        dataMatNAFiltered <- filterContrastData(
          dataMat,
          conditionCombsOVR,
          nObservations = subsetLimpaNObservations(dataMat)
        )

        if (nrow(dataMatNAFiltered) == 0) {
          cli::cli_abort(
            c(
              "No rows remained after NA-filtering for one-vs-rest comparison {.val {compName}}.",
              i = "{.arg condCol}: {.val {condCol}}",
              i = "{.arg batchCol}: {.val {batchCol}}",
              i = "Consider reducing {.arg leastRepCount} (lower limit for non-missing values within each condition-level combination)."
            ),
            class = "normalyzerde_error",
            call = NULL
          )
        }

        naFilterContrast <- rownames(dataMat) %in% rownames(dataMatNAFiltered)

        statResults <- calculateContrastStatistics(
          dataMatNAFiltered,
          groupHeader,
          c(groupLabel, restLabel),
          designDfCurrent = designDfOVR,
          backendKey = compName
        )

        compLists <- storeContrastResults(
          compLists,
          statResults,
          compName,
          naFilterContrast
        )
      }
    } else {
      if (is.null(comparisons)) {
        cli::cli_abort(
            "Provide {.arg comparisons}, or set {.arg oneVsRest}={.val TRUE}.",
          class = "normalyzerde_error",
          call = NULL
        )
      }

      comparisons <- as.character(comparisons)
      comparisons(nst) <- comparisons
      verifyContrasts(
        sampleReplicateGroupsStrings,
        comparisons,
        splitter = splitter
      )

      if (!subsetByComparison) {
        dataMatNAFiltered <- filterContrastData(
          dataMat,
          conditionCombs,
          nObservations = subsetLimpaNObservations(dataMat)
        )

        if (nrow(dataMatNAFiltered) == 0) {
          cli::cli_abort(
            c(
              "No rows remained after NA-filtering for {.arg condCol}={.val {condCol}}.",
              i = "{.arg batchCol}: {.val {batchCol}}",
              i = "Consider reducing {.arg leastRepCount} (lower limit for non-missing values within each condition-level combination).",
              i = "You could also try running without {.arg batchCol} and check if there is enough data per condition."
            ),
            class = "normalyzerde_error",
            call = NULL
          )
        }

        naFilterContrast <- rownames(dataMat) %in% rownames(dataMatNAFiltered)

        limmaState <- if (type %in% c("limma", "limma_intensity", "limpa")) {
          prepareLimmaContrastState(
            dataMatNAFiltered,
            designDf,
            backendKey = ".global"
          )
        } else {
          NULL
        }

        for (comp in comparisons) {
          compLevels <- parseContrastLevels(comp, splitter)
          assertContrastLevelsPresent(compLevels, sampleReplicateGroupsStrings)

          statResults <- calculateContrastStatistics(
            dataMatNAFiltered,
            sampleReplicateGroupsStrings,
            compLevels,
            limmaState = limmaState
          )

          compLists <- storeContrastResults(
            compLists,
            statResults,
            comp,
            naFilterContrast
          )
        }
      } else {
        for (comp in comparisons) {
          compLevels <- parseContrastLevels(comp, splitter)
          assertContrastLevelsPresent(compLevels, sampleReplicateGroupsStrings)
          level1 <- compLevels[1]
          level2 <- compLevels[2]

          groupMatch <- sampleReplicateGroupsStrings %in% c(level1, level2)
          designDfComp <- designDf[groupMatch, , drop = FALSE]
          dataMatComp <- dataMat[, groupMatch, drop = FALSE]

          if (!is.null(batchCol)) {
            conditionCombsComp <- paste(
              designDfComp[, condCol],
              designDfComp[, batchCol],
              sep = "_"
            )
          } else {
            conditionCombsComp <- designDfComp[, condCol]
          }

          dataMatNAFiltered <- filterContrastData(
            dataMatComp,
            conditionCombsComp,
            nObservations = subsetLimpaNObservations(dataMatComp)
          )

          if (nrow(dataMatNAFiltered) == 0) {
            cli::cli_abort(
              c(
                "No rows remained after NA-filtering for comparison {.val {comp}}.",
                i = "{.arg condCol}: {.val {condCol}}",
                i = "{.arg batchCol}: {.val {batchCol}}",
                i = "Consider reducing {.arg leastRepCount} (lower limit for non-missing values within each condition-level combination)."
              ),
              class = "normalyzerde_error",
              call = NULL
            )
          }

          naFilterContrast <- rownames(dataMat) %in% rownames(dataMatNAFiltered)

          statResults <- calculateContrastStatistics(
            dataMatNAFiltered,
            as.character(designDfComp[, condCol]),
            compLevels,
            designDfCurrent = designDfComp,
            backendKey = comp
          )

          compLists <- storeContrastResults(
            compLists,
            statResults,
            comp,
            naFilterContrast
          )
        }
      }
    }

    pairwiseCompsP(nst) <- compLists[["P"]]
    pairwiseCompsFdr(nst) <- compLists[["FDR"]]
    pairwiseCompsAve(nst) <- compLists[["Ave"]]
    pairwiseCompsFold(nst) <- compLists[["Fold"]]

    if (!is.null(limpaBackend)) {
      existing <- backendData(nst)
      existing[["limpa"]] <- limpaBackend
      backendData(nst) <- existing
    }

    nst
  }
)


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
verifyContrasts <- function(designLevels, contrasts, splitter = "-") {
  for (contrast in contrasts) {
    parts <- unlist(strsplit(contrast, splitter, fixed = TRUE))

    if (length(parts) != 2) {
      cli::cli_abort(
        c(
          "A contrast string delimited by one splitter ({.val {splitter}}) was expected.",
          i = "Got: {.val {contrast}}."
        ),
        class = "normalyzerde_error",
        call = NULL
      )
    }

    if (!all(parts %in% designLevels)) {
      cli::cli_abort(
        c(
          "There were issues in your contrast.",
          i = "All contrasts: {.val {paste(contrasts, collapse = \", \")}}.",
          i = "Problematic contrast: {.val {contrast}}.",
          i = "Levels present in design: {.val {paste(unique(designLevels), collapse = \", \")}}."
        ),
        class = "normalyzerde_error",
        call = NULL
      )
    }
  }
}

.getContrastSplitter <- function(nst, default = "-") {
  if (!methods::is(nst, "NormalyzerStatistics")) {
    return(default)
  }

  sep <- tryCatch(contrastSplitter(nst), error = function(e) character())
  if (length(sep) == 0 || is.na(sep[1]) || !nzchar(sep[1])) {
    return(default)
  }

  sep[1]
}

setupModelFromDesign <- function(
  designDf,
  condCol,
  batchCol = NULL,
  type = "limma"
) {
  if (is.null(batchCol)) {
    Variable <- as.factor(designDf[, condCol])
    model <- ~ 0 + Variable
  } else {
    if (!(type %in% c("limma", "limma_intensity", "limpa"))) {
      cli::cli_abort(
        "Batch compensation is only compatible with {.val limma} and {.val limpa}, got {.val {type}}.",
        class = "normalyzerde_error",
        call = NULL
      )
    }
    Variable <- as.factor(designDf[, condCol])
    Batch <- as.factor(designDf[, batchCol])
    model <- ~ 0 + Variable + Batch
  }
  model
}

parseContrastLevels <- function(comp, splitter) {
  compSplit <- unlist(strsplit(comp, splitter, fixed = TRUE))

  if (length(compSplit) != 2) {
    cli::cli_abort(
      c(
        "{.arg comparisons} entries must be in format {.val cond1}{.val {splitter}}{.val cond2}.",
        i = "Split result: {.val {paste(compSplit, collapse = ' ')}}."
      ),
      class = "normalyzerde_error",
      call = NULL
    )
  }

  compSplit
}

assertContrastLevelsPresent <- function(levels, sampleGroups) {
  level1 <- levels[1]
  level2 <- levels[2]

  if (!any(sampleGroups %in% level1)) {
    cli::cli_abort(
      c(
        "No samples matching condition {.val {level1}}.",
        i = "Conditions present: {.val {paste(sampleGroups, collapse = ' ')}}."
      ),
      class = "normalyzerde_error",
      call = NULL
    )
  }

  if (!any(sampleGroups %in% level2)) {
    cli::cli_abort(
      c(
        "No samples matching condition {.val {level2}}.",
        i = "Conditions present: {.val {paste(sampleGroups, collapse = ' ')}}."
      ),
      class = "normalyzerde_error",
      call = NULL
    )
  }

  invisible(levels)
}

initializeContrastResultLists <- function(statMeasures) {
  out <- stats::setNames(vector("list", length(statMeasures)), statMeasures)
  for (statMeasure in statMeasures) {
    out[[statMeasure]] <- list()
  }
  out
}

storeContrastResults <- function(
  compLists,
  statResults,
  compName,
  naFilterContrast,
  statMeasures = names(compLists)
) {
  for (statMeasure in statMeasures) {
    compLists[[statMeasure]][[compName]] <- c()
    compLists[[statMeasure]][[compName]][naFilterContrast] <- statResults[[
      statMeasure
    ]]
  }

  compLists
}

sanitizeLimmaDesign <- function(limmaDesign) {
  rawNames <- colnames(limmaDesign)
  safeNames <- base::make.names(rawNames, unique = TRUE)
  colnames(limmaDesign) <- safeNames
  list(design = limmaDesign, coefMap = stats::setNames(safeNames, rawNames))
}

chooseOneVsRestLabel <- function(
  existingLabels,
  candidates = c("rest", "others", "all_other", "all_others")
) {
  if (length(candidates) == 0) {
    cli::cli_abort(
      "Expected at least one candidate label.",
      class = "normalyzerde_error",
      call = NULL
    )
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

  doTTest <- function(c1Vals, c2Vals, default = NA) {
    if (
      length(stats::na.omit(c1Vals)) > 1 &&
        length(stats::na.omit(c2Vals)) > 1
    ) {
      tryCatch(stats::t.test(c1Vals, c2Vals)[[3]], error = function(x) NA)
    } else {
      NA
    }
  }

  welchPValCol <- apply(
    dataMat,
    1,
    function(row) doTTest(row[s1cols], row[s2cols], default = NA)
  )
  welchFDRCol <- stats::p.adjust(welchPValCol, method = "BH")

  statResults <- list()

  statResults[["P"]] <- welchPValCol
  statResults[["FDR"]] <- welchFDRCol
  statResults[["Ave"]] <- rowMeans(dataMat, na.rm = TRUE)
  statResults[["Fold"]] <- rowMeans(
    dataMat[, s1cols, drop = FALSE],
    na.rm = TRUE
  ) -
    rowMeans(dataMat[, s2cols, drop = FALSE], na.rm = TRUE)

  statResults
}

calculateLimmaContrast <- function(
  dataMat,
  limmaDesign,
  limmaFit,
  levels,
  useIntensityTrend,
  coefMap = NULL
) {
  coefLevel1 <- paste0("Variable", levels[1])
  coefLevel2 <- paste0("Variable", levels[2])

  if (!is.null(coefMap)) {
    if (!(coefLevel1 %in% names(coefMap))) {
      cli::cli_abort(
        c(
          "Could not find limma coefficient name for level {.val {levels[1]}}.",
          i = "Expected: {.val {coefLevel1}}."
        ),
        class = "normalyzerde_error",
        call = NULL
      )
    }
    if (!(coefLevel2 %in% names(coefMap))) {
      cli::cli_abort(
        c(
          "Could not find limma coefficient name for level {.val {levels[2]}}.",
          i = "Expected: {.val {coefLevel2}}."
        ),
        class = "normalyzerde_error",
        call = NULL
      )
    }

    coefLevel1 <- coefMap[[coefLevel1]]
    coefLevel2 <- coefMap[[coefLevel2]]
  }

  myContrast <- paste0(coefLevel1, "-", coefLevel2)
  contrastMatrix <- limma::makeContrasts(
    contrasts = c(myContrast),
    levels = limmaDesign
  )
  fitContrasts <- limma::contrasts.fit(limmaFit, contrastMatrix)
  fitBayes <- limma::eBayes(fitContrasts, trend = useIntensityTrend)
  limmaTable <- limma::topTable(
    fitBayes,
    coef = 1,
    number = Inf,
    sort.by = "none"
  )

  statResults <- list()
  statResults[["P"]] <- limmaTable$P.Value
  statResults[["FDR"]] <- limmaTable$adj.P.Val
  statResults[["Ave"]] <- limmaTable$AveExpr
  statResults[["Fold"]] <- limmaTable$logFC
  statResults
}
