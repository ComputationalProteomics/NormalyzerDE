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
#' @slot comparisons Spot for saving vector of last used contrasts
#' @slot condCol Column containing last used conditions
#' @slot batchCol Column containing last used batch conditions
#' @slot splitter Character dividing contrast conditions
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
    batchCol = "numeric",
    splitter = "character"
  )
)

#' Constructor for NormalyzerStatistics
#'
#' @param experimentObj Instance of SummarizedExperiment containing matrix
#'   and design information as column data
#' @param logTrans Whether the input data should be log transformed. When
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
      warning(
        "Non-finite values produced by log2 transform (e.g. zeros or negative values) ",
        "were treated as missing (set to NA)."
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
#' the completed expression matrix (as a limma \code{EList}) and/or the fitted
#' \code{MArrayLM} object(s) when \code{limpaKeep} is enabled.
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
#' @param splitter Character dividing contrast conditions.
#' @param type Type of statistical test ("limma", "limma_intensity", "welch" or
#'   "limpa"). "limpa" uses the optional Bioconductor package \pkg{limpa} to
#'   handle missing values via a detection probability curve (DPC) model.
#' @param leastRepCount Least replicates in each group to be retained for
#'   contrast calculations. For \code{type="limpa"}, a feature is retained if at
#'   least one group has \code{leastRepCount} observed samples (and features
#'   entirely missing across all samples are removed).
#' @param impute Whether to impute values (ignored for \code{type="limpa"}).
#' @param imputeMinFraction Minimum fraction non-NA values for an analyte in any
#'   group to impute in other groups (ignored for \code{type="limpa"}).
#' @param subsetByComparison If TRUE, subset data and design to each comparison
#'   before NA-filtering, imputation and model fitting.
#' @param oneVsRest If TRUE, compute one-vs-rest contrasts for each group in
#'   \code{condCol} (or the subset in \code{oneVsRestGroups}).
#' @param oneVsRestGroups Optional character vector specifying which groups in
#'   \code{condCol} to compare against all other samples.
#' @param limpaProteinIdCol For \code{type="limpa"}, optionally summarize
#'   peptide/precursor rows to protein-level using \code{limpa::dpcQuant()}.
#'   Set to a column name in the row annotation (for example \code{"Protein.Group"})
#'   to use as the protein identifier. Use \code{"auto"} (default) to try common
#'   identifiers. If the chosen column contains duplicate identifiers, the data
#'   are summarized once across all samples and the output rows correspond to
#'   proteins. Set to \code{NULL} to disable protein summarization and treat each
#'   row as one feature (recommended for PTM-level data such as phosphoproteomics
#'   where each row corresponds to a modified site).
#' @param limpaByRow For \code{type="limpa"}, treat each input row as a separate
#'   feature and always use \code{limpa::dpcQuantByRow()} instead of summarizing
#'   via \code{limpa::dpcQuant()}. This is recommended for PTM-level matrices
#'   (e.g., phosphosites). Equivalent to setting \code{limpaProteinIdCol=NULL}.
#' @param limpaDpc For \code{type="limpa"}, optional DPC parameters to pass to
#'   \code{limpa::dpcQuant()} / \code{limpa::dpcQuantByRow()}. Can be a list as
#'   returned by \code{limpa::dpc()}, or a numeric vector \code{c(beta0, beta1)}.
#' @param limpaDpcMethod For \code{type="limpa"}, optional method to estimate the
#'   DPC parameters from the data when \code{limpaDpc} is not supplied.
#'   \code{"none"} (default) uses a fixed slope (\code{limpaQuantArgs$dpc.slope},
#'   default \code{0.8}) and lets limpa estimate the intercept internally.
#'   \code{"dpc"} estimates both DPC parameters from the observed-normal model
#'   via \code{limpa::dpc()}. \code{"dpcCN"} estimates the DPC from the
#'   complete-normal model via \code{limpa::dpcCN()}, which can be more robust
#'   for datasets with very large fold-changes.
#' @param limpaDpcArgs For \code{type="limpa"}, optional named list of additional
#'   arguments forwarded to \code{limpa::dpc()} or \code{limpa::dpcCN()} when
#'   \code{limpaDpcMethod} is not \code{"none"}. Argument \code{y} is ignored.
#'   For \code{limpaDpcMethod="dpcCN"}, \code{dpc.slope.start} defaults to
#'   \code{limpaQuantArgs$dpc.slope}.
#' @param limpaQuantArgs For \code{type="limpa"}, optional named list of
#'   additional arguments forwarded to \code{limpa::dpcQuant()} /
#'   \code{limpa::dpcQuantByRow()}. Use this to set \code{dpc.slope} (default
#'   \code{0.8}), \code{chunk} (default \code{1000L}), and \code{verbose} (default
#'   \code{FALSE}), plus any additional \code{...} arguments supported by limpa.
#'   Arguments \code{y}, \code{protein.id}, and \code{dpc} are ignored.
#' @param limpaQuantifiedRds For \code{type="limpa"}, optional path to an RDS
#'   file containing a quantified \code{EList} object (as returned by
#'   \code{limpa::dpcQuant()} or \code{limpa::dpcQuantByRow()}). When provided,
#'   NormalyzerDE skips the internal \code{dpcQuant*()} step and reuses the
#'   cached quantification (including \code{standard.error} and
#'   \code{n.observations}) for differential expression. This is required to
#'   preserve quantification uncertainty when you first ran \code{dpcQuant*()}
#'   upstream (for example via \code{normalyzer(preQuant=\"limpa\")}).
#' @param limpaPostQuantNorm For \code{type="limpa"}, optional between-sample
#'   normalization applied to the quantified expression matrix after
#'   \code{limpa::dpcQuant()} / \code{limpa::dpcQuantByRow()} and before
#'   \code{limpa::dpcDE()}. One of \code{"none"} (default), \code{"GI"},
#'   \code{"median"}, \code{"mean"}, \code{"Quantile"} (or \code{"quantile"}),
#'   \code{"CycLoess"}, or \code{"RLR"}.
#'   Avoid double-normalization if your input was already normalized upstream.
#' @param limpaDEArgs For \code{type="limpa"}, optional named list of additional
#'   arguments forwarded to \code{limpa::dpcDE()} (and then to
#'   \code{limpa::voomaLmFitWithImputation()}). To enable limma sample weights,
#'   set \code{limpaDEArgs = list(sample.weights = TRUE)}.
#' @param limpaKeep For \code{type="limpa"}, optionally store intermediate limpa
#'   objects in \code{backendData(nst)$limpa} for reuse/debugging. Set to
#'   \code{"elist"} to keep the \code{EList} object (completed expression matrix
#'   plus uncertainty estimates), \code{"fit"} to keep the fitted \code{MArrayLM}
#'   object(s), or \code{"all"} to keep both. Default is \code{"none"}.
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
    limpaProteinIdCol = "auto",
    limpaDpc = NULL,
    limpaDpcMethod = c("none", "dpc", "dpcCN"),
    limpaDpcArgs = NULL,
    limpaQuantArgs = NULL,
    limpaQuantifiedRds = NULL,
    limpaDEArgs = NULL,
    limpaKeep = c("none", "elist", "fit", "all"),
    limpaByRow = FALSE,
    limpaPostQuantNorm = c(
      "none",
      "GI",
      "median",
      "mean",
      "Quantile",
      "CycLoess",
      "RLR",
      "quantile"
    )
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
    limpaProteinIdCol = "auto",
    limpaDpc = NULL,
    limpaDpcMethod = c("none", "dpc", "dpcCN"),
    limpaDpcArgs = NULL,
    limpaQuantArgs = NULL,
    limpaQuantifiedRds = NULL,
    limpaDEArgs = NULL,
    limpaKeep = c("none", "elist", "fit", "all"),
    limpaByRow = FALSE,
    limpaPostQuantNorm = c(
      "none",
      "GI",
      "median",
      "mean",
      "Quantile",
      "CycLoess",
      "RLR",
      "quantile"
    )
  ) {
    dataMat <- dataMat(nst)
    dataMatOriginalRowNames <- rownames(dataMat)
    designDf <- designDf(nst)

    contrastSplitter(nst) <- splitter
    condCol(nst) <- as.character(designDf[, condCol])

    if (
      !is.null(batchCol) && !(type %in% c("limma", "limma_intensity", "limpa"))
    ) {
      stop(
        "Batch compensation only compatible with Limma and limpa, got: ",
        type
      )
    }

    if (!is.null(batchCol)) {
      conditionCombs <- paste(
        designDf[, condCol],
        designDf[, batchCol],
        sep = "_"
      )
      batchCol(nst) <- as.factor(designDf[, batchCol])
    } else {
      conditionCombs <- designDf[, condCol]
      batchCol(nst) <- numeric()
    }

    rownames(dataMat) <- seq_len(nrow(dataMat))

    sampleReplicateGroupsStrings <- as.character(designDf[, condCol])
    statMeasures <- c("P", "FDR", "Ave", "Fold")

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
          stop(
            "Batch compensation only compatible with Limma and limpa, got: ",
            type
          )
        }
        Variable <- as.factor(designDf[, condCol])
        Batch <- as.factor(designDf[, batchCol])
        model <- ~ 0 + Variable + Batch
      }
      model
    }

    warnIfLimpaInputLooksNotLog2 <- function(dataMat, threshold = 50) {
      finiteVals <- dataMat[is.finite(dataMat)]
      if (length(finiteVals) == 0) {
        return(invisible(NULL))
      }

      maxVal <- max(finiteVals)
      if (is.finite(maxVal) && maxVal > threshold) {
        warning(
          "For type='limpa', the input data should be on the log2 scale (missing values as NA). ",
          "The values look large for log2 data (max finite value = ",
          format(signif(maxVal, 4), trim = TRUE),
          "). If your matrix is on the linear scale, set `logTrans=TRUE` in ",
          "`normalyzerDE()` (or log2-transform upstream before calling `calculateContrasts()`).",
          call. = FALSE
        )
      }

      invisible(NULL)
    }

    limpaPostQuantNormUse <- "none"
    if (type == "limpa") {
      requireLimpaPackageInternal("Statistics type 'limpa'")
      warnIfLimpaInputLooksNotLog2(dataMat)

      limpaDpcMethod <- match.arg(limpaDpcMethod)
      limpaKeep <- match.arg(limpaKeep)
      limpaPostQuantNormUse <- match.arg(limpaPostQuantNorm)

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
          stop(
            "limpaQuantifiedRds must be a length-1 character file path (or NULL)."
          )
        }
        limpaQuantifiedRds <- as.character(limpaQuantifiedRds)[1]
        if (is.na(limpaQuantifiedRds) || !nzchar(limpaQuantifiedRds)) {
          stop("limpaQuantifiedRds must be a non-empty file path (or NULL).")
        }
        if (!file.exists(limpaQuantifiedRds)) {
          stop("limpaQuantifiedRds file does not exist: ", limpaQuantifiedRds)
        }
      }

      validateLimpaDpc(limpaDpc)
    }

    if (type != "limpa" && !is.null(limpaQuantifiedRds)) {
      stop("limpaQuantifiedRds is only supported for type='limpa'.")
    }

    limpaDpcArgsUse <- if (type == "limpa") {
      sanitizeLimpaDpcArgs(limpaDpcArgs)
    } else {
      list()
    }

    limpaQuantArgsUse <- if (type == "limpa") {
      applyLimpaQuantDefaultsAndValidate(sanitizeLimpaQuantArgs(limpaQuantArgs))
    } else {
      list()
    }

    limpaVerboseUse <- if (type == "limpa") {
      isTRUE(limpaQuantArgsUse[["verbose"]])
    } else {
      FALSE
    }
    limpaDpcSlopeUse <- if (type == "limpa") {
      as.numeric(limpaQuantArgsUse[["dpc.slope"]])
    } else {
      0.8
    }

    if (type == "limpa" && leastRepCount > 1 && isTRUE(limpaVerboseUse)) {
      message(
        "For type='limpa', leastRepCount>1 can drop informative sparse features. ",
        "Consider leastRepCount=1 unless you explicitly want to filter sparse rows."
      )
    }

    limpaDpcUse <- if (type == "limpa") {
      if (!is.null(limpaDpc)) {
        if (limpaDpcMethod != "none" && isTRUE(limpaVerboseUse)) {
          message(
            "limpaDpc was supplied; ignoring limpaDpcMethod='",
            limpaDpcMethod,
            "'."
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

    if (type == "limpa" && limpaKeep != "none") {
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
        elists = list()
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

      if (length(proteinId) != nrow(dataMat)) {
        stop(
          "Row annotation column '",
          proteinIdCol,
          "' does not match the number of rows in the data matrix."
        )
      }

      if (anyNA(proteinId) || any(proteinId == "")) {
        stop(
          "Row annotation column '",
          proteinIdCol,
          "' contains missing or empty protein identifiers. ",
          "Remove these rows or choose another column."
        )
      }

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
      stableCols <- inferStableProteinAnnotationCols(
        genesDf = genesInput,
        proteinId = proteinId,
        proteinIdCol = proteinIdCol
      )
      genesInput <- genesInput[,
        unique(c(proteinIdCol, stableCols)),
        drop = FALSE
      ]
      yPeptide <- methods::new("EList", list(E = dataMat, genes = genesInput))
      yProtein <- do.call(
        limpa::dpcQuant,
        c(
          list(
            y = yPeptide,
            protein.id = proteinIdCol,
            dpc = limpaDpcUse
          ),
          limpaQuantArgsUse
        )
      )

      proteinIds <- rownames(yProtein$E)
      rowIds <- as.character(seq_len(nrow(yProtein$E)))

      rownames(yProtein$E) <- rowIds

      if (!is.null(yProtein$other$n.observations)) {
        rownames(yProtein$other$n.observations) <- rowIds
      }
      if (!is.null(yProtein$other$standard.error)) {
        rownames(yProtein$other$standard.error) <- rowIds
      }

      genes <- if (!is.null(yProtein$genes)) {
        as.data.frame(yProtein$genes, check.names = FALSE)
      } else {
        data.frame(check.names = FALSE)
      }
      if (!(proteinIdCol %in% colnames(genes))) {
        genes[[proteinIdCol]] <- proteinIds
      }
      genes <- genes[,
        c(proteinIdCol, setdiff(names(genes), proteinIdCol)),
        drop = FALSE
      ]
      rownames(genes) <- rowIds
      yProtein$genes <- genes

      yProtein
    }

    filterLowRepLimpa <- function(
      df,
      groups,
      leastRep = 1,
      nObservations = NULL
    ) {
      obsMat <- if (!is.null(nObservations)) {
        if (!is.matrix(nObservations)) {
          stop("nObservations must be a matrix when provided.")
        }
        if (!all(dim(nObservations) == dim(df))) {
          stop("nObservations must have the same dimensions as df.")
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
        warning(
          "limpaProteinIdCol '",
          proteinIdCol,
          "' contains no duplicated identifiers, so no peptide/precursor-to-protein summarization ",
          "will be performed. The analysis will proceed with the input rows as features (this may ",
          "indicate the input is already protein-level).",
          call. = FALSE
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

      if (
        length(dataRowIds) == length(cachedRowIds) &&
          anyDuplicated(cachedRowIds) == 0
      ) {
        warning(
          "Could not match limpaQuantifiedRds row identifiers by name; assuming row order matches the input matrix.",
          call. = FALSE
        )
        return(stats::setNames(cachedRowIds, dataRowIds))
      }

      stop(
        "limpaQuantifiedRds row identifiers could not be matched to the data matrix. ",
        "Ensure the quantified EList was generated from the same matrix (rows), ",
        "or provide matching row names."
      )
    }

    limpaCachedRowMap <- NULL
    limpaQuantifiedCached <- NULL
    if (type == "limpa" && !is.null(limpaQuantifiedRds)) {
      limpaQuantifiedCached <- readRDS(limpaQuantifiedRds)
      if (!inherits(limpaQuantifiedCached, "EList")) {
        stop(
          "limpaQuantifiedRds must contain a limma EList object, got: ",
          paste(class(limpaQuantifiedCached), collapse = ", ")
        )
      }
      if (
        is.null(limpaQuantifiedCached$E) || !is.matrix(limpaQuantifiedCached$E)
      ) {
        stop("limpaQuantifiedRds EList must contain a matrix element 'E'.")
      }
      if (is.null(colnames(limpaQuantifiedCached$E))) {
        stop("limpaQuantifiedRds EList$E must contain sample column names.")
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
        stop(
          "limpaQuantifiedRds EList$E must have unique, non-empty row names."
        )
      }
      rownames(limpaQuantifiedCached$E) <- cachedRowIds

      if (
        is.null(limpaQuantifiedCached$other$n.observations) ||
          !is.matrix(limpaQuantifiedCached$other$n.observations)
      ) {
        stop(
          "limpaQuantifiedRds EList must contain other$n.observations (matrix). ",
          "Use an EList returned by limpa::dpcQuant() or limpa::dpcQuantByRow()."
        )
      }
      if (
        !all(
          dim(limpaQuantifiedCached$other$n.observations) ==
            dim(limpaQuantifiedCached$E)
        )
      ) {
        stop(
          "limpaQuantifiedRds other$n.observations must have the same dimensions as EList$E."
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
        stop(
          "limpaQuantifiedRds EList must contain other$standard.error (matrix). ",
          "Use an EList returned by limpa::dpcQuant() or limpa::dpcQuantByRow()."
        )
      }
      if (
        !all(
          dim(limpaQuantifiedCached$other$standard.error) ==
            dim(limpaQuantifiedCached$E)
        )
      ) {
        stop(
          "limpaQuantifiedRds other$standard.error must have the same dimensions as EList$E."
        )
      }
      rownames(limpaQuantifiedCached$other$standard.error) <- cachedRowIds
      if (is.null(colnames(limpaQuantifiedCached$other$standard.error))) {
        colnames(limpaQuantifiedCached$other$standard.error) <- colnames(
          limpaQuantifiedCached$E
        )
      }

      if (anyNA(dataMat)) {
        stop(
          "limpaQuantifiedRds requires a completed expression matrix without NA values. ",
          "Use the post-quant matrix (for example the 'log2' output from normalyzer(preQuant='limpa'))."
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
        stop(
          "limpaQuantifiedRds is missing sample columns required by the data matrix: ",
          paste(missingCols, collapse = ", ")
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
        stop(
          "limpaPostQuantNorm='",
          limpaPostQuantNormUse,
          "' requires a completed expression matrix without NA values."
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
        stop("Unknown limpaPostQuantNorm value: ", limpaPostQuantNormUse)
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
          limpa::dpcQuantByRow,
          c(
            list(
              y = dataMat,
              dpc = limpaDpcUse
            ),
            limpaQuantArgsUse
          )
        )
        yImputed <- applyLimpaPostQuantNorm(yImputed)

        limpaDEArgsUse <-
          applyLimpaDEDefaultsAndValidate(sanitizeLimpaDEArgs(limpaDEArgs))
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
        stop(
          "Failed to map data matrix rows to limpaQuantifiedRds rows. ",
          "Please ensure both inputs come from the same quantified matrix."
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
          stop(
            "limpaQuantifiedRds genes rows could not be aligned to EList rows."
          )
        }
        rownames(genesDf) <- rows
        yUse$genes <- genesDf
      }

      if (!is.null(limpaQuantifiedCached)) {
        if (anyNA(dataMat)) {
          stop(
            "limpaQuantifiedRds requires a completed expression matrix without NA values. ",
            "Use the post-quant matrix (for example the 'log2' output from normalyzer(preQuant='limpa'))."
          )
        }
        yUse$E <- dataMat
      }

      yUse <- applyLimpaPostQuantNorm(yUse)

      limpaDEArgsUse <-
        applyLimpaDEDefaultsAndValidate(sanitizeLimpaDEArgs(limpaDEArgs))
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
        stop(
          "Failed to map filtered data rows to limpaQuantifiedRds n.observations rows."
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

      missingGroups <- setdiff(
        targetGroups,
        unique(sampleReplicateGroupsStrings)
      )
      if (length(missingGroups) > 0) {
        stop(
          "Some groups in oneVsRestGroups were not found in condCol '",
          condCol,
          "': ",
          paste(missingGroups, collapse = ", ")
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

        dataMatNAFiltered <- if (type == "limpa") {
          filterLowRepLimpa(
            dataMat,
            conditionCombsOVR,
            leastRep = leastRepCount,
            nObservations = subsetLimpaNObservations(dataMat)
          )
        } else {
          filterLowRep(
            dataMat,
            conditionCombsOVR,
            leastRep = leastRepCount
          )
        }

        if (nrow(dataMatNAFiltered) == 0) {
          stop(
            "No rows remained after NA-filtering for one-vs-rest comparison: '",
            compName,
            "' (condCol: '",
            condCol,
            "', batchCol: '",
            batchCol,
            "' (if empty then batchCol is not specified))\n",
            "Consider whether you can reduce the 'leastRepCount' setting which sets the lower limit ",
            "of number of NA values in each condition-level combination"
          )
        }

        naFilterContrast <- rownames(dataMat) %in% rownames(dataMatNAFiltered)

        if (type != "limpa" && leastRepCount == 0 && impute) {
          dataMatNAFiltered <- imputeGroupValues(
            dataMatNAFiltered,
            conditionCombsOVR,
            minFraction = imputeMinFraction
          )
        }

        if (type == "welch") {
          statResults <- calculateWelch(
            dataMatNAFiltered,
            groupHeader,
            c(groupLabel, restLabel)
          )
        } else if (type %in% c("limma", "limma_intensity", "limpa")) {
          model <- setupModelFromDesign(
            designDfOVR,
            condCol,
            batchCol = batchCol,
            type = type
          )
          limmaDesignRaw <- stats::model.matrix(model)
          limmaPrepared <- sanitizeLimmaDesign(limmaDesignRaw)
          limmaDesign <- limmaPrepared$design
          limmaCoefMap <- limmaPrepared$coefMap
          limmaFit <- if (type == "limpa") {
            calculateLimpaFit(
              dataMatNAFiltered,
              limmaDesign,
              backendKey = compName
            )
          } else {
            limma::lmFit(dataMatNAFiltered, limmaDesign)
          }

          statResults <- calculateLimmaContrast(
            dataMatNAFiltered,
            limmaDesign,
            limmaFit,
            c(groupLabel, restLabel),
            useIntensityTrend = type == "limma_intensity",
            coefMap = limmaCoefMap
          )
        } else {
          stop("Unknown statistics type: ", type)
        }

        for (statMeasure in statMeasures) {
          compLists[[statMeasure]][[compName]] <- c()
          compLists[[statMeasure]][[compName]][
            naFilterContrast
          ] <- statResults[[statMeasure]]
        }
      }
    } else {
      if (is.null(comparisons)) {
        stop("Argument 'comparisons' must be provided unless oneVsRest=TRUE.")
      }

      comparisons <- as.character(comparisons)
      comparisons(nst) <- comparisons
      verifyContrasts(
        sampleReplicateGroupsStrings,
        comparisons,
        splitter = splitter
      )

      if (!subsetByComparison) {
        dataMatNAFiltered <- if (type == "limpa") {
          filterLowRepLimpa(
            dataMat,
            conditionCombs,
            leastRep = leastRepCount,
            nObservations = subsetLimpaNObservations(dataMat)
          )
        } else {
          filterLowRep(
            dataMat,
            conditionCombs,
            leastRep = leastRepCount
          )
        }

        if (nrow(dataMatNAFiltered) == 0) {
          stop(
            "No rows remained after NA-filtering for condition: '",
            condCol,
            "' and batchCol: '",
            batchCol,
            "' (if empty then batchCol is not specified)\n",
            "Consider whether you can reduce the 'leastRepCount' setting which sets the lower limit ",
            "of number of NA values in each condition-level combination ",
            "You could also try running without batchCol and see if there is enough data per condition then"
          )
        }

        naFilterContrast <- rownames(dataMat) %in% rownames(dataMatNAFiltered)

        if (type != "limpa" && leastRepCount == 0 && impute) {
          dataMatNAFiltered <- imputeGroupValues(
            dataMatNAFiltered,
            conditionCombs,
            minFraction = imputeMinFraction
          )
        }

        model <- setupModel(nst, condCol, batchCol = batchCol, type = type)

        if (type %in% c("limma", "limma_intensity", "limpa")) {
          limmaDesignRaw <- stats::model.matrix(model)
          limmaPrepared <- sanitizeLimmaDesign(limmaDesignRaw)
          limmaDesign <- limmaPrepared$design
          limmaCoefMap <- limmaPrepared$coefMap
          limmaFit <- if (type == "limpa") {
            calculateLimpaFit(dataMatNAFiltered, limmaDesign)
          } else {
            limma::lmFit(dataMatNAFiltered, limmaDesign)
          }
        }

        for (comp in comparisons) {
          compSplit <- unlist(strsplit(comp, splitter, fixed = TRUE))

          if (length(compSplit) != 2) {
            stop(
              "Comparison should be in format cond1-cond2 ",
              "here the split product was: ",
              paste(compSplit, collapse = " ")
            )
          }

          level1 <- compSplit[1]
          level2 <- compSplit[2]

          if (!any(sampleReplicateGroupsStrings %in% level1)) {
            stop(
              "No samples matching condition ",
              level1,
              " found in conditions: ",
              paste(sampleReplicateGroupsStrings, collapse = " ")
            )
          }

          if (!any(sampleReplicateGroupsStrings %in% level2)) {
            stop(
              "No samples matching condition ",
              level2,
              " found in conditions: ",
              paste(sampleReplicateGroupsStrings, collapse = " ")
            )
          }

          if (type == "welch") {
            statResults <- calculateWelch(
              dataMatNAFiltered,
              sampleReplicateGroupsStrings,
              c(level1, level2)
            )
          } else if (type %in% c("limma", "limma_intensity", "limpa")) {
            statResults <- calculateLimmaContrast(
              dataMatNAFiltered,
              limmaDesign,
              limmaFit,
              c(level1, level2),
              useIntensityTrend = type == "limma_intensity",
              coefMap = limmaCoefMap
            )
          } else {
            stop("Unknown statistics type: ", type)
          }

          for (statMeasure in statMeasures) {
            compLists[[statMeasure]][[comp]] <- c()
            compLists[[statMeasure]][[comp]][naFilterContrast] <- statResults[[
              statMeasure
            ]]
          }
        }
      } else {
        for (comp in comparisons) {
          compSplit <- unlist(strsplit(comp, splitter, fixed = TRUE))

          if (length(compSplit) != 2) {
            stop(
              "Comparison should be in format cond1-cond2 ",
              "here the split product was: ",
              paste(compSplit, collapse = " ")
            )
          }

          level1 <- compSplit[1]
          level2 <- compSplit[2]

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

          dataMatNAFiltered <- if (type == "limpa") {
            filterLowRepLimpa(
              dataMatComp,
              conditionCombsComp,
              leastRep = leastRepCount,
              nObservations = subsetLimpaNObservations(dataMatComp)
            )
          } else {
            filterLowRep(
              dataMatComp,
              conditionCombsComp,
              leastRep = leastRepCount
            )
          }

          if (nrow(dataMatNAFiltered) == 0) {
            stop(
              "No rows remained after NA-filtering for comparison: '",
              comp,
              "' (condCol: '",
              condCol,
              "', batchCol: '",
              batchCol,
              "' (if empty then batchCol is not specified))\n",
              "Consider whether you can reduce the 'leastRepCount' setting which sets the lower limit ",
              "of number of NA values in each condition-level combination"
            )
          }

          naFilterContrast <- rownames(dataMat) %in% rownames(dataMatNAFiltered)

          if (type != "limpa" && leastRepCount == 0 && impute) {
            dataMatNAFiltered <- imputeGroupValues(
              dataMatNAFiltered,
              conditionCombsComp,
              minFraction = imputeMinFraction
            )
          }

          if (type == "welch") {
            statResults <- calculateWelch(
              dataMatNAFiltered,
              as.character(designDfComp[, condCol]),
              c(level1, level2)
            )
          } else if (type %in% c("limma", "limma_intensity", "limpa")) {
            model <- setupModelFromDesign(
              designDfComp,
              condCol,
              batchCol = batchCol,
              type = type
            )
            limmaDesignRaw <- stats::model.matrix(model)
            limmaPrepared <- sanitizeLimmaDesign(limmaDesignRaw)
            limmaDesign <- limmaPrepared$design
            limmaCoefMap <- limmaPrepared$coefMap
            limmaFit <- if (type == "limpa") {
              calculateLimpaFit(
                dataMatNAFiltered,
                limmaDesign,
                backendKey = comp
              )
            } else {
              limma::lmFit(dataMatNAFiltered, limmaDesign)
            }

            statResults <- calculateLimmaContrast(
              dataMatNAFiltered,
              limmaDesign,
              limmaFit,
              c(level1, level2),
              useIntensityTrend = type == "limma_intensity",
              coefMap = limmaCoefMap
            )
          } else {
            stop("Unknown statistics type: ", type)
          }

          for (statMeasure in statMeasures) {
            compLists[[statMeasure]][[comp]] <- c()
            compLists[[statMeasure]][[comp]][naFilterContrast] <- statResults[[
              statMeasure
            ]]
          }
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
      stop(
        "A contrast string delimited by one splitter (",
        splitter,
        ") was expected. Instead following was found: ",
        contrast
      )
    }

    if (!all(parts %in% designLevels)) {
      stop(
        "There were issues in your contrast. \n",
        "All contrasts: ",
        paste(contrasts, collapse = ", "),
        "\n",
        "Part with issue: ",
        contrast,
        "\n",
        "Not all parts was found in the design column levels. Levels present in design: \n",
        paste(unique(designLevels), collapse = ", ")
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

setupModel <- function(nst, condCol, batchCol = NULL, type = "limma") {
  if (is.null(batchCol)) {
    Variable <- as.factor(designDf(nst)[, condCol])
    model <- ~ 0 + Variable
  } else {
    if (!(type %in% c("limma", "limma_intensity", "limpa"))) {
      stop(
        "Batch compensation only compatible with Limma and limpa, got: ",
        type
      )
    }
    Variable <- as.factor(designDf(nst)[, condCol])
    Batch <- as.factor(designDf(nst)[, batchCol])
    model <- ~ 0 + Variable + Batch
  }
  model
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
      stop(
        "Could not find limma coefficient name for level '",
        levels[1],
        "' (expected: '",
        coefLevel1,
        "')"
      )
    }
    if (!(coefLevel2 %in% names(coefMap))) {
      stop(
        "Could not find limma coefficient name for level '",
        levels[2],
        "' (expected: '",
        coefLevel2,
        "')"
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
