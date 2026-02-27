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
#'   row as one protein (recommended for PTM-level data such as phosphoproteomics
#'   where each row corresponds to a modified site).
#' @param limpaByRow For \code{type="limpa"}, treat each input row as a separate
#'   protein and always use \code{limpa::dpcQuantByRow()} instead of summarizing
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
    limpaDEArgs = NULL,
    limpaKeep = c("none", "elist", "fit", "all"),
    limpaByRow = FALSE
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
    limpaDEArgs = NULL,
    limpaKeep = c("none", "elist", "fit", "all"),
    limpaByRow = FALSE
  ) {
    dataMat <- dataMat(nst)
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

    requireLimpaPackage <- function() {
      if (!requireNamespace("limpa", quietly = TRUE)) {
        stop(
          "Statistics type 'limpa' requires the optional Bioconductor package 'limpa'.\n",
          "Install it with `BiocManager::install(\"limpa\")`."
        )
      }
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

    if (type == "limpa") {
      requireLimpaPackage()
      warnIfLimpaInputLooksNotLog2(dataMat)

      limpaDpcMethod <- match.arg(limpaDpcMethod)
      limpaKeep <- match.arg(limpaKeep)

      limpaByRow <- as.logical(limpaByRow)[1]
      if (is.na(limpaByRow)) {
        stop("limpaByRow must be TRUE or FALSE.")
      }

      if (isTRUE(limpaByRow)) {
        proteinIdColValue <- if (is.null(limpaProteinIdCol)) {
          NULL
        } else {
          as.character(limpaProteinIdCol)[1]
        }
        if (!is.null(proteinIdColValue) && !identical(proteinIdColValue, "auto")) {
          stop(
            "limpaByRow=TRUE is incompatible with a non-default limpaProteinIdCol. ",
            "Set limpaProteinIdCol=NULL (or leave it as 'auto')."
          )
        }
        limpaProteinIdCol <- NULL
      }

      if (!is.null(limpaDpcArgs) && !is.list(limpaDpcArgs)) {
        stop("limpaDpcArgs must be a named list (or NULL).")
      }
      if (!is.null(limpaDpcArgs) && length(limpaDpcArgs) > 0) {
        if (is.null(names(limpaDpcArgs))) {
          stop("limpaDpcArgs must be a named list.")
        }
      }

      if (!is.null(limpaQuantArgs) && !is.list(limpaQuantArgs)) {
        stop("limpaQuantArgs must be a named list (or NULL).")
      }
      if (!is.null(limpaQuantArgs) && length(limpaQuantArgs) > 0) {
        if (is.null(names(limpaQuantArgs))) {
          stop("limpaQuantArgs must be a named list.")
        }
      }

      if (
        !is.null(limpaDpc) &&
          !(is.list(limpaDpc) ||
            (is.numeric(limpaDpc) && length(limpaDpc) == 2))
      ) {
        stop(
          "limpaDpc must be NULL, a list returned by limpa::dpc(), or a numeric vector c(beta0, beta1)."
        )
      }

      if (!is.null(limpaDEArgs) && !is.list(limpaDEArgs)) {
        stop("limpaDEArgs must be a list (or NULL).")
      }
    }

    sanitizeLimpaDpcArgs <- function(args) {
      if (is.null(args)) {
        return(list())
      }
      if (length(args) == 0) {
        return(list())
      }
      if (is.null(names(args))) {
        stop("limpaDpcArgs must be a named list.")
      }
      args[["y"]] <- NULL
      args
    }

    limpaDpcArgsUse <- if (type == "limpa") {
      sanitizeLimpaDpcArgs(limpaDpcArgs)
    } else {
      list()
    }

    sanitizeLimpaDEArgs <- function(args) {
      if (is.null(args)) {
        return(list())
      }
      if (length(args) == 0) {
        return(list())
      }
      if (is.null(names(args))) {
        stop("limpaDEArgs must be a named list.")
      }

      forbidden <- c("y", "design", "plot")
      args[forbidden] <- NULL
      args
    }

    applyLimpaDEDefaultsAndValidate <- function(args) {
      if (is.null(args[["sample.weights"]])) {
        args[["sample.weights"]] <- FALSE
      }
      sampleWeights <- as.logical(args[["sample.weights"]])[1]
      if (is.na(sampleWeights)) {
        stop("limpaDEArgs$sample.weights must be TRUE or FALSE.")
      }
      args[["sample.weights"]] <- sampleWeights
      args
    }

    sanitizeLimpaQuantArgs <- function(args) {
      if (is.null(args)) {
        return(list())
      }
      if (length(args) == 0) {
        return(list())
      }
      if (is.null(names(args))) {
        stop("limpaQuantArgs must be a named list.")
      }

      forbidden <- c("y", "protein.id", "dpc")
      args[forbidden] <- NULL
      args
    }

    applyLimpaQuantDefaultsAndValidate <- function(args) {
      if (!("dpc.slope" %in% names(args))) {
        args[["dpc.slope"]] <- 0.8
      }
      if (!("chunk" %in% names(args))) {
        args[["chunk"]] <- 1000L
      }
      if (!("verbose" %in% names(args))) {
        args[["verbose"]] <- FALSE
      }

      slope <- as.numeric(args[["dpc.slope"]])[1]
      if (is.na(slope) || !is.finite(slope) || slope <= 0) {
        stop("limpaQuantArgs$dpc.slope must be a single positive numeric value.")
      }
      args[["dpc.slope"]] <- slope

      chunk <- as.integer(args[["chunk"]])[1]
      if (is.na(chunk) || chunk < 1) {
        stop("limpaQuantArgs$chunk must be a positive integer.")
      }
      args[["chunk"]] <- chunk

      verbose <- as.logical(args[["verbose"]])[1]
      if (is.na(verbose)) {
        stop("limpaQuantArgs$verbose must be TRUE or FALSE.")
      }
      args[["verbose"]] <- verbose

      args
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

    validateLimpaArgsByFormals <- function(args, fn, methodLabel) {
      if (length(args) == 0) {
        return(args)
      }

      allowed <- names(formals(fn))
      allowed <- allowed[!is.na(allowed) & nzchar(allowed)]
      allowed <- setdiff(allowed, "y")

      unknown <- setdiff(names(args), allowed)
      if (length(unknown) > 0) {
        stop(
          "limpaDpcArgs contains unsupported argument names for limpaDpcMethod='",
          methodLabel,
          "': ",
          paste(unknown, collapse = ", "),
          ". Allowed: ",
          paste(allowed, collapse = ", "),
          "."
        )
      }

      args
    }

    estimateLimpaDpc <- function(dataMat, method = c("none", "dpc", "dpcCN")) {
      method <- match.arg(method)
      if (method == "none") {
        return(NULL)
      }

      dpcArgsUse <- limpaDpcArgsUse

      if (method == "dpc") {
        dpcArgsUse <- validateLimpaArgsByFormals(
          dpcArgsUse,
          fn = limpa::dpc,
          methodLabel = method
        )
        dpcCall <- c(list(y = dataMat), dpcArgsUse)
        if (!isTRUE(limpaVerboseUse)) {
          return(suppressMessages(do.call(limpa::dpc, dpcCall)))
        }
        return(do.call(limpa::dpc, dpcCall))
      }

      dpcArgsUse <- validateLimpaArgsByFormals(
        dpcArgsUse,
        fn = limpa::dpcCN,
        methodLabel = method
      )

      if (!("dpc.slope.start" %in% names(dpcArgsUse))) {
        dpcArgsUse[["dpc.slope.start"]] <- limpaDpcSlopeUse
      }
      if (!("verbose" %in% names(dpcArgsUse))) {
        dpcArgsUse[["verbose"]] <- isTRUE(limpaVerboseUse)
      }

      dpcCall <- c(list(y = dataMat), dpcArgsUse)
      do.call(limpa::dpcCN, dpcCall)
    }

    inferLimpaProteinIdCol <- function(annotationMat, proteinIdCol) {
      if (is.null(proteinIdCol)) {
        return(NULL)
      }

      proteinIdCol <- as.character(proteinIdCol)[1]
      annotationCols <- colnames(annotationMat)
      if (is.null(annotationCols)) {
        annotationCols <- character()
      }

      if (identical(proteinIdCol, "auto")) {
        candidates <- c("Protein.Group", "Protein")
        proteinIdCol <- candidates[candidates %in% annotationCols][1]
        if (is.na(proteinIdCol) || is.null(proteinIdCol)) {
          return(NULL)
        }
        return(proteinIdCol)
      }

      if (!(proteinIdCol %in% annotationCols)) {
        stop(
          "limpaProteinIdCol '",
          proteinIdCol,
          "' was not found in the row annotation.\n",
          "Available columns: ",
          paste(annotationCols, collapse = ", ")
        )
      }

      proteinIdCol
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
        estimateLimpaDpc(dataMat, method = limpaDpcMethod)
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

      inferStableProteinAnnotationCols <- function(
        genesDf,
        proteinId,
        proteinIdCol
      ) {
        colNames <- colnames(genesDf)
        if (is.null(colNames) || length(colNames) == 0) {
          return(character())
        }

        candidates <- setdiff(colNames, proteinIdCol)
        if (length(candidates) == 0) {
          return(character())
        }

        proteinId <- as.character(proteinId)
        idxByProtein <- split(seq_along(proteinId), proteinId)

        isStable <- function(values) {
          values <- as.character(values)
          values[values == ""] <- NA_character_

          if (!any(!is.na(values))) {
            return(FALSE)
          }

          all(vapply(
            idxByProtein,
            function(idx) {
              x <- values[idx]
              x <- x[!is.na(x)]
              length(unique(x)) <= 1
            },
            logical(1)
          ))
        }

        stable <- vapply(
          candidates,
          function(col) isStable(genesDf[[col]]),
          logical(1)
        )
        candidates[stable]
      }

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
      genesInput <- genesInput[, unique(c(proteinIdCol, stableCols)), drop = FALSE]
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

    filterLowRepLimpa <- function(df, groups, leastRep = 1) {
      hasAnyObs <- rowSums(!is.na(df)) > 0
      df <- df[hasAnyObs, , drop = FALSE]
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
        counts <- rowSums(!is.na(df[, cols, drop = FALSE]))
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

    limpaProteinIdColUsed <- if (type == "limpa") {
      inferLimpaProteinIdCol(annotMat(nst), limpaProteinIdCol)
    } else {
      NULL
    }

    if (type == "limpa") {
      warnIfLimpaProteinSummarizationWillBeSkipped(
        annotMat(nst),
        limpaProteinIdColUsed
      )
    }

    limpaQuantified <- if (type == "limpa" && !is.null(limpaProteinIdColUsed)) {
      prepareLimpaQuantified(dataMat, annotMat(nst), limpaProteinIdColUsed)
    } else {
      NULL
    }

    if (!is.null(limpaBackend)) {
      limpaBackend$proteinIdCol <- limpaProteinIdColUsed
      limpaBackend$usedDpcQuant <- !is.null(limpaQuantified)
      if (!is.null(limpaQuantified) && limpaKeep %in% c("elist", "all")) {
        limpaBackend$quantifiedEList <- limpaQuantified
      }
    }

    if (!is.null(limpaQuantified)) {
      dataMat <- as.matrix(limpaQuantified$E)
      slot(nst, "dataMat") <- dataMat
      slot(nst, "annotMat") <- as.matrix(limpaQuantified$genes)
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
        yUse$other$n.observations <- yUse$other$n.observations[, cols, drop = FALSE]
      }
      if (!is.null(yUse$other$standard.error)) {
        yUse$other$standard.error <- yUse$other$standard.error[, cols, drop = FALSE]
      }

      rows <- rownames(dataMat)
      yUse$E <- yUse$E[rows, , drop = FALSE]

      if (!is.null(yUse$other$n.observations)) {
        yUse$other$n.observations <- yUse$other$n.observations[
          rows,
          ,
          drop = FALSE
        ]
      }
      if (!is.null(yUse$other$standard.error)) {
        yUse$other$standard.error <- yUse$other$standard.error[
          rows,
          ,
          drop = FALSE
        ]
      }
      if (!is.null(yUse$genes)) {
        yUse$genes <- yUse$genes[rows, , drop = FALSE]
      }

      limpaDEArgsUse <-
        applyLimpaDEDefaultsAndValidate(sanitizeLimpaDEArgs(limpaDEArgs))
      fit <- do.call(
        limpa::dpcDE,
        c(list(y = yUse, design = limmaDesign, plot = FALSE), limpaDEArgsUse)
      )
      recordLimpaBackend(backendKey, y = yUse, fit = fit, design = limmaDesign)
      fit
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
            leastRep = leastRepCount
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
            leastRep = leastRepCount
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
              leastRep = leastRepCount
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
