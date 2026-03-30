requireLimpaPackageInternal <- function(context = "type='limpa'") {
  if (!requireNamespace("limpa", quietly = TRUE)) {
    cli::cli_abort(
      c(
        "NormalyzerDE requires the Bioconductor package {.pkg limpa}, but it is not available for {context}.",
        i = "Reinstall NormalyzerDE and its Bioconductor dependencies, including {.pkg limpa}."
      ),
      class = "normalyzerde_error",
      call = NULL
    )
  }
}

#' Create validated limpa options
#'
#' Helper to construct a grouped options object for `preQuant = "limpa"` and
#' `type = "limpa"` workflows.
#'
#' @param proteinIdCol Optional protein identifier column used for summarizing
#'   peptide/precursor rows to proteins. Use `"auto"` (default) to infer a
#'   common identifier, or `NULL` to keep each row separate.
#' @param byRow Whether to quantify each row separately with
#'   `limpa::dpcQuantByRow()`.
#' @param dpc Optional DPC parameters to pass to `limpa`.
#' @param dpcMethod DPC estimation method when `dpc` is not supplied.
#' @param dpcArgs Optional named list of extra arguments for DPC estimation.
#' @param quantArgs Optional named list of extra arguments for
#'   `limpa::dpcQuant()` / `limpa::dpcQuantByRow()`.
#' @param quantifiedRds Optional path to a quantified `EList` RDS file for
#'   reuse with `normalyzerDE(type = "limpa")`.
#' @param deArgs Optional named list of extra arguments for `limpa::dpcDE()`.
#' @param keep Which intermediate `limpa` objects to retain in
#'   `backendData(nst)$limpa`.
#' @param postQuantNorm Optional between-sample normalization applied after
#'   quantification and before `limpa::dpcDE()`.
#' @return A validated list suitable to pass as `limpaOptions`.
#' @export
#' @examples
#' limpaOptions(
#'   byRow = TRUE,
#'   quantArgs = list(chunk = 1000L, verbose = FALSE)
#' )
limpaOptions <- function(
  proteinIdCol = "auto",
  byRow = FALSE,
  dpc = NULL,
  dpcMethod = c("none", "dpc", "dpcON", "dpcCN"),
  dpcArgs = NULL,
  quantArgs = NULL,
  quantifiedRds = NULL,
  deArgs = NULL,
  keep = c("none", "elist", "fit", "all"),
  postQuantNorm = c(
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
  call <- sys.call()
  argNames <- names(call)[-1]
  argNames <- argNames[nzchar(argNames)]
  unknown <- setdiff(argNames, names(formals(sys.function())))
  if (length(unknown) > 0) {
    cli::cli_abort(
      "Unknown argument(s): {paste(unknown, collapse = ', ')}.",
      class = "normalyzerde_error",
      call = NULL
    )
  }

  buildLimpaOptionsInternal(
    proteinIdCol = proteinIdCol,
    byRow = byRow,
    dpc = dpc,
    dpcMethod = dpcMethod,
    dpcArgs = dpcArgs,
    quantArgs = quantArgs,
    quantifiedRds = quantifiedRds,
    deArgs = deArgs,
    keep = keep,
    postQuantNorm = postQuantNorm
  )
}

buildLimpaOptionsInternal <- function(
  proteinIdCol = "auto",
  byRow = FALSE,
  dpc = NULL,
  dpcMethod = c("none", "dpc", "dpcON", "dpcCN"),
  dpcArgs = NULL,
  quantArgs = NULL,
  quantifiedRds = NULL,
  deArgs = NULL,
  keep = c("none", "elist", "fit", "all"),
  postQuantNorm = c(
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
  if (!is.null(proteinIdCol)) {
    proteinIdCol <- as.character(proteinIdCol)
    if (
      length(proteinIdCol) != 1 ||
        is.na(proteinIdCol) ||
        !nzchar(proteinIdCol)
    ) {
      cli::cli_abort(
        "{.arg proteinIdCol} must be a single non-empty character value (or {.val NULL}).",
        class = "normalyzerde_error",
        call = NULL
      )
    }
  }

  byRow <- as.logical(byRow)[1]
  if (is.na(byRow)) {
    cli::cli_abort(
      "{.arg byRow} must be TRUE or FALSE.",
      class = "normalyzerde_error",
      call = NULL
    )
  }

  validateLimpaDpc(dpc)
  dpcMethod <- match.arg(dpcMethod)
  keep <- match.arg(keep)
  postQuantNorm <- match.arg(postQuantNorm)
  dpcArgs <- sanitizeLimpaDpcArgs(dpcArgs)
  quantArgs <- sanitizeLimpaQuantArgs(quantArgs)
  deArgs <- sanitizeLimpaDEArgs(deArgs)

  if (!is.null(quantifiedRds)) {
    quantifiedRds <- as.character(quantifiedRds)
    if (
      length(quantifiedRds) != 1 ||
        is.na(quantifiedRds) ||
        !nzchar(quantifiedRds)
    ) {
      cli::cli_abort(
        "{.arg quantifiedRds} must be a single non-empty file path (or {.val NULL}).",
        class = "normalyzerde_error",
        call = NULL
      )
    }
  }

  structure(
    list(
      proteinIdCol = proteinIdCol,
      byRow = byRow,
      dpc = dpc,
      dpcMethod = dpcMethod,
      dpcArgs = dpcArgs,
      quantArgs = quantArgs,
      quantifiedRds = quantifiedRds,
      deArgs = deArgs,
      keep = keep,
      postQuantNorm = postQuantNorm
    ),
    class = c("normalyzerde_limpa_options", "list")
  )
}

normalizeLimpaOptionsInput <- function(limpaOptions) {
  if (is.null(limpaOptions)) {
    return(NULL)
  }
  if (!is.list(limpaOptions)) {
    cli::cli_abort(
      "{.arg limpaOptions} must be a list created by {.fn limpaOptions} (or {.val NULL}).",
      class = "normalyzerde_error",
      call = NULL
    )
  }
  if (length(limpaOptions) == 0) {
    return(buildLimpaOptionsInternal())
  }
  if (is.null(names(limpaOptions)) || anyNA(names(limpaOptions)) || any(names(limpaOptions) == "")) {
    cli::cli_abort(
      "{.arg limpaOptions} must be a named list.",
      class = "normalyzerde_error",
      call = NULL
    )
  }

  allowed <- names(formals(buildLimpaOptionsInternal))
  unknown <- setdiff(names(limpaOptions), allowed)
  if (length(unknown) > 0) {
    cli::cli_abort(
      c(
        "{.arg limpaOptions} contains unknown field(s): {paste(unknown, collapse = ', ')}.",
        i = "Create it with {.fn limpaOptions} to validate names automatically."
      ),
      class = "normalyzerde_error",
      call = NULL
    )
  }

  do.call(buildLimpaOptionsInternal, limpaOptions)
}

resolveLimpaOptions <- function(limpaOptions = NULL) {
  limpaOptions <- normalizeLimpaOptionsInput(limpaOptions)
  if (is.null(limpaOptions)) {
    return(buildLimpaOptionsInternal())
  }
  limpaOptions
}

resolveLimpaQuantByRowFn <- function(limpaNs = NULL) {
  if (is.null(limpaNs)) {
    requireLimpaPackageInternal()
    limpaNs <- asNamespace("limpa")
  }

  if (exists("dpcQuantByRow", envir = limpaNs, inherits = FALSE)) {
    return(get("dpcQuantByRow", envir = limpaNs, inherits = FALSE))
  }
  if (exists("dpcImpute", envir = limpaNs, inherits = FALSE)) {
    return(get("dpcImpute", envir = limpaNs, inherits = FALSE))
  }

  cli::cli_abort(
    c(
      "{.pkg limpa} does not export a row-wise quantification function compatible with NormalyzerDE.",
      i = "Expected {.code dpcQuantByRow()} or the older {.code dpcImpute()}."
    ),
    class = "normalyzerde_error",
    call = NULL
  )
}

validateLimpaDpc <- function(limpaDpc) {
  if (
    !is.null(limpaDpc) &&
      !(is.list(limpaDpc) || (is.numeric(limpaDpc) && length(limpaDpc) == 2))
  ) {
    cli::cli_abort(
      "{.arg dpc} must be NULL, a list returned by limpa::dpc()/dpcON()/dpcCN(), or a numeric vector c(beta0, beta1).",
      class = "normalyzerde_error",
      call = NULL
    )
  }

  invisible(limpaDpc)
}

sanitizeLimpaDpcArgs <- function(args) {
  if (is.null(args) || length(args) == 0) {
    return(list())
  }
  if (!is.list(args)) {
    cli::cli_abort(
      "{.arg dpcArgs} must be a named list (or NULL).",
      class = "normalyzerde_error",
      call = NULL
    )
  }
  if (is.null(names(args))) {
    cli::cli_abort(
      "{.arg dpcArgs} must be a named list.",
      class = "normalyzerde_error",
      call = NULL
    )
  }

  args[["y"]] <- NULL
  args
}

sanitizeLimpaQuantArgs <- function(args) {
  if (is.null(args) || length(args) == 0) {
    return(list())
  }
  if (!is.list(args)) {
    cli::cli_abort(
      "{.arg quantArgs} must be a named list (or NULL).",
      class = "normalyzerde_error",
      call = NULL
    )
  }
  if (is.null(names(args))) {
    cli::cli_abort(
      "{.arg quantArgs} must be a named list.",
      class = "normalyzerde_error",
      call = NULL
    )
  }

  forbidden <- c("y", "protein.id", "dpc")
  args[forbidden] <- NULL
  args
}

sanitizeLimpaDEArgs <- function(args) {
  if (is.null(args) || length(args) == 0) {
    return(list())
  }
  if (!is.list(args)) {
    cli::cli_abort(
      "{.arg deArgs} must be a list (or NULL).",
      class = "normalyzerde_error",
      call = NULL
    )
  }
  if (is.null(names(args))) {
    cli::cli_abort(
      "{.arg deArgs} must be a named list.",
      class = "normalyzerde_error",
      call = NULL
    )
  }

  forbidden <- c("y", "design", "plot")
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
    cli::cli_abort(
      "{.arg quantArgs$dpc.slope} must be a single positive numeric value.",
      class = "normalyzerde_error",
      call = NULL
    )
  }
  args[["dpc.slope"]] <- slope

  chunk <- as.integer(args[["chunk"]])[1]
  if (is.na(chunk) || chunk < 1) {
    cli::cli_abort(
      "{.arg quantArgs$chunk} must be a positive integer.",
      class = "normalyzerde_error",
      call = NULL
    )
  }
  args[["chunk"]] <- chunk

  verbose <- as.logical(args[["verbose"]])[1]
  if (is.na(verbose)) {
    cli::cli_abort(
      "{.arg quantArgs$verbose} must be TRUE or FALSE.",
      class = "normalyzerde_error",
      call = NULL
    )
  }
  args[["verbose"]] <- verbose

  args
}

applyLimpaDEDefaultsAndValidate <- function(args) {
  if (is.null(args[["sample.weights"]])) {
    args[["sample.weights"]] <- FALSE
  }
  sampleWeights <- as.logical(args[["sample.weights"]])[1]
  if (is.na(sampleWeights)) {
    cli::cli_abort(
      "{.arg deArgs$sample.weights} must be TRUE or FALSE.",
      class = "normalyzerde_error",
      call = NULL
    )
  }
  args[["sample.weights"]] <- sampleWeights
  args
}

normalizeLimpaQuantifiedEList <- function(
  yQuant,
  proteinIdCol = NULL,
  proteinIds = NULL
) {
  rowIds <- as.character(seq_len(nrow(yQuant$E)))
  rownames(yQuant$E) <- rowIds

  if (!is.null(yQuant$other$n.observations)) {
    rownames(yQuant$other$n.observations) <- rowIds
  }
  if (!is.null(yQuant$other$standard.error)) {
    rownames(yQuant$other$standard.error) <- rowIds
  }

  if (!is.null(yQuant$genes)) {
    genes <- as.data.frame(yQuant$genes, check.names = FALSE)
  } else {
    genes <- data.frame(check.names = FALSE)
  }

  if (!is.null(proteinIdCol)) {
    if (!(proteinIdCol %in% colnames(genes)) && !is.null(proteinIds)) {
      genes[[proteinIdCol]] <- proteinIds
    }
    if (proteinIdCol %in% colnames(genes)) {
      genes <- genes[,
        c(proteinIdCol, setdiff(names(genes), proteinIdCol)),
        drop = FALSE
      ]
    }
  }

  if (ncol(genes) > 0) {
    rownames(genes) <- rowIds
    yQuant$genes <- genes
  }

  yQuant
}

quantifyLimpaByRowInternal <- function(
  dataMat,
  genesDf,
  dpc,
  quantArgs,
  limpaQuantByRow = NULL
) {
  if (is.null(limpaQuantByRow)) {
    limpaQuantByRow <- resolveLimpaQuantByRowFn()
  }

  yInput <- methods::new("EList", list(E = dataMat, genes = genesDf))
  yQuant <- do.call(
    limpaQuantByRow,
    c(list(y = yInput, dpc = dpc), quantArgs)
  )

  normalizeLimpaQuantifiedEList(yQuant)
}

quantifyLimpaByProteinInternal <- function(
  dataMat,
  genesDf,
  proteinIdCol,
  dpc,
  quantArgs
) {
  proteinId <- as.character(genesDf[[proteinIdCol]])

  if (length(proteinId) != nrow(dataMat)) {
    cli::cli_abort(
      "Row annotation column {.val {proteinIdCol}} does not match the number of rows in the data matrix.",
      class = "normalyzerde_error",
      call = NULL
    )
  }

  if (anyNA(proteinId) || any(proteinId == "")) {
    cli::cli_abort(
      c(
        "Row annotation column {.val {proteinIdCol}} contains missing or empty protein identifiers.",
        i = "Remove these rows or choose another column."
      ),
      class = "normalyzerde_error",
      call = NULL
    )
  }

  stableCols <- inferStableProteinAnnotationCols(
    genesDf = genesDf,
    proteinId = proteinId,
    proteinIdCol = proteinIdCol
  )
  genesForQuant <- genesDf[,
    unique(c(proteinIdCol, stableCols)),
    drop = FALSE
  ]

  yPeptide <- methods::new("EList", list(E = dataMat, genes = genesForQuant))
  yProtein <- do.call(
    limpa::dpcQuant,
    c(
      list(
        y = yPeptide,
        protein.id = proteinIdCol,
        dpc = dpc
      ),
      quantArgs
    )
  )

  normalizeLimpaQuantifiedEList(
    yProtein,
    proteinIdCol = proteinIdCol,
    proteinIds = rownames(yProtein$E)
  )
}

extractLimpaSampleWeightsFromFit <- function(fit) {
  if (
    is.null(fit) || is.null(fit$targets) || is.null(fit$targets$sample.weight)
  ) {
    return(NULL)
  }

  sampleWeights <- as.numeric(fit$targets$sample.weight)
  sampleNames <- NULL

  if (!is.null(fit$EList) && !is.null(fit$EList$E)) {
    sampleNames <- colnames(fit$EList$E)
  }
  if (is.null(sampleNames)) {
    sampleNames <- names(fit$targets$sample.weight)
  }
  if (
    is.null(sampleNames) &&
      !is.null(rownames(fit$targets)) &&
      length(rownames(fit$targets)) == length(sampleWeights)
  ) {
    sampleNames <- rownames(fit$targets)
  }

  if (!is.null(sampleNames) && length(sampleNames) == length(sampleWeights)) {
    names(sampleWeights) <- as.character(sampleNames)
  }

  sampleWeights
}

collectLimpaSampleWeights <- function(object, comparison = NULL) {
  backend <- backendData(object)[["limpa"]]
  if (is.null(backend)) {
    return(NULL)
  }

  sampleWeights <- backend$sampleWeights
  if (is.null(sampleWeights) || length(sampleWeights) == 0) {
    fits <- backend$fits
    if (!is.null(fits) && length(fits) > 0) {
      sampleWeights <- lapply(fits, extractLimpaSampleWeightsFromFit)
      sampleWeights <- sampleWeights[
        !vapply(
          sampleWeights,
          is.null,
          logical(1)
        )
      ]
    }
  }

  if (is.null(sampleWeights) || length(sampleWeights) == 0) {
    return(NULL)
  }

  keys <- names(sampleWeights)
  if (is.null(keys)) {
    keys <- rep.int(".global", length(sampleWeights))
  }
  keep <- !vapply(sampleWeights, function(x) length(x) == 0, logical(1))
  sampleWeights <- sampleWeights[keep]
  keys <- keys[keep]

  if (!is.null(comparison)) {
    comparison <- as.character(comparison)
    keep <- keys %in% comparison
    sampleWeights <- sampleWeights[keep]
    keys <- keys[keep]
  }

  if (length(sampleWeights) == 0) {
    return(data.frame(
      comparison = character(),
      sample = character(),
      sampleWeight = numeric(),
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  }

  rows <- lapply(seq_along(sampleWeights), function(i) {
    weights <- sampleWeights[[i]]
    sampleNames <- names(weights)
    if (is.null(sampleNames)) {
      sampleNames <- as.character(seq_along(weights))
    }

    data.frame(
      comparison = rep(keys[[i]], length(weights)),
      sample = as.character(sampleNames),
      sampleWeight = as.numeric(weights),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })

  do.call(rbind, rows)
}

#' Extract estimated limpa sample weights
#'
#' Returns the sample-specific quality weights estimated by
#' \code{limpa::dpcDE()} / \code{limpa::voomaLmFitWithImputation()} when
#' \code{limpaOptions(deArgs = list(sample.weights = TRUE))} was used.
#'
#' @param object A \code{NormalyzerStatistics} object returned by
#'   \code{\link{calculateContrasts}} with \code{type = "limpa"}.
#' @param comparison Optional comparison/backend key to extract. If \code{NULL}
#'   (default), weights for all available keys are returned.
#' @return A data frame with columns \code{comparison}, \code{sample}, and
#'   \code{sampleWeight}.
#' @export
getLimpaSampleWeights <- function(object, comparison = NULL) {
  if (!methods::is(object, "NormalyzerStatistics")) {
    cli::cli_abort(
      "{.arg object} must be a {.cls NormalyzerStatistics} object.",
      class = "normalyzerde_error",
      call = NULL
    )
  }

  out <- collectLimpaSampleWeights(object, comparison = comparison)
  if (is.null(out) || nrow(out) == 0) {
    cli::cli_abort(
      c(
        "No limpa sample weights were found in {.arg object}.",
        i = "Run {.fn calculateContrasts} with {.arg type}={.val limpa} and {.code limpaOptions = limpaOptions(deArgs = list(sample.weights = TRUE))}."
      ),
      class = "normalyzerde_error",
      call = NULL
    )
  }

  out
}

validateLimpaArgsByFormals <- function(
  args,
  fn,
  methodLabel,
  argsLabel = "dpcArgs"
) {
  if (length(args) == 0) {
    return(args)
  }

  allowed <- names(formals(fn))
  allowed <- allowed[!is.na(allowed) & nzchar(allowed)]
  allowed <- setdiff(allowed, "y")

  unknown <- setdiff(names(args), allowed)
  if (length(unknown) > 0) {
    cli::cli_abort(
      c(
        "{.arg {argsLabel}} contains unsupported argument names for {.arg dpcMethod}={.val {methodLabel}}.",
        i = "Unsupported: {paste(unknown, collapse = \", \")}.",
        i = "Allowed: {paste(allowed, collapse = \", \")}."
      ),
      class = "normalyzerde_error",
      call = NULL
    )
  }

  args
}

normalizeLimpaByRowConfig <- function(limpaByRow, limpaProteinIdCol) {
  limpaByRowUse <- as.logical(limpaByRow)[1]
  if (is.na(limpaByRowUse)) {
    cli::cli_abort(
      "{.arg byRow} must be TRUE or FALSE.",
      class = "normalyzerde_error",
      call = NULL
    )
  }

  if (isTRUE(limpaByRowUse)) {
    proteinIdColValue <- if (is.null(limpaProteinIdCol)) {
      NULL
    } else {
      as.character(limpaProteinIdCol)[1]
    }
    if (!is.null(proteinIdColValue) && !identical(proteinIdColValue, "auto")) {
      cli::cli_abort(
        "{.arg byRow}=TRUE is incompatible with a non-default {.arg proteinIdCol}. Set {.arg proteinIdCol}=NULL (or leave it as {.val auto}).",
        class = "normalyzerde_error",
        call = NULL
      )
    }
    limpaProteinIdCol <- NULL
  }

  list(
    limpaByRow = limpaByRowUse,
    limpaProteinIdCol = limpaProteinIdCol
  )
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
    cli::cli_abort(
      c(
        "{.arg proteinIdCol} {.val {proteinIdCol}} was not found in the row annotation.",
        i = "Available columns: {paste(annotationCols, collapse = \", \")}."
      ),
      class = "normalyzerde_error",
      call = NULL
    )
  }

  proteinIdCol
}

inferStableProteinAnnotationCols <- function(genesDf, proteinId, proteinIdCol) {
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

estimateLimpaDpcFromData <- function(
  dataMat,
  limpaDpcMethod = c("none", "dpc", "dpcON", "dpcCN"),
  limpaDpcArgs = list(),
  dpcSlope = 0.8,
  verbose = FALSE
) {
  method <- match.arg(limpaDpcMethod)
  if (identical(method, "none")) {
    return(NULL)
  }

  if (identical(method, "dpc")) {
    dpcArgsUse <- validateLimpaArgsByFormals(
      limpaDpcArgs,
      fn = limpa::dpc,
      methodLabel = method,
      argsLabel = "dpcArgs"
    )
    dpcCall <- c(list(y = dataMat), dpcArgsUse)
    if (!isTRUE(verbose)) {
      return(suppressMessages(do.call(limpa::dpc, dpcCall)))
    }
    return(do.call(limpa::dpc, dpcCall))
  }

  if (identical(method, "dpcON")) {
    limpaNs <- asNamespace("limpa")
    if (!exists("dpcON", envir = limpaNs, inherits = FALSE)) {
      cli::cli_abort(
        c(
          "{.arg dpcMethod}={.val dpcON} requires limpa::dpcON(), but it was not found in the installed limpa version.",
          i = "Update limpa or use {.arg dpcMethod}={.val dpc}."
        ),
        class = "normalyzerde_error",
        call = NULL
      )
    }
    dpcON <- get("dpcON", envir = limpaNs, inherits = FALSE)

    dpcArgsUse <- validateLimpaArgsByFormals(
      limpaDpcArgs,
      fn = dpcON,
      methodLabel = method,
      argsLabel = "dpcArgs"
    )
    if (!("dpc.slope.start" %in% names(dpcArgsUse))) {
      dpcArgsUse[["dpc.slope.start"]] <- dpcSlope
    }
    if (!("verbose" %in% names(dpcArgsUse))) {
      dpcArgsUse[["verbose"]] <- isTRUE(verbose)
    }

    dpcCall <- c(list(y = dataMat), dpcArgsUse)
    return(do.call(dpcON, dpcCall))
  }

  dpcArgsUse <- validateLimpaArgsByFormals(
    limpaDpcArgs,
    fn = limpa::dpcCN,
    methodLabel = method,
    argsLabel = "dpcArgs"
  )
  if (!("dpc.slope.start" %in% names(dpcArgsUse))) {
    dpcArgsUse[["dpc.slope.start"]] <- dpcSlope
  }
  if (!("verbose" %in% names(dpcArgsUse))) {
    dpcArgsUse[["verbose"]] <- isTRUE(verbose)
  }

  dpcCall <- c(list(y = dataMat), dpcArgsUse)
  do.call(limpa::dpcCN, dpcCall)
}
