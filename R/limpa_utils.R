requireLimpaPackageInternal <- function(context = "type='limpa'") {
  if (!requireNamespace("limpa", quietly = TRUE)) {
    cli::cli_abort(
      c(
        "{context} requires the optional Bioconductor package {.pkg limpa}.",
        i = "Install it with {.code BiocManager::install('limpa')}."
      ),
      class = "normalyzerde_error",
      call = NULL
    )
  }
}

validateLimpaDpc <- function(limpaDpc) {
  if (
    !is.null(limpaDpc) &&
      !(is.list(limpaDpc) || (is.numeric(limpaDpc) && length(limpaDpc) == 2))
  ) {
    cli::cli_abort(
      "{.arg limpaDpc} must be NULL, a list returned by limpa::dpc()/dpcON()/dpcCN(), or a numeric vector c(beta0, beta1).",
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
      "{.arg limpaDpcArgs} must be a named list (or NULL).",
      class = "normalyzerde_error",
      call = NULL
    )
  }
  if (is.null(names(args))) {
    cli::cli_abort(
      "{.arg limpaDpcArgs} must be a named list.",
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
      "{.arg limpaQuantArgs} must be a named list (or NULL).",
      class = "normalyzerde_error",
      call = NULL
    )
  }
  if (is.null(names(args))) {
    cli::cli_abort(
      "{.arg limpaQuantArgs} must be a named list.",
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
      "{.arg limpaDEArgs} must be a list (or NULL).",
      class = "normalyzerde_error",
      call = NULL
    )
  }
  if (is.null(names(args))) {
    cli::cli_abort(
      "{.arg limpaDEArgs} must be a named list.",
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
      "{.arg limpaQuantArgs$dpc.slope} must be a single positive numeric value.",
      class = "normalyzerde_error",
      call = NULL
    )
  }
  args[["dpc.slope"]] <- slope

  chunk <- as.integer(args[["chunk"]])[1]
  if (is.na(chunk) || chunk < 1) {
    cli::cli_abort(
      "{.arg limpaQuantArgs$chunk} must be a positive integer.",
      class = "normalyzerde_error",
      call = NULL
    )
  }
  args[["chunk"]] <- chunk

  verbose <- as.logical(args[["verbose"]])[1]
  if (is.na(verbose)) {
    cli::cli_abort(
      "{.arg limpaQuantArgs$verbose} must be TRUE or FALSE.",
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
      "{.arg limpaDEArgs$sample.weights} must be TRUE or FALSE.",
      class = "normalyzerde_error",
      call = NULL
    )
  }
  args[["sample.weights"]] <- sampleWeights
  args
}

validateLimpaArgsByFormals <- function(
  args,
  fn,
  methodLabel,
  argsLabel = "limpaDpcArgs"
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
        "{.arg {argsLabel}} contains unsupported argument names for {.arg limpaDpcMethod}={.val {methodLabel}}.",
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
      "{.arg limpaByRow} must be TRUE or FALSE.",
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
        "{.arg limpaByRow}=TRUE is incompatible with a non-default {.arg limpaProteinIdCol}. Set {.arg limpaProteinIdCol}=NULL (or leave it as {.val auto}).",
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
        "{.arg limpaProteinIdCol} {.val {proteinIdCol}} was not found in the row annotation.",
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
      argsLabel = "limpaDpcArgs"
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
          "{.arg limpaDpcMethod}={.val dpcON} requires limpa::dpcON(), but it was not found in the installed limpa version.",
          i = "Update limpa or use {.arg limpaDpcMethod}={.val dpc}."
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
      argsLabel = "limpaDpcArgs"
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
    argsLabel = "limpaDpcArgs"
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

autoDetectLimpaQuantifiedRds <- function(dataPath) {
  if (
    is.null(dataPath) ||
      !is.character(dataPath) ||
      length(dataPath) != 1 ||
      is.na(dataPath) ||
      !nzchar(dataPath)
  ) {
    return(NULL)
  }

  dataDir <- dirname(dataPath)
  if (!dir.exists(dataDir)) {
    return(NULL)
  }

  candidates <- list.files(
    dataDir,
    pattern = "_limpa_quantified\\.rds$",
    full.names = TRUE
  )
  if (length(candidates) == 0) {
    return(NULL)
  }

  canonical <- file.path(
    dataDir,
    paste0(basename(dataDir), "_limpa_quantified.rds")
  )
  if (canonical %in% candidates) {
    return(canonical)
  }

  if (length(candidates) > 1) {
    cli::cli_abort(
      c(
        "Multiple {.path '*_limpa_quantified.rds'} files were found in {.path {dataDir}}, and no canonical file named {.path {basename(canonical)}} was present.",
        i = "Please specify {.arg limpaQuantifiedRds} explicitly."
      ),
      class = "normalyzerde_error",
      call = NULL
    )
  }

  cli::cli_warn(
    c(
      "Found a single {.path '*_limpa_quantified.rds'} file in {.path {dataDir}} ({.path {basename(candidates[[1]])}}), but it does not match the expected canonical filename {.path {basename(canonical)}}.",
      i = "Auto-detection was skipped; please set {.arg limpaQuantifiedRds} explicitly."
    ),
    class = "normalyzerde_warning",
    call = NULL
  )
  NULL
}
