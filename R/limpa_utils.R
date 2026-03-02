requireLimpaPackageInternal <- function(context = "type='limpa'") {
  if (!requireNamespace("limpa", quietly = TRUE)) {
    stop(
      context,
      " requires the optional Bioconductor package 'limpa'.\n",
      "Install it with `BiocManager::install(\"limpa\")`."
    )
  }
}

validateLimpaDpc <- function(limpaDpc) {
  if (
    !is.null(limpaDpc) &&
      !(is.list(limpaDpc) || (is.numeric(limpaDpc) && length(limpaDpc) == 2))
  ) {
    stop(
      "limpaDpc must be NULL, a list returned by limpa::dpc(), or a numeric vector c(beta0, beta1)."
    )
  }

  invisible(limpaDpc)
}

sanitizeLimpaDpcArgs <- function(args) {
  if (is.null(args) || length(args) == 0) {
    return(list())
  }
  if (!is.list(args)) {
    stop("limpaDpcArgs must be a named list (or NULL).")
  }
  if (is.null(names(args))) {
    stop("limpaDpcArgs must be a named list.")
  }

  args[["y"]] <- NULL
  args
}

sanitizeLimpaQuantArgs <- function(args) {
  if (is.null(args) || length(args) == 0) {
    return(list())
  }
  if (!is.list(args)) {
    stop("limpaQuantArgs must be a named list (or NULL).")
  }
  if (is.null(names(args))) {
    stop("limpaQuantArgs must be a named list.")
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
    stop("limpaDEArgs must be a list (or NULL).")
  }
  if (is.null(names(args))) {
    stop("limpaDEArgs must be a named list.")
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
    stop(
      argsLabel,
      " contains unsupported argument names for limpaDpcMethod='",
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

normalizeLimpaByRowConfig <- function(limpaByRow, limpaProteinIdCol) {
  limpaByRowUse <- as.logical(limpaByRow)[1]
  if (is.na(limpaByRowUse)) {
    stop("limpaByRow must be TRUE or FALSE.")
  }

  if (isTRUE(limpaByRowUse)) {
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
  limpaDpcMethod = c("none", "dpc", "dpcCN"),
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
    stop(
      "Multiple '*_limpa_quantified.rds' files were found in ",
      dataDir,
      ", and no canonical file named '",
      basename(canonical),
      "' was present. Please specify `limpaQuantifiedRds` explicitly."
    )
  }

  warning(
    "Found a single '*_limpa_quantified.rds' file in ",
    dataDir,
    " (",
    basename(candidates[[1]]),
    "), but it does not match the expected canonical filename '",
    basename(canonical),
    "'. Auto-detection was skipped; please set `limpaQuantifiedRds` explicitly.",
    call. = FALSE
  )
  NULL
}
