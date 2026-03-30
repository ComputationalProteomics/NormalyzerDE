diannIsParquet <- function(filePath) {
  grepl("\\.parquet$", filePath, ignore.case = TRUE)
}

diannReadHeader <- function(filePath, sep = "\t") {
  headerLine <- readLines(filePath, n = 1, warn = FALSE)
  if (length(headerLine) < 1) {
    stop("DIANN input file was empty: ", filePath)
  }
  strsplit(headerLine, sep, fixed = TRUE)[[1]]
}

diannIsReportHeader <- function(header) {
  any(c("Run", "File.Name") %in% header)
}

diannSelectFirstPresent <- function(candidates, available) {
  first <- candidates[candidates %in% available][1]
  if (is.na(first) || is.null(first)) {
    NULL
  } else {
    first
  }
}

diannResolveMinPositive <- function(diannMinPositive) {
  if (is.null(diannMinPositive)) {
    return(0)
  }

  diannMinPositive <- as.numeric(diannMinPositive)[1]
  if (is.na(diannMinPositive) || diannMinPositive < 0) {
    stop("diannMinPositive must be a single non-negative numeric value.")
  }

  diannMinPositive
}

diannNormalizeInputOptions <- function(inputOptions, sep = "\t") {
  if (is.null(inputOptions)) {
    inputOptions <- list()
  }
  if (!is.list(inputOptions)) {
    stop("inputOptions must be a list (or NULL).")
  }

  level <- inputOptions$level
  if (is.null(level) || length(level) == 0) {
    level <- "auto"
  }
  level <- match.arg(as.character(level)[1], c("auto", "protein", "precursor"))

  if (!is.null(inputOptions$sep)) {
    sep <- as.character(inputOptions$sep)[1]
  }

  columns <- inputOptions$columns
  if (is.null(columns)) {
    columns <- list()
  }
  if (!is.list(columns)) {
    stop("inputOptions$columns must be a list (or NULL).")
  }

  filters <- inputOptions$filters
  if (is.null(filters)) {
    filters <- list()
  }
  if (!is.list(filters)) {
    stop("inputOptions$filters must be a list (or NULL).")
  }

  filterDecoy <- filters$decoy
  if (is.null(filterDecoy)) {
    filterDecoy <- TRUE
  }
  filterDecoy <- isTRUE(filterDecoy)

  q <- filters$q
  filterQValue <- TRUE
  qCols <- NULL
  qCutoffs <- 0.01
  if (is.null(q)) {
    filterQValue <- TRUE
  } else if (is.logical(q)) {
    filterQValue <- isTRUE(q)
  } else if (is.list(q)) {
    if (!is.null(q$enable)) {
      filterQValue <- isTRUE(q$enable)
    }
    if (!is.null(q$cols)) {
      qCols <- q$cols
    }
    if (!is.null(q$cutoffs)) {
      qCutoffs <- q$cutoffs
    }
  } else {
    stop("inputOptions$filters$q must be logical, a list, or NULL.")
  }

  minPositive <- filters$min_positive
  if (is.null(minPositive)) {
    minPositive <- filters$minPositive
  }
  if (is.null(minPositive)) {
    minPositive <- 0
  }
  minPositive <- diannResolveMinPositive(minPositive)

  rt <- inputOptions$rt
  rtCol <- NULL
  if (is.null(rt)) {
    rtCol <- columns$rt
    if (is.null(rtCol)) {
      rtCol <- "RT"
    }
  } else if (is.logical(rt)) {
    if (isTRUE(rt)) {
      rtCol <- columns$rt
      if (is.null(rtCol)) {
        rtCol <- "RT"
      }
    } else {
      rtCol <- NULL
    }
  } else if (is.list(rt)) {
    rtCol <- rt$col
    if (is.null(rtCol)) {
      rtCol <- columns$rt
    }
    if (is.null(rtCol)) {
      rtCol <- "RT"
    }
  } else if (is.character(rt)) {
    rtCol <- as.character(rt)[1]
  } else {
    stop("inputOptions$rt must be a list, logical, character, or NULL.")
  }

  list(
    sep = sep,
    level = level,
    sampleCol = columns$sample,
    featureCol = columns$feature,
    quantityCol = columns$quantity,
    extraCols = columns$extra,
    filterDecoy = filterDecoy,
    filterQValue = filterQValue,
    qCols = qCols,
    qCutoffs = qCutoffs,
    minPositive = minPositive,
    rtCol = rtCol
  )
}

#' Create validated DIA-NN input options
#'
#' Helper to construct an `inputOptions` list for `inputFormat = "diann"`. This
#' keeps DIA-NN-specific parameters out of the main entrypoint signatures while
#' still providing a typo-safe, validated API.
#'
#' @param level One of `"auto"`, `"protein"`, or `"precursor"`. If `"auto"`,
#'   NormalyzerDE tries to infer whether the DIA-NN file is protein- or
#'   precursor-level.
#' @param sep Field separator for DIA-NN TSV files.
#' @param sampleCol Optional DIA-NN report sample column name. If `NULL`,
#'   NormalyzerDE uses `"Run"` or `"File.Name"` when present.
#' @param featureCol Optional DIA-NN report feature column name.
#' @param quantityCol Optional DIA-NN report quantity column name.
#' @param extraCols Optional character vector of additional DIA-NN columns to
#'   carry along as row annotation.
#' @param decoy Whether to filter out decoys when the DIA-NN report includes a
#'   `"Decoy"` column.
#' @param qEnable Whether to filter rows by q-values when available.
#' @param qCols Character vector of q-value columns to use. Use `NULL` or
#'   `"auto"` to select reasonable defaults, or `"none"` to disable q-value
#'   filtering.
#' @param qCutoffs Numeric vector of q-value thresholds. If a single value is
#'   provided it is recycled to match `qCols`.
#' @param minPositive Minimum positive value threshold. Non-zero values below
#'   this are converted to `NA` before analysis.
#' @param rt Retention time column configuration. Use `NULL` (default behavior),
#'   `FALSE` to disable, `TRUE` to enable default `"RT"`, a character column
#'   name, or `list(col = "<name>")`.
#' @return A list suitable to pass as `inputOptions`.
#' @export
#' @examples
#' diannInputOptions(
#'   level = "precursor",
#'   quantityCol = "Precursor.Quantity",
#'   qCols = "Q.Value",
#'   qCutoffs = 0.01,
#'   rt = "RT"
#' )
diannInputOptions <- function(
  level = c("auto", "protein", "precursor"),
  sep = "\t",
  sampleCol = NULL,
  featureCol = NULL,
  quantityCol = NULL,
  extraCols = NULL,
  decoy = TRUE,
  qEnable = TRUE,
  qCols = NULL,
  qCutoffs = 0.01,
  minPositive = 0,
  rt = NULL
) {
  call <- sys.call()
  argNames <- names(call)[-1]
  argNames <- argNames[nzchar(argNames)]
  unknown <- setdiff(argNames, names(formals(sys.function())))
  if (length(unknown) > 0) {
    stop("Unknown argument(s): ", paste(unknown, collapse = ", "))
  }

  level <- match.arg(level)
  if (!is.character(sep) || length(sep) != 1 || is.na(sep) || sep == "") {
    stop("sep must be a single non-empty character value.")
  }

  validateColName <- function(x, argName) {
    if (is.null(x)) {
      return(NULL)
    }
    x <- as.character(x)
    if (length(x) != 1 || is.na(x) || x == "") {
      stop(argName, " must be a single non-empty character value (or NULL).")
    }
    x
  }

  sampleCol <- validateColName(sampleCol, "sampleCol")
  featureCol <- validateColName(featureCol, "featureCol")
  quantityCol <- validateColName(quantityCol, "quantityCol")

  if (!is.null(extraCols)) {
    extraCols <- as.character(extraCols)
    if (anyNA(extraCols) || any(extraCols == "")) {
      stop(
        "extraCols must be a character vector of non-empty column names (or NULL)."
      )
    }
  }

  if (!is.logical(decoy) || length(decoy) != 1 || is.na(decoy)) {
    stop("decoy must be TRUE or FALSE.")
  }

  if (!is.logical(qEnable) || length(qEnable) != 1 || is.na(qEnable)) {
    stop("qEnable must be TRUE or FALSE.")
  }
  if (!is.null(qCols)) {
    qCols <- as.character(qCols)
    if (anyNA(qCols) || any(qCols == "")) {
      stop(
        "qCols must be a character vector of non-empty column names (or NULL)."
      )
    }
  }
  if (!is.numeric(qCutoffs)) {
    stop("qCutoffs must be numeric.")
  }

  minPositive <- as.numeric(minPositive)[1]
  if (is.na(minPositive) || minPositive < 0) {
    stop("minPositive must be a single non-negative numeric value.")
  }

  if (!is.null(rt)) {
    if (is.logical(rt)) {
      if (length(rt) != 1 || is.na(rt)) {
        stop("rt must be a single logical value (or NULL).")
      }
    } else if (is.character(rt)) {
      if (length(rt) != 1 || is.na(rt) || rt == "") {
        stop("rt must be a single non-empty character value (or NULL).")
      }
    } else if (is.list(rt)) {
      allowed <- c("col")
      extra <- setdiff(names(rt), allowed)
      if (length(extra) > 0) {
        stop("rt list only supports the field: col")
      }
      if (!is.null(rt$col)) {
        rt$col <- validateColName(rt$col, "rt$col")
      }
    } else {
      stop("rt must be NULL, logical, character, or list(col = \"...\").")
    }
  }

  q <- if (isTRUE(qEnable)) {
    list(enable = TRUE, cols = qCols, cutoffs = qCutoffs)
  } else {
    FALSE
  }

  list(
    level = level,
    sep = sep,
    columns = list(
      sample = sampleCol,
      feature = featureCol,
      quantity = quantityCol,
      extra = extraCols
    ),
    filters = list(
      decoy = decoy,
      q = q,
      min_positive = minPositive
    ),
    rt = rt
  )
}

diannStripPath <- function(paths) {
  sub("^.*[\\\\/]", "", paths)
}

diannStripKnownExtension <- function(fileNames) {
  sub("\\.(raw|mzml|d|wiff)$", "", fileNames, ignore.case = TRUE)
}

diannCleanSampleName <- function(sampleNames) {
  diannStripKnownExtension(diannStripPath(sampleNames))
}

diannCandidateSampleColumns <- function(columnNames) {
  grepl("[\\\\/]|\\.(raw|mzml|d|wiff)$", columnNames, ignore.case = TRUE)
}

diannRenameSampleColumnsForDesign <- function(dataFrame, designSampleNames) {
  if (is.null(designSampleNames) || length(designSampleNames) == 0) {
    return(dataFrame)
  }

  designSampleNames <- as.character(designSampleNames)
  if (all(designSampleNames %in% colnames(dataFrame))) {
    return(dataFrame)
  }

  candidateCols <- diannCandidateSampleColumns(colnames(dataFrame))
  if (!any(candidateCols)) {
    return(dataFrame)
  }

  cleaned <- colnames(dataFrame)
  cleaned[candidateCols] <- diannCleanSampleName(cleaned[candidateCols])

  dupNames <- unique(cleaned[candidateCols][duplicated(cleaned[candidateCols])])
  if (length(dupNames) > 0) {
    stop(
      "DIANN sample columns are not unique after stripping paths/extensions. ",
      "Duplicate sample names include: ",
      paste(utils::head(dupNames, 10), collapse = ", "),
      "\n",
      "Provide unique sample names in DIA-NN export, or use full file paths in the design matrix."
    )
  }

  if (all(designSampleNames %in% cleaned)) {
    colnames(dataFrame) <- cleaned
  }

  dataFrame
}

diannChooseReportSpec <- function(
  reportColumns,
  diannLevel = c("auto", "protein", "precursor"),
  diannSampleCol = NULL,
  diannFeatureCol = NULL,
  diannQuantityCol = NULL
) {
  diannLevel <- match.arg(diannLevel)

  sampleCol <- if (!is.null(diannSampleCol)) {
    if (!(diannSampleCol %in% reportColumns)) {
      stop(
        "DIANN report file is missing requested sample column: ",
        diannSampleCol
      )
    }
    diannSampleCol
  } else if ("Run" %in% reportColumns) {
    "Run"
  } else if ("File.Name" %in% reportColumns) {
    "File.Name"
  } else {
    stop("DIANN report file is missing both 'Run' and 'File.Name' columns")
  }

  proteinFeatureCol <- "Protein.Group"
  precursorFeatureCol <- "Precursor.Id"
  proteinQuantityCandidates <- c("PG.MaxLFQ", "PG.Normalised", "PG.Quantity")
  precursorQuantityCandidates <- c(
    "Precursor.Normalised",
    "Precursor.Quantity",
    "Precursor.Translated"
  )

  inferQuantity <- function(featureCol) {
    if (identical(featureCol, proteinFeatureCol)) {
      diannSelectFirstPresent(proteinQuantityCandidates, reportColumns)
    } else if (identical(featureCol, precursorFeatureCol)) {
      diannSelectFirstPresent(precursorQuantityCandidates, reportColumns)
    } else {
      NULL
    }
  }

  featureCol <- NULL
  quantityCol <- NULL

  if (!is.null(diannFeatureCol)) {
    if (!(diannFeatureCol %in% reportColumns)) {
      stop(
        "DIANN report file is missing requested feature column: ",
        diannFeatureCol
      )
    }
    featureCol <- diannFeatureCol
    if (!is.null(diannQuantityCol)) {
      if (!(diannQuantityCol %in% reportColumns)) {
        stop(
          "DIANN report file is missing requested quantity column: ",
          diannQuantityCol
        )
      }
      quantityCol <- diannQuantityCol
    } else {
      quantityCol <- inferQuantity(featureCol)
      if (is.null(quantityCol)) {
        stop(
          "Could not infer DIANN quantity column for feature column '",
          featureCol,
          "'. ",
          "Provide diannQuantityCol explicitly."
        )
      }
    }
  } else if (!is.null(diannQuantityCol)) {
    if (!(diannQuantityCol %in% reportColumns)) {
      stop(
        "DIANN report file is missing requested quantity column: ",
        diannQuantityCol
      )
    }

    quantityCol <- diannQuantityCol
    featureCol <- if (identical(diannLevel, "protein")) {
      proteinFeatureCol
    } else if (identical(diannLevel, "precursor")) {
      precursorFeatureCol
    } else if (proteinFeatureCol %in% reportColumns) {
      proteinFeatureCol
    } else {
      precursorFeatureCol
    }

    if (!(featureCol %in% reportColumns)) {
      stop(
        "DIANN report file is missing requested feature column: ",
        featureCol
      )
    }
  } else {
    if (identical(diannLevel, "protein") || identical(diannLevel, "auto")) {
      if (proteinFeatureCol %in% reportColumns) {
        featureCol <- proteinFeatureCol
        quantityCol <- diannSelectFirstPresent(
          proteinQuantityCandidates,
          reportColumns
        )
      }
    }

    if (
      is.null(quantityCol) &&
        (identical(diannLevel, "precursor") || identical(diannLevel, "auto"))
    ) {
      if (precursorFeatureCol %in% reportColumns) {
        featureCol <- precursorFeatureCol
        quantityCol <- diannSelectFirstPresent(
          precursorQuantityCandidates,
          reportColumns
        )
      }
    }

    if (is.null(quantityCol) || !(featureCol %in% reportColumns)) {
      stop(
        "Could not infer DIANN report feature/quantity columns. ",
        "Expected protein-level columns like 'Protein.Group' + 'PG.Quantity'/'PG.MaxLFQ' ",
        "or precursor-level columns like 'Precursor.Id' + 'Precursor.Quantity'."
      )
    }
  }

  extraColsCandidates <- if (identical(featureCol, proteinFeatureCol)) {
    c(
      "Protein.Group",
      "Protein.Ids",
      "Protein.Names",
      "Genes",
      "First.Protein.Description"
    )
  } else if (identical(featureCol, precursorFeatureCol)) {
    c(
      "Precursor.Id",
      "Modified.Sequence",
      "Stripped.Sequence",
      "Precursor.Charge",
      "Protein.Group",
      "Protein.Names",
      "Genes",
      "Proteotypic"
    )
  } else {
    c(featureCol)
  }

  decoyCol <- if ("Decoy" %in% reportColumns) "Decoy" else NULL

  list(
    sampleCol = sampleCol,
    featureCol = featureCol,
    quantityCol = quantityCol,
    extraCols = intersect(extraColsCandidates, reportColumns),
    decoyCol = decoyCol
  )
}

diannReadReportTSV <- function(filePath, sep = "\t", selectCols) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    suppressWarnings(
      data.table::fread(
        filePath,
        sep = sep,
        select = selectCols,
        data.table = FALSE,
        showProgress = FALSE
      )
    )
  } else {
    header <- diannReadHeader(filePath, sep = sep)
    colClasses <- rep("NULL", length(header))
    keep <- header %in% selectCols
    colClasses[keep] <- "character"
    utils::read.table(
      filePath,
      sep = sep,
      header = TRUE,
      quote = "",
      comment.char = "",
      check.names = FALSE,
      stringsAsFactors = FALSE,
      colClasses = colClasses,
      na.strings = c("NA", "null", "")
    )
  }
}

diannReadReportParquet <- function(filePath, selectCols) {
  if (!requireNamespace("arrow", quietly = TRUE)) {
    stop(
      "Reading DIANN parquet files requires the optional 'arrow' package. ",
      "Install it, or export DIANN output as TSV instead."
    )
  }
  reader <- arrow::ParquetFileReader$create(filePath)
  schemaNames <- reader$GetSchema()$names

  selectCols <- unique(selectCols)
  selectIndices <- match(selectCols, schemaNames)
  missing <- selectCols[is.na(selectIndices)]
  if (length(missing) > 0) {
    stop(
      "DIANN parquet file is missing expected columns: ",
      paste(utils::head(missing, 10), collapse = ", ")
    )
  }

  arrow::read_parquet(
    filePath,
    col_select = as.integer(selectIndices),
    as_data_frame = TRUE
  )
}

diannFilterDecoys <- function(reportDf, decoyCol = "Decoy") {
  if (is.null(decoyCol) || !(decoyCol %in% colnames(reportDf))) {
    return(reportDf)
  }

  decoy <- reportDf[[decoyCol]]
  keep <- rep(TRUE, length(decoy))

  if (is.logical(decoy)) {
    keep <- is.na(decoy) | !decoy
  } else if (is.numeric(decoy)) {
    keep <- is.na(decoy) | decoy == 0
  } else {
    decoyStr <- tolower(as.character(decoy))
    keep <- is.na(decoyStr) | !(decoyStr %in% c("1", "true", "t", "yes", "y"))
  }

  reportDf[keep, , drop = FALSE]
}

diannDefaultQValueCols <- function(featureCol) {
  if (identical(featureCol, "Precursor.Id")) {
    return(c("Q.Value", "Lib.Q.Value", "Lib.PG.Q.Value"))
  }
  if (identical(featureCol, "Protein.Group")) {
    return(c("PG.Q.Value", "Global.PG.Q.Value", "Lib.PG.Q.Value"))
  }
  c("Q.Value", "Lib.Q.Value")
}

diannResolveQValueCols <- function(diannQValueCols, reportColumns, featureCol) {
  if (is.null(diannQValueCols)) {
    candidates <- diannDefaultQValueCols(featureCol)
    return(intersect(candidates, reportColumns))
  }

  diannQValueCols <- as.character(diannQValueCols)
  if (length(diannQValueCols) == 0) {
    return(character())
  }
  if (
    length(diannQValueCols) == 1 && identical(tolower(diannQValueCols), "auto")
  ) {
    candidates <- diannDefaultQValueCols(featureCol)
    return(intersect(candidates, reportColumns))
  }
  if (
    length(diannQValueCols) == 1 && identical(tolower(diannQValueCols), "none")
  ) {
    return(character())
  }

  missing <- setdiff(diannQValueCols, reportColumns)
  if (length(missing) > 0) {
    stop(
      "DIANN report file is missing requested q-value columns: ",
      paste(utils::head(missing, 10), collapse = ", ")
    )
  }
  diannQValueCols
}

diannResolveQValueCutoffs <- function(diannQValueCols, diannQValueCutoffs) {
  if (length(diannQValueCols) == 0) {
    return(numeric())
  }
  if (is.null(diannQValueCutoffs) || length(diannQValueCutoffs) == 0) {
    return(rep(0.01, length(diannQValueCols)))
  }
  diannQValueCutoffs <- as.numeric(diannQValueCutoffs)
  diannQValueCutoffs <- diannQValueCutoffs[!is.na(diannQValueCutoffs)]
  if (length(diannQValueCutoffs) == 0) {
    return(rep(0.01, length(diannQValueCols)))
  }
  if (length(diannQValueCutoffs) != length(diannQValueCols)) {
    diannQValueCutoffs <- rep_len(
      diannQValueCutoffs[1],
      length(diannQValueCols)
    )
  }
  diannQValueCutoffs
}

diannFilterByQValue <- function(reportDf, qCols, qCutoffs) {
  if (length(qCols) == 0) {
    return(reportDf)
  }

  keep <- rep(TRUE, nrow(reportDf))
  for (i in seq_along(qCols)) {
    q <- reportDf[[qCols[i]]]
    if (!is.numeric(q)) {
      q <- suppressWarnings(as.numeric(q))
    }
    keep[q > qCutoffs[i]] <- FALSE
  }

  reportDf[keep, , drop = FALSE]
}

diannResolveRTCol <- function(diannRTCol, reportColumns, featureCol) {
  if (!identical(featureCol, "Precursor.Id")) {
    return(NULL)
  }
  if (is.null(diannRTCol) || length(diannRTCol) == 0) {
    return(NULL)
  }
  diannRTCol <- as.character(diannRTCol)[1]
  if (is.na(diannRTCol) || diannRTCol == "") {
    return(NULL)
  }
  if (diannRTCol %in% reportColumns) {
    return(diannRTCol)
  }
  NULL
}

diannWarnIfAutoInferenceIsAmbiguous <- function(reportColumns, opts, spec) {
  if (
    is.null(opts) ||
      is.null(spec) ||
      !identical(opts$level, "auto") ||
      !is.null(opts$featureCol) ||
      !is.null(opts$quantityCol)
  ) {
    return(invisible(NULL))
  }

  proteinFeatureCol <- "Protein.Group"
  precursorFeatureCol <- "Precursor.Id"

  proteinQuantityCandidates <- c("PG.MaxLFQ", "PG.Normalised", "PG.Quantity")
  precursorQuantityCandidates <- c(
    "Precursor.Normalised",
    "Precursor.Quantity",
    "Precursor.Translated"
  )

  proteinPossible <- proteinFeatureCol %in% reportColumns &&
    any(proteinQuantityCandidates %in% reportColumns)
  precursorPossible <- precursorFeatureCol %in% reportColumns &&
    any(precursorQuantityCandidates %in% reportColumns)

  if (!proteinPossible || !precursorPossible) {
    return(invisible(NULL))
  }

  inferred <- paste0(spec$featureCol, " + ", spec$quantityCol)

  warning(
    "DIANN report contains both protein-level and precursor-level quantities. ",
    "Using inferred columns by default: ",
    inferred,
    ". If you intended precursor-level analysis, set `inputOptions=diannInputOptions(level=\"precursor\", quantityCol=\"Precursor.Quantity\")` ",
    "(or set level/feature/quantity explicitly).",
    call. = FALSE
  )

  invisible(NULL)
}

diannReportToWide <- function(
  reportDf,
  sampleCol,
  featureCol,
  quantityCol,
  extraCols,
  designSampleNames = NULL,
  diannMinPositive = 0,
  rtCol = NULL
) {
  reportDf[[sampleCol]] <- as.character(reportDf[[sampleCol]])
  reportDf[[featureCol]] <- as.character(reportDf[[featureCol]])

  diannMinPositive <- diannResolveMinPositive(diannMinPositive)

  quantities <- reportDf[[quantityCol]]
  if (!is.numeric(quantities)) {
    quantities <- suppressWarnings(as.numeric(quantities))
  }
  quantities[quantities == 0] <- NA_real_
  if (diannMinPositive > 0) {
    quantities[!is.na(quantities) & quantities < diannMinPositive] <- NA_real_
  }

  samples <- unique(reportDf[[sampleCol]])
  if (!is.null(designSampleNames) && !all(designSampleNames %in% samples)) {
    cleanedSamples <- diannCleanSampleName(samples)
    sampleMap <- stats::setNames(cleanedSamples, samples)
    mapped <- sampleMap[reportDf[[sampleCol]]]
    mapped <- as.character(mapped)
    mapped[is.na(mapped)] <- reportDf[[sampleCol]][is.na(mapped)]

    mappedUnique <- unique(mapped)
    if (all(designSampleNames %in% mappedUnique)) {
      reportDf[[sampleCol]] <- mapped
      samples <- unique(reportDf[[sampleCol]])
    }
  }

  features <- unique(reportDf[[featureCol]])
  wide <- matrix(NA_real_, nrow = length(features), ncol = length(samples))
  colnames(wide) <- samples

  sampleIndex <- match(reportDf[[sampleCol]], samples)
  featureIndex <- match(reportDf[[featureCol]], features)
  linearIndex <- featureIndex + (sampleIndex - 1L) * length(features)
  keep <- !is.na(linearIndex) & !is.na(quantities)
  linearIndex <- linearIndex[keep]
  quantities <- quantities[keep]

  if (length(linearIndex) > 0) {
    if (anyDuplicated(linearIndex)) {
      maxByIndex <- tapply(quantities, linearIndex, max, na.rm = TRUE)
      wide[as.integer(names(maxByIndex))] <- maxByIndex
    } else {
      wide[linearIndex] <- quantities
    }
  }

  dedup <- !duplicated(reportDf[[featureCol]])
  annotation <- reportDf[dedup, extraCols, drop = FALSE]
  annotation <- annotation[
    match(features, annotation[[featureCol]]),
    ,
    drop = FALSE
  ]

  if (!is.null(rtCol) && (rtCol %in% colnames(reportDf))) {
    rtValues <- reportDf[[rtCol]]
    if (!is.numeric(rtValues)) {
      rtValues <- suppressWarnings(as.numeric(rtValues))
    }
    rtByFeature <- tapply(
      rtValues,
      reportDf[[featureCol]],
      stats::median,
      na.rm = TRUE
    )
    annotation[["RT"]] <- as.numeric(rtByFeature[features])
  }

  data.frame(
    annotation,
    as.data.frame(wide, check.names = FALSE),
    check.names = FALSE
  )
}

readDiannToDataFrame <- function(
  filePath,
  sep = "\t",
  designSampleNames = NULL,
  inputOptions = NULL
) {
  opts <- diannNormalizeInputOptions(inputOptions, sep = sep)
  sep <- opts$sep

  if (diannIsParquet(filePath)) {
    if (!requireNamespace("arrow", quietly = TRUE)) {
      stop(
        "Reading DIANN parquet files requires the optional 'arrow' package. ",
        "Install it, or export DIANN output as TSV instead."
      )
    }

    reportColumns <- arrow::ParquetFileReader$create(filePath)$GetSchema()$names
    spec <- diannChooseReportSpec(
      reportColumns,
      diannLevel = opts$level,
      diannSampleCol = opts$sampleCol,
      diannFeatureCol = opts$featureCol,
      diannQuantityCol = opts$quantityCol
    )
    if (isTRUE(getOption("NormalyzerDE.warnDiannAutoAmbiguous", FALSE))) {
      diannWarnIfAutoInferenceIsAmbiguous(reportColumns, opts, spec)
    }

    extraCols <- spec$extraCols
    if (!is.null(opts$extraCols)) {
      requested <- as.character(opts$extraCols)
      if (length(requested) == 0) {
        extraCols <- character()
      } else {
        missingExtra <- setdiff(requested, reportColumns)
        if (length(missingExtra) > 0) {
          warning(
            "DIANN report file is missing some requested extra columns: ",
            paste(utils::head(missingExtra, 10), collapse = ", ")
          )
        }
        extraCols <- intersect(requested, reportColumns)
      }
      extraCols <- unique(c(spec$featureCol, extraCols))
    }

    qCols <- if (isTRUE(opts$filterQValue)) {
      diannResolveQValueCols(opts$qCols, reportColumns, spec$featureCol)
    } else {
      character()
    }
    qCutoffs <- diannResolveQValueCutoffs(qCols, opts$qCutoffs)
    rtCol <- diannResolveRTCol(opts$rtCol, reportColumns, spec$featureCol)

    selectCols <- unique(c(
      spec$sampleCol,
      spec$featureCol,
      spec$quantityCol,
      extraCols,
      spec$decoyCol,
      qCols,
      rtCol
    ))
    reportDf <- diannReadReportParquet(filePath, selectCols = selectCols)
    if (isTRUE(opts$filterDecoy)) {
      reportDf <- diannFilterDecoys(reportDf, spec$decoyCol)
    }
    if (length(qCols) > 0) {
      reportDf <- diannFilterByQValue(
        reportDf,
        qCols = qCols,
        qCutoffs = qCutoffs
      )
    }
    return(
      diannReportToWide(
        reportDf,
        sampleCol = spec$sampleCol,
        featureCol = spec$featureCol,
        quantityCol = spec$quantityCol,
        extraCols = extraCols,
        designSampleNames = designSampleNames,
        diannMinPositive = opts$minPositive,
        rtCol = rtCol
      )
    )
  }

  header <- diannReadHeader(filePath, sep = sep)
  if (diannIsReportHeader(header)) {
    spec <- diannChooseReportSpec(
      header,
      diannLevel = opts$level,
      diannSampleCol = opts$sampleCol,
      diannFeatureCol = opts$featureCol,
      diannQuantityCol = opts$quantityCol
    )
    if (isTRUE(getOption("NormalyzerDE.warnDiannAutoAmbiguous", FALSE))) {
      diannWarnIfAutoInferenceIsAmbiguous(header, opts, spec)
    }

    extraCols <- spec$extraCols
    if (!is.null(opts$extraCols)) {
      requested <- as.character(opts$extraCols)
      if (length(requested) == 0) {
        extraCols <- character()
      } else {
        missingExtra <- setdiff(requested, header)
        if (length(missingExtra) > 0) {
          warning(
            "DIANN report file is missing some requested extra columns: ",
            paste(utils::head(missingExtra, 10), collapse = ", ")
          )
        }
        extraCols <- intersect(requested, header)
      }
      extraCols <- unique(c(spec$featureCol, extraCols))
    }

    qCols <- if (isTRUE(opts$filterQValue)) {
      diannResolveQValueCols(opts$qCols, header, spec$featureCol)
    } else {
      character()
    }
    qCutoffs <- diannResolveQValueCutoffs(qCols, opts$qCutoffs)
    rtCol <- diannResolveRTCol(opts$rtCol, header, spec$featureCol)

    selectCols <- unique(c(
      spec$sampleCol,
      spec$featureCol,
      spec$quantityCol,
      extraCols,
      spec$decoyCol,
      qCols,
      rtCol
    ))
    reportDf <- diannReadReportTSV(filePath, sep = sep, selectCols = selectCols)
    if (isTRUE(opts$filterDecoy)) {
      reportDf <- diannFilterDecoys(reportDf, spec$decoyCol)
    }
    if (length(qCols) > 0) {
      reportDf <- diannFilterByQValue(
        reportDf,
        qCols = qCols,
        qCutoffs = qCutoffs
      )
    }
    return(
      diannReportToWide(
        reportDf,
        sampleCol = spec$sampleCol,
        featureCol = spec$featureCol,
        quantityCol = spec$quantityCol,
        extraCols = extraCols,
        designSampleNames = designSampleNames,
        diannMinPositive = opts$minPositive,
        rtCol = rtCol
      )
    )
  }

  matrixDf <- utils::read.table(
    filePath,
    sep = sep,
    header = TRUE,
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  matrixDf <- diannRenameSampleColumnsForDesign(matrixDf, designSampleNames)
  sampleCols <- intersect(as.character(designSampleNames), colnames(matrixDf))
  if (length(sampleCols) > 0) {
    for (colName in sampleCols) {
      if (is.numeric(matrixDf[[colName]])) {
        matrixDf[[colName]][matrixDf[[colName]] == 0] <- NA_real_
        if (opts$minPositive > 0) {
          matrixDf[[colName]][
            !is.na(matrixDf[[colName]]) & matrixDf[[colName]] < opts$minPositive
          ] <- NA_real_
        }
      }
    }
  }

  matrixDf
}
