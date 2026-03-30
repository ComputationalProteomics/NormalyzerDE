#' Load raw data into dataframe
#'
#' General function which allows specifying different types of input data
#' including "proteios", "maxquantpep" (peptide output from MaxQuant) and
#' "maxquantprot" (protein output from MaxQuant) formats.
#'
#' @param dataPath File path to data matrix.
#' @param inputFormat If input is given in standard NormalyzerDE format,
#' Proteios format or in MaxQuant protein or peptide format
#' @param inputOptions Optional list of input-reader options. For example,
#'   \code{defaultInputOptions(sep=",")} to read comma-separated files.
#' @return rawData Raw data loaded into data frame
#' @export
#' @examples
#' data_path <- system.file(package="NormalyzerDE", "extdata", "tiny_data.tsv")
#' df <- loadData(data_path)
#'
#' @seealso \code{\link{defaultInputOptions}}, \code{\link{proteiosInputOptions}},
#'   and \code{\link{maxQuantInputOptions}} for delimiter configuration.
loadData <- function(dataPath, inputFormat = "default", inputOptions = NULL) {
  sep <- "\t"
  if (is.list(inputOptions) && !is.null(inputOptions$sep)) {
    sep <- as.character(inputOptions$sep)[1]
  }

  if (inputFormat == "default") {
    rawData <- loadRawDataFromFile(dataPath, sep = sep)
  } else if (inputFormat == "proteios") {
    rawData <- proteiosToNormalyzer(dataPath, sep = sep)
  } else if (inputFormat == "maxquantpep") {
    rawData <- maxQuantToNormalyzer(dataPath, protLevel = FALSE, sep = sep)
  } else if (inputFormat == "maxquantprot") {
    rawData <- maxQuantToNormalyzer(dataPath, protLevel = TRUE, sep = sep)
  } else {
    valids <- c("default", "proteios", "maxquantpep", "maxquantprot")
    cli::cli_abort(
      c(
        "Unknown {.arg inputFormat}: {.val {inputFormat}}.",
        i = "Valid values: {paste(valids, collapse = \", \")}."
      ),
      class = "normalyzerde_error",
      call = NULL
    )
  }

  rawData
}

#' Create input options for default (delimited) matrices
#'
#' Helper to construct an \code{inputOptions} list for \code{inputFormat =
#' "default"}. Currently this supports selecting the input delimiter.
#'
#' @param sep Field separator used when reading the data matrix.
#' @return A list suitable to pass as \code{inputOptions}.
#' @export
#' @examples
#' defaultInputOptions(sep = ",")
defaultInputOptions <- function(sep = "\t") {
  call <- sys.call()
  argNames <- names(call)[-1]
  argNames <- argNames[nzchar(argNames)]
  unknown <- setdiff(argNames, names(formals(sys.function())))
  if (length(unknown) > 0) {
    cli::cli_abort(
      "Unknown argument(s): {paste(unknown, collapse = \", \")}.",
      class = "normalyzerde_error",
      call = NULL
    )
  }

  if (!is.character(sep) || length(sep) != 1 || is.na(sep) || sep == "") {
    cli::cli_abort(
      "{.arg sep} must be a single non-empty character value.",
      class = "normalyzerde_error",
      call = NULL
    )
  }
  list(sep = sep)
}

#' Load raw design into dataframe
#'
#' Takes a design path, loads the matrix and ensures that the sample column
#' is in character format and that the group column is in factor format.
#'
#' @param designPath File path to design matrix.
#' @param sampleCol Column name for column containing sample names.
#' @param groupCol Column name for column containing condition levels.
#' @return designMatrix Design data loaded into data frame
#' @export
#' @examples
#' design_path <- system.file(package="NormalyzerDE", "extdata", "tiny_design.tsv")
#' df <- loadDesign(design_path)
loadDesign <- function(designPath, sampleCol = "sample", groupCol = "group") {
  designMatrix <- utils::read.table(
    designPath,
    sep = "\t",
    stringsAsFactors = FALSE,
    header = TRUE,
    comment.char = "",
    check.names = FALSE
  )

  if (
    !(sampleCol %in% colnames(designMatrix)) ||
      !(groupCol %in% colnames(designMatrix))
  ) {
    cli::cli_abort(
      c(
        "Both {.arg sampleCol} and {.arg groupCol} must be present in the design matrix header.",
        i = "{.arg sampleCol}: {.val {sampleCol}}",
        i = "{.arg groupCol}: {.val {groupCol}}",
        i = "Design matrix header: {paste(colnames(designMatrix), collapse = \", \")}."
      ),
      class = "normalyzerde_error",
      call = NULL
    )
  }

  designMatrix[, sampleCol] <- as.character(designMatrix[, sampleCol])
  designMatrix[, groupCol] <- as.factor(as.character(designMatrix[, groupCol]))
  designMatrix
}

#' Prepare a SummarizedExperiment object for normalization
#'
#' @param dataPath File path to data matrix.
#' @param designPath File path to design matrix.
#' @param inputFormat Type of matrix for data, can be either 'default',
#'   'proteios', 'maxquantprot', 'maxquantpep' or 'diann'
#' @param zeroToNA If TRUE zeroes in the data is automatically converted to
#'   NA values
#' @param sampleColName Column name for column containing sample names
#' @param groupColName Column name for column containing condition levels
#' @param inputOptions Optional list of input-reader options. For DIA-NN, use
#'   \code{\link{diannInputOptions}}. For delimited inputs, use
#'   \code{\link{defaultInputOptions}}, \code{\link{proteiosInputOptions}}, or
#'   \code{\link{maxQuantInputOptions}} to configure the delimiter.
#' @return experimentObj SummarizedExperiment object containing the data, design
#'   and annotation information
#' @export
#' @examples
#' data_path <- system.file(package="NormalyzerDE", "extdata", "tiny_data.tsv")
#' design_path <- system.file(package="NormalyzerDE", "extdata", "tiny_design.tsv")
#' df <- setupRawDataObject(data_path, design_path)
setupRawDataObject <- function(
  dataPath,
  designPath,
  inputFormat = "default",
  zeroToNA = FALSE,
  sampleColName = "sample",
  groupColName = "group",
  inputOptions = NULL
) {
  rawDesign <- loadDesign(
    designPath,
    sampleCol = sampleColName,
    groupCol = groupColName
  )

  if (identical(inputFormat, "diann")) {
    fullDf <- readDiannToDataFrame(
      dataPath,
      designSampleNames = rawDesign[[sampleColName]],
      inputOptions = inputOptions
    )
    rawData <- as.matrix(rbind(colnames(fullDf), fullDf))
  } else {
    rawData <- loadData(
      dataPath,
      inputFormat = inputFormat,
      inputOptions = inputOptions
    )
  }

  rdf <- rawData[2:nrow(rawData), ]
  colnames(rdf) <- rawData[1, ]

  verifyDesignMatrix(rdf, rawDesign, sampleColName)

  rawDesign[[sampleColName]] <- as.character(rawDesign[[sampleColName]])

  sdf <- rdf[, as.character(rawDesign[[sampleColName]])]
  adf <- rdf[,
    !(colnames(rdf) %in% as.character(rawDesign[[sampleColName]])),
    drop = FALSE
  ]

  experimentObj <- SummarizedExperiment::SummarizedExperiment(
    assays = list(raw = as.matrix(sdf)),
    rowData = adf,
    colData = rawDesign,
    metadata = list(sample = sampleColName, group = groupColName)
  )
  experimentObj
}

#' Prepare SummarizedExperiment object for statistics data
#'
#' @param dataPath Path to raw data matrix
#' @param designPath Path to design matrix
#' @param sampleColName Name for column in design matrix containing sample names
#' @param inputFormat Type of input format for \code{dataPath}. Supports
#'   \code{"default"} and \code{"diann"}.
#' @param inputOptions Optional list of input-reader options. Used when
#'   \code{inputFormat="diann"} to select columns and filters (q-value cutoffs,
#'   decoys, min-positive threshold, RT summarization, etc.).
#' @return experimentObj Prepared instance of SummarizedExperiment
#' @export
#' @examples
#' data_path <- system.file(package="NormalyzerDE", "extdata", "tiny_data.tsv")
#' design_path <- system.file(package="NormalyzerDE", "extdata", "tiny_design.tsv")
#' sumExpObj <- setupRawContrastObject(data_path, design_path, "sample")
setupRawContrastObject <- function(
  dataPath,
  designPath,
  sampleColName,
  inputFormat = "default",
  inputOptions = NULL
) {
  designDf <- tryCatch(
    utils::read.csv(
      designPath,
      sep = "\t",
      stringsAsFactors = FALSE,
      quote = "",
      comment.char = "",
      check.names = FALSE
    ),
    error = function(e) {
      cli::cli_abort(
        c(
          "Failed to read {.arg designPath}: {.path {designPath}}.",
          i = "{conditionMessage(e)}"
        ),
        class = "normalyzerde_error",
        call = NULL
      )
    }
  )

  if (identical(inputFormat, "diann")) {
    fullDf <- tryCatch(
      readDiannToDataFrame(
        dataPath,
        designSampleNames = designDf[[sampleColName]],
        inputOptions = inputOptions
      ),
      error = function(e) {
        cli::cli_abort(
          c(
            "Failed to read {.arg dataPath}: {.path {dataPath}}.",
            i = "{conditionMessage(e)}"
          ),
          class = "normalyzerde_error",
          call = NULL
        )
      }
    )
  } else {
    fullDf <- tryCatch(
      utils::read.csv(
        dataPath,
        sep = "\t",
        stringsAsFactors = FALSE,
        quote = "",
        comment.char = "",
        check.names = FALSE
      ),
      error = function(e) {
        cli::cli_abort(
          c(
            "Failed to read {.arg dataPath}: {.path {dataPath}}.",
            i = "{conditionMessage(e)}"
          ),
          class = "normalyzerde_error",
          call = NULL
        )
      }
    )
  }

  verifyDesignMatrix(fullDf, designDf, sampleColName)

  sdf <- fullDf[, designDf[[sampleColName]]]
  adf <- fullDf[,
    !(colnames(fullDf) %in% as.character(designDf[[sampleColName]])),
    drop = FALSE
  ]

  experimentObj <- SummarizedExperiment::SummarizedExperiment(
    assays = list(raw = as.matrix(sdf)),
    colData = designDf,
    rowData = adf,
    metadata = list(sample = sampleColName)
  )
  experimentObj
}

#' Verify that input data is in correct format, and if so, return a generated
#'  NormalyzerDE data object from that input data
#'
#' This function performs a number of checks on the input data and provides
#' informative error messages if the data isn't fulfilling the required format.
#' Checks include verifying that the design matrix matches to the data matrix,
#' that the data matrix contains valid numbers and that samples have enough
#' values for analysis
#'
#' @param jobName Name of ongoing run.
#' @param summarizedExp Summarized experiment input object
#' @param threshold Minimum number of features.
#' @param omitSamples Automatically omit invalid samples from analysis.
#' @param requireReplicates Require there to be at least to samples per
#'        condition
#' @param quiet Don't print output messages during processing
#' @param noLogTransform Don't log-transform the provided data
#' @param tinyRunThres If less features in run, a limited run is performed
#'
#' @return Normalyzer data object representing verified input data.
#' @export
#' @examples
#' data(example_summarized_experiment)
#' normObj <- getVerifiedNormalyzerObject("job_name", example_summarized_experiment)
getVerifiedNormalyzerObject <- function(
  jobName,
  summarizedExp,
  threshold = 15,
  omitSamples = FALSE,
  requireReplicates = TRUE,
  quiet = FALSE,
  noLogTransform = FALSE,
  tinyRunThres = 50
) {
  SummarizedExperiment::assay(summarizedExp) <- preprocessData(
    SummarizedExperiment::assay(summarizedExp),
    quiet = quiet
  )
  summarizedExp <- filterOnlyNARows(summarizedExp)

  metadata <- S4Vectors::metadata(summarizedExp)
  groupCol <- metadata$group
  sampleCol <- metadata$sample
  designMatrix <- as.data.frame(
    SummarizedExperiment::colData(summarizedExp),
    optional = TRUE
  )

  if (!groupCol %in% colnames(designMatrix)) {
    cli::cli_abort(
      "Given {.arg groupCol} {.val {groupCol}} was not present among design matrix columns.",
      class = "normalyzerde_error",
      call = NULL
    )
  }

  groups <- SummarizedExperiment::colData(summarizedExp)[[groupCol]]
  samples <- SummarizedExperiment::colData(summarizedExp)[[sampleCol]]
  dataMatrix <- SummarizedExperiment::assay(summarizedExp)
  annotationMatrix <- as.matrix(data.frame(
    lapply(
      SummarizedExperiment::rowData(summarizedExp),
      as.character
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  ))

  verifyDesignMatrix(dataMatrix, designMatrix, sampleCol)
  verifyValidNumbers(
    dataMatrix,
    groups,
    noLogTransform = noLogTransform,
    quiet = quiet
  )

  processedRawData <- getReplicateSortedData(dataMatrix, groups)
  # repSortedRawData <- getReplicateSortedData(dataMatrix, groups)
  # processedRawData <- preprocessData(repSortedRawData, quiet=quiet)

  lowCountSampleFiltered <- getLowCountSampleFiltered(
    processedRawData,
    groups,
    threshold = threshold,
    stopIfTooFew = !omitSamples
  )

  designMatrix <- designMatrix[
    designMatrix[[sampleCol]] %in% colnames(lowCountSampleFiltered),
  ]

  # If no samples left after omitting, stop
  verifyMultipleSamplesPresent(
    lowCountSampleFiltered,
    groups,
    requireReplicates = requireReplicates,
    quiet = quiet
  )

  validateSampleReplication(
    lowCountSampleFiltered,
    groups,
    requireReplicates = requireReplicates,
    quiet = quiet
  )

  nds <- NormalyzerDataset(
    jobName = jobName,
    designMatrix = designMatrix,
    rawData = processedRawData,
    annotationData = annotationMatrix,
    sampleNameCol = sampleCol,
    groupNameCol = groupCol,
    tinyRunThres = tinyRunThres,
    quiet = quiet
  )

  nds
}


filterOnlyNARows <- function(summarizedExp) {
  dataMatrix <- SummarizedExperiment::assay(summarizedExp)

  nonFullNAContr <- rowSums(is.na(SummarizedExperiment::assay(
    summarizedExp
  ))) !=
    ncol(summarizedExp)
  omittedCount <- sum(!nonFullNAContr)
  if (omittedCount > 0) {
    cli::cli_inform(
      c(i = "{omittedCount} entries with only NA values omitted")
    )
    summarizedExp <- summarizedExp[nonFullNAContr, ]
  }

  summarizedExp
}

#' Try reading raw Normalyzer matrix from provided filepath
#'
#' @param inputPath Path to Normalyzer data.
#' @return Table containing raw data from input file.
#' @keywords internal
loadRawDataFromFile <- function(inputPath, sep = "\t") {
  warningEnv <- new.env(parent = emptyenv())
  warningEnv$condition <- NULL

  rawData <- withCallingHandlers(
    try(
      as.matrix(
        utils::read.table(
          inputPath,
          header = FALSE,
          sep = sep,
          stringsAsFactors = FALSE,
          quote = "",
          comment.char = ""
        )
      ),
      silent = TRUE
    ),
    warning = function(w) {
      warningEnv$condition <- w
      invokeRestart("muffleWarning")
    }
  )

  if (!is.null(warningEnv$condition)) {
    cli::cli_abort(
      c(
        "An issue was encountered when attempting to load {.path {inputPath}}.",
        i = "{conditionMessage(warningEnv$condition)}",
        i = "Please investigate and provide a valid input file."
      ),
      class = "normalyzerde_error",
      call = NULL
    )
  }

  if (inherits(rawData, "try-error")) {
    errorCondition <- attr(rawData, "condition")
    errorMessage <- if (inherits(errorCondition, "condition")) {
      conditionMessage(errorCondition)
    } else {
      as.character(rawData)
    }

    cli::cli_abort(
      c(
        "Failed to read input file {.path {inputPath}}.",
        i = "{errorMessage}",
        i = "Please provide a valid input file."
      ),
      class = "normalyzerde_error",
      call = NULL
    )
  }

  rawData
}


#' Verify that input fields conform to the expected formats
#'
#' @param rawDataOnly Data frame with input data.
#' @param groups Condition levels for comparisons.
#' @return None
#' @keywords internal
verifyValidNumbers <- function(
  rawDataOnly,
  groups,
  noLogTransform = FALSE,
  quiet = FALSE
) {
  numericPattern <- if (isTRUE(noLogTransform)) {
    "[\\+\\-]?\\d+(\\.\\d+)?([eE][\\+\\-]?\\d+)?"
  } else {
    "\\d+(\\.\\d+)?([eE][\\+\\-]?\\d+)?"
  }

  # Fields expected to contain numbers in decimal or scientific notation, or containing NA or null
  validPatterns <- c(
    numericPattern,
    "NA",
    "\"NA\"",
    "null",
    ""
  )

  regexPattern <- sprintf("^(%s)$", paste(validPatterns, collapse = "|"))
  nonMatchIndices <- grep(
    regexPattern,
    rawDataOnly,
    perl = TRUE,
    ignore.case = TRUE,
    invert = TRUE
  )
  naIndices <- which(is.na(rawDataOnly))
  invalidNonNAIndices <- nonMatchIndices[!nonMatchIndices %in% naIndices]
  rowsWithIssues <- unique((invalidNonNAIndices - 1) %% nrow(rawDataOnly) + 1)

  if (length(invalidNonNAIndices) > 0) {
    invalidValues <- unique(rawDataOnly[unique(invalidNonNAIndices)])
    rowsPreview <- utils::head(rowsWithIssues, 10)
    firstIssueRow <- rawDataOnly[rowsWithIssues[1], ]

    cli::cli_abort(
      c(
        "Invalid values encountered in input data.",
        x = "Expected numeric values (dot-decimal, not comma) or NA/null fields.",
        x = "Invalid field values: {.val {invalidValues}}",
        i = "Rows with issues (showing up to 10): {.val {rowsPreview}}",
        i = "Content of the first row with issues: {.val {firstIssueRow}}"
      ),
      class = "normalyzerde_error",
      call = NULL
    )
  }

  if (!noLogTransform) {
    numericVals <- suppressWarnings(as.numeric(rawDataOnly))
    belowOneMatches <- which(
      is.finite(numericVals) &
        numericVals > 0 &
        numericVals < 1
    )

    if (length(belowOneMatches) > 0) {
      rowsWithIssues <- unique((belowOneMatches - 1) %% nrow(rawDataOnly) + 1)
      rowsPreview <- utils::head(rowsWithIssues, 10)
      firstIssueRow <- rawDataOnly[rowsWithIssues[1], ]

      cli::cli_abort(
        c(
          "Encountered below-one values in raw data.",
          x = "NormalyzerDE log2-transforms raw input by default; values below 1 suggest the matrix is already transformed or on an unexpected scale.",
          x = "Rows with issues (showing up to 10): {.val {rowsPreview}}",
          i = "Content of the first row with issues (row {.val {rowsWithIssues[1]}}): {.val {firstIssueRow}}",
          i = "If your input is already on the log2 scale, set {.arg noLogTransform}={.val TRUE}.",
          i = "Otherwise, confirm the expected raw-value scale before continuing."
        ),
        class = "normalyzerde_error",
        call = NULL
      )
    }

    warnIfLooksAlreadyLog2 <- function(mat, threshold = 50) {
      numericVals <- suppressWarnings(as.numeric(mat))
      finiteVals <- numericVals[is.finite(numericVals)]
      if (length(finiteVals) == 0) {
        return(invisible(NULL))
      }

      maxVal <- max(finiteVals)
      if (is.finite(maxVal) && maxVal < threshold) {
        maxValDisp <- signif(maxVal, 4)
        cli::cli_warn(
          c(
            "Input values may already be on the log2 scale (max finite value = {.val {maxValDisp}}).",
            i = "NormalyzerDE will log2-transform input by default.",
            i = "If your input is already on the log2 scale, set {.arg noLogTransform}={.val TRUE}."
          ),
          class = "normalyzerde_warning",
          call = NULL
        )
      }

      invisible(NULL)
    }

    warnIfLooksAlreadyLog2(rawDataOnly)
  }

  if (!quiet) {
    cli::cli_inform(c(v = "Input data checked. All fields are valid."))
  }
}


#' Verify a SummarizedExperiment contains matching samples
#'
#' Checks that the sample names (columns) present in a \code{SummarizedExperiment}
#' object match the sample IDs in its \code{colData}.
#'
#' @param summarizedExp SummarizedExperiment object to validate.
#' @param sampleCol Column in \code{colData} containing sample IDs.
#' @return None
#' @keywords internal
verifySummarizedExperiment <- function(summarizedExp, sampleCol) {
  fullMatrix <- cbind(
    data.frame(SummarizedExperiment::rowData(summarizedExp)),
    SummarizedExperiment::assay(summarizedExp)
  )

  designMatrix <- data.frame(SummarizedExperiment::colData(summarizedExp))

  verifyDesignMatrix(
    fullMatrix,
    designMatrix,
    sampleCol
  )
}


#' Verify that design matrix setup matches the data matrix
#'
#' @param fullMatrix Data frame with input data.
#' @param designMatrix Data frame with design setup.
#' @param sampleCol Column in design matrix containing sample IDs.
#'
#' @return None
#' @keywords internal
verifyDesignMatrix <- function(fullMatrix, designMatrix, sampleCol) {
  if (!(sampleCol %in% colnames(designMatrix))) {
    cli::cli_abort(
      c(
        "Design matrix header must contain {.arg sampleCol} name.",
        i = "Provided {.arg sampleCol}: {.val {sampleCol}}",
        i = "Design matrix header: {paste(colnames(designMatrix), collapse = \", \")}."
      ),
      class = "normalyzerde_error",
      call = NULL
    )
  }

  designColnames <- as.character(designMatrix[, sampleCol])

  if (!all(designColnames %in% colnames(fullMatrix))) {
    missing <- base::setdiff(designColnames, colnames(fullMatrix))
    cli::cli_abort(
      c(
        "Not all samples in the design matrix are present in the data matrix.",
        x = "Missing from data matrix header: {.val {missing}}",
        i = "Check that the data matrix column names match {.arg sampleCol} in the design matrix."
      ),
      class = "normalyzerde_error",
      call = NULL
    )
  }

  dataMatrix <- fullMatrix[, designColnames, drop = FALSE]
  dataColumns <- dataMatrix[, designColnames, drop = FALSE]

  if (length(designColnames) != ncol(dataColumns)) {
    cli::cli_abort(
      c(
        "Number of samples does not match the number of selected columns.",
        x = "Found {ncol(dataColumns)} column{?s}.",
        x = "Expected {length(designColnames)} column{?s}.",
        i = "Are all columns in the design matrix present in the data matrix?"
      ),
      class = "normalyzerde_error",
      call = NULL
    )
  }

  if (length(unique(designColnames)) != length(designColnames)) {
    duplicatedSamples <- unique(designColnames[duplicated(designColnames)])
    cli::cli_abort(
      c(
        "Sample labels must be unique.",
        x = "Duplicated sample labels: {.val {duplicatedSamples}}",
        i = "Found {length(unique(designColnames))} unique labels; expected {length(designColnames)}."
      ),
      class = "normalyzerde_error",
      call = NULL
    )
  }
}


#' Replace empty values (0 or empty field) with NA in input data
#'
#' @param dataMatrix Matrix with raw data.
#' @param quiet Don't show diagnostic messages
#' @return Parsed rawdata where 0 values are replaced with NA
#' @keywords internal
preprocessData <- function(dataMatrix, quiet = FALSE) {
  zeroFields <- length(dataMatrix[!is.na(dataMatrix) & dataMatrix == 0])
  emptyFields <- length(dataMatrix[!is.na(dataMatrix) & dataMatrix == ""])
  nullFields <- length(dataMatrix[!is.na(dataMatrix) & dataMatrix == "null"])

  if (zeroFields != 0) {
    if (!quiet) {
      cli::cli_inform(
        c(i = "{zeroFields} fields with '0' were replaced by 'NA'")
      )
    }
    dataMatrix[dataMatrix == 0] <- NA
  }

  if (emptyFields != 0) {
    if (!quiet) {
      cli::cli_inform(
        c(i = "{emptyFields} empty fields were replaced by 'NA'")
      )
    }
    dataMatrix[dataMatrix == ""] <- NA
  }

  if (nullFields != 0) {
    if (!quiet) {
      cli::cli_inform(
        c(i = "{nullFields} 'null' fields were replaced by 'NA'")
      )
    }
    dataMatrix[dataMatrix == "null"] <- NA
  }

  dataMatrix
}

#' Verify that samples contain at least a lowest number of values
#'
#' @param dataMatrix Data frame with processed input data.
#' @param groups Vector containing condition levels.
#' @param threshold Lowest number of allowed values in a column.
#' @param stopIfTooFew Abort run if lower than threshold number of values in
#'        column
#' @return None
#' @keywords internal
getLowCountSampleFiltered <- function(
  dataMatrix,
  groups,
  threshold = 15,
  stopIfTooFew = TRUE
) {
  sampleIndices <- seq_along(groups)
  sampleLabels <- colnames(dataMatrix)
  if (is.null(sampleLabels)) {
    sampleLabels <- as.character(sampleIndices)
  }

  numberOfValues <- colSums(!is.na(dataMatrix))
  notPassingThreshold <- which(numberOfValues < threshold)

  if (length(notPassingThreshold) == length(numberOfValues)) {
    failingSamples <- sampleLabels[notPassingThreshold]
    failingCounts <- numberOfValues[notPassingThreshold]
    countSummary <- paste0(failingSamples, "=", failingCounts, collapse = ", ")

    cli::cli_abort(
      c(
        "None of the samples had enough valid non-NA values.",
        x = "Threshold: {.val {threshold}} non-NA value{?s} per sample.",
        x = "Non-NA counts (sample=count): {countSummary}.",
        i = "You can try lowering the threshold via {.arg sampleAbundThres}.",
        i = "Be aware that this may lead to downstream crashes."
      ),
      class = "normalyzerde_error",
      call = NULL
    )
  } else if (length(notPassingThreshold) > 0) {
    if (stopIfTooFew) {
      failingSamples <- sampleLabels[notPassingThreshold]
      failingCounts <- numberOfValues[notPassingThreshold]
      countSummary <- paste0(failingSamples, "=", failingCounts, collapse = ", ")

      cli::cli_abort(
        c(
          "Some samples do not contain enough non-NA values.",
          x = "Threshold: {.val {threshold}} non-NA value{?s} per sample.",
          x = "Non-NA counts (sample=count): {countSummary}.",
          i = "You can force processing without these samples by setting {.arg omitLowAbundSamples}={.val TRUE}."
        ),
        class = "normalyzerde_error",
        call = NULL
      )
    } else {
      failingSamples <- sampleLabels[notPassingThreshold]
      failingCounts <- numberOfValues[notPassingThreshold]
      countSummary <- paste0(failingSamples, "=", failingCounts, collapse = ", ")

      cli::cli_warn(
        c(
          "Some samples do not contain enough non-NA values.",
          "!" = "Threshold: {.val {threshold}} non-NA value{?s} per sample.",
          "!" = "Non-NA counts (sample=count): {countSummary}.",
          i = "You can force processing without these samples by setting {.arg omitLowAbundSamples}={.val TRUE}."
        ),
        class = "normalyzerde_warning",
        call = NULL
      )
    }
  }

  if (length(notPassingThreshold) > 0) {
    naSamplesOmittedDf <- dataMatrix[, -sampleIndices[notPassingThreshold]]
  } else {
    dataMatrix
  }
}


#' Check whether all samples have replicates
#'
#' @param dataMatrix Prepared matrix containing expression data.
#' @param groups Vector containing condition levels
#' @param requireReplicates By default stops processing if not all samples
#'  have replicates
#' @return None
#' @keywords internal
validateSampleReplication <- function(
  dataMatrix,
  groups,
  requireReplicates = TRUE,
  quiet = FALSE
) {
  headerCounts <- table(groups)
  nonReplicatedSamples <- names(headerCounts[headerCounts == 1])

  if (length(nonReplicatedSamples) > 0) {
    if (requireReplicates) {
      cli::cli_abort(
        c(
          "Some group conditions have no replicates.",
          x = "Group conditions without replicates: {.val {nonReplicatedSamples}}",
          i = "Set {.arg requireReplicates}={.val FALSE} to continue, but replicate-dependent metrics and downstream steps may be unavailable."
        ),
        class = "normalyzerde_error",
        call = NULL
      )
    } else if (!quiet) {
      cli::cli_warn(
        c(
          "Some group conditions have no replicates.",
          "!" = "Group conditions without replicates: {.val {nonReplicatedSamples}}",
          i = "Continuing because {.arg requireReplicates}={.val FALSE}; replicate-dependent metrics and downstream steps may be unavailable."
        ),
        class = "normalyzerde_warning",
        call = NULL
      )
    }
  } else {
    if (!quiet) {
      cli::cli_inform(
        c(v = "Sample replication check: All samples have replicates")
      )
    }
  }
}


#' Check whether more than one sample is present
#'
#' @param dataMatrix Prepared dataframe.
#' @param groups Vector containing condition levels
#' @param requireReplicates By default stops processing if not all samples
#'  have replicates
#' @return None
#' @keywords internal
verifyMultipleSamplesPresent <- function(
  dataMatrix,
  groups,
  requireReplicates = TRUE,
  quiet = FALSE
) {
  samples <- groups[as.numeric(as.factor(groups)) > 0]
  distinctSamples <- unique(samples)

  if (length(samples) < 2) {
    cli::cli_abort(
      c(
        "At least two samples are required to run Normalyzer.",
        x = "Found {length(samples)} sample(s): {.val {samples}}."
      ),
      class = "normalyzerde_error",
      call = NULL
    )
  }

  if (length(distinctSamples) == 1) {
    if (requireReplicates) {
      cli::cli_abort(
        c(
          "Less than two distinct sample groups found.",
          x = "Found group: {.val {distinctSamples}}.",
          i = "For full processing, two or more sample groups are required.",
          i = "Set {.arg requireReplicates}={.val FALSE} to continue, but condition-comparison steps will be unavailable."
        ),
        class = "normalyzerde_error",
        call = NULL
      )
    } else if (!quiet) {
      cli::cli_warn(
        c(
          "Less than two distinct sample groups found.",
          "!" = "Found group: {.val {distinctSamples}}.",
          i = "Continuing because {.arg requireReplicates}={.val FALSE}; condition-comparison steps will be unavailable."
        ),
        class = "normalyzerde_warning",
        call = NULL
      )
    }
  } else if (length(distinctSamples) == 0) {
    cli::cli_abort(
      "No replicate groups found. Double check your input file and that your data haven't been filtered out in preceding input validation steps.",
      class = "normalyzerde_error",
      call = NULL
    )
  } else {
    if (!quiet) {
      cli::cli_inform(c(v = "Sample check: More than one sample group found"))
    }
  }
}
