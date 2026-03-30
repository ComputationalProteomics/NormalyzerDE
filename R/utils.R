#' Create empty directory for run
#'
#' Creates a directory at provided path named to the jobname.
#'
#' @param jobName Name of the run.
#' @param outputDir Path to directory where to create the output directory.
#' @param reuseOutputDir Reuse an existing non-empty output directory.
#' @return Path to newly created directory.
#' @export
#' @examples
#' setupJobDir("job_name", "path/to/outdir")
setupJobDir <- function(jobName, outputDir, reuseOutputDir = FALSE) {
  sanitizedJobName <- sanitizeJobName(jobName)

  if (is.null(outputDir)) {
    jobDir <- file.path(getwd(), sanitizedJobName)
  } else {
    jobDir <- file.path(outputDir, sanitizedJobName)
  }

  if (
    dir.exists(jobDir) &&
      length(list.files(jobDir, all.files = TRUE, no.. = TRUE)) > 0
  ) {
    if (!isTRUE(reuseOutputDir)) {
      cli::cli_abort(
        c(
          "Output directory already exists and is not empty: {.path {jobDir}}.",
          i = "Choose a new {.arg jobName} or {.arg outputDir}, or set {.arg reuseOutputDir}={.val TRUE} to reuse it explicitly."
        ),
        class = "normalyzerde_error",
        call = NULL
      )
    }

    cli::cli_warn(
      c(
        "Reusing existing output directory: {.path {jobDir}}.",
        i = "Existing files may be mixed with new outputs."
      ),
      class = "normalyzerde_warning",
      call = NULL
    )
  }
  createDirectory(jobDir)

  jobDir
}

sanitizeJobName <- function(jobName) {
  if (length(jobName) < 1 || is.null(jobName) || is.na(jobName[1])) {
    cli::cli_abort(
      "Invalid {.arg jobName}: must be a non-empty character value.",
      class = "normalyzerde_error",
      call = NULL
    )
  }

  originalJobName <- as.character(jobName[1])
  sanitizedJobName <- basename(originalJobName)

  sanitizedJobName <- gsub("[[:cntrl:]]", "_", sanitizedJobName)
  sanitizedJobName <- gsub("[<>:\"/\\\\|?*]", "_", sanitizedJobName)
  sanitizedJobName <- trimws(sanitizedJobName)
  sanitizedJobName <- sub("[. ]+$", "", sanitizedJobName)

  if (!nzchar(sanitizedJobName)) {
    cli::cli_abort(
      "Invalid {.arg jobName}: must contain at least one character after sanitization.",
      class = "normalyzerde_error",
      call = NULL
    )
  }

  if (!identical(originalJobName, sanitizedJobName)) {
    cli::cli_warn(
      "{.arg jobName} was sanitized from {.val {originalJobName}} to {.val {sanitizedJobName}} for filesystem compatibility.",
      class = "normalyzerde_warning",
      call = NULL
    )
  }

  sanitizedJobName
}

#' Get dataframe with raw data column sorted on replicates
#'
#' @param rawDataOnly Data frame with unparsed input data matrix.
#' @param groups Vector containing condition levels.
#' @return rawData sorted on replicate
#' @keywords internal
getReplicateSortedData <- function(rawDataOnly, groups) {
  indexList <- getIndexList(groups)
  orderedIndices <- unlist(indexList[sort(names(indexList))])

  rawDataOnly[, orderedIndices, drop = FALSE]
}

#' Create directory, or return error if already present
#'
#' @param targetPath Path where to attempt to create directory
#' @return None
#' @keywords internal
createDirectory <- function(targetPath) {
  if (dir.exists(targetPath)) {
    return(invisible(NULL))
  }

  if (file.exists(targetPath)) {
    cli::cli_abort(
      "Path already exists and is not a directory: {.path {targetPath}}.",
      class = "normalyzerde_error",
      call = NULL
    )
  }

  created <- dir.create(targetPath, recursive = TRUE, showWarnings = FALSE)
  if (!isTRUE(created) && !dir.exists(targetPath)) {
    cli::cli_abort(
      "Failed to create directory: {.path {targetPath}}.",
      class = "normalyzerde_error",
      call = NULL
    )
  }

  invisible(NULL)
}

#' Get number of seconds between two Sys.time() objects
#'
#' @param start Start-time object
#' @param end End-time object
#' @return None
#' @keywords internal
elapsedSecondsBetweenSystimes <- function(start, end) {
  startSecond <- strtoi(format(start, "%s"))
  endSecond <- strtoi(format(end, "%s"))
  elapsed <- end - start
  elapsed
}

#' Return list containing vector positions of values in string
#'
#' @param targetVector Vector of values to index (for example, condition labels).
#' @return indexList List where key is condition level and values are indices
#'   for the condition
#' @keywords internal
getIndexList <- function(targetVector) {
  indexList <- list()
  uniqVals <- unique(targetVector)
  for (val in uniqVals) {
    indexList[[toString(val)]] <- which(targetVector == val)
  }
  indexList
}

#' Get contrast vector (TRUE/FALSE-values) indicating whether both at least
#' half values are present, and each sample has at least one non-NA value
#'
#' @param dataMatrix Matrix with expression values for entities in replicate
#'  samples.
#' @param replicateHeader Header showing how samples in matrix are replicated.
#' @param minCount Minimum number of required values present in samples.
#' @return Contrast vector
#' @keywords internal
getRowNAFilterContrast <- function(dataMatrix, replicateHeader, minCount = 1) {
  replicatesHaveData <- rep(TRUE, nrow(dataMatrix))
  indexList <- getIndexList(replicateHeader)

  for (sampleIndex in seq_along(names(indexList))) {
    repVal <- names(indexList)[sampleIndex]
    cols <- indexList[[repVal]]

    nbrNAperReplicate <- rowSums(is.na(dataMatrix[, cols, drop = FALSE]))
    nbrReplicates <- length(cols)
    nbrNonNA <- nbrReplicates - nbrNAperReplicate
    replicatesHaveData <- (nbrNonNA >= minCount & replicatesHaveData)
  }

  replicatesHaveData
}

#' Generate a random test dataset with features, sample values and retention times
#'
#' @param nSamples Number of samples
#' @param nFeatures Number of features
#' @param rtMin Minimum retention time
#' @param rtMax Maximum retention time
#' @param mean Mean value for sample intensities
#' @param sd Standard deviation for sample intensities
#' @return Test dataset
#' @export
#' @examples
#' df <- setupTestData(6, 20)
#' df <- setupTestData(6, 20, mean=15, sd=1)
#' @keywords internal
setupTestData <- function(
  nSamples,
  nFeatures,
  rtMin = 40,
  rtMax = 80,
  mean = 20,
  sd = 4
) {
  featureNames <- paste0("feature_", seq(1, nFeatures))
  sampleData <- matrix(
    stats::rnorm(nSamples * nFeatures, mean, sd),
    nFeatures,
    nSamples
  )
  rtData <- stats::runif(nFeatures, rtMin, rtMax)

  df <- data.frame(
    feature = featureNames,
    RT = rtData,
    as.data.frame(sampleData)
  )

  colnames(df) <- c("feature", "RT", paste0("S", seq(1, nSamples)))
  df
}

#' General function for calculating percentage difference of average column
#' means in matrix
#'
#' @param targetMat Matrix for which column means should be compared
#' @return percDiffVector Vector with percentage difference, where first element
#'   always will be 100
#' @keywords internal
calculatePercentageAvgDiffInMat <- function(targetMat) {
  calculatePercDiff <- function(sampleIndex, mat) {
    mean(mat[, sampleIndex]) * 100 / mean(mat[, 1])
  }

  percDiffVector <- vapply(
    seq_len(ncol(targetMat)),
    calculatePercDiff,
    0,
    mat = targetMat
  )

  percDiffVector
}

#' Filter rows with lower than given number of replicates for any condition
#'
#' @param df Data frame with expression data to filter
#' @param groups Condition groups header
#' @param leastRep Minimum number of replicates in each group
#'   to retain
#' @return collDesignDf Reduced design matrix
#' @keywords internal
filterLowRep <- function(df, groups, leastRep = 2) {
  allReplicatesHaveValuesContrast <- function(row, groups, minCount) {
    names(row) <- groups
    repCounts <- table(names(stats::na.omit(row)))
    length(repCounts) == length(unique(groups)) &&
      min(repCounts) >= minCount ||
      minCount == 0
  }

  rowMeetThresContrast <- apply(
    df,
    1,
    allReplicatesHaveValuesContrast,
    groups = groups,
    minCount = leastRep
  )

  filteredDf <- df[rowMeetThresContrast, , drop = FALSE]
  filteredDf
}

#' Impute missing groups when another group has sufficient observations
#'
#' Imputes one low value in groups containing only missing values when at least
#' one group for the same feature has a sufficient fraction of observed values.
#'
#' @param df Data frame with expression data
#' @param groups Condition-group labels
#' @param minFraction Minimum observed fraction required in one group before
#'   imputing a low value in other groups
#' @return imputedDf Imputed data frame
#' @keywords internal
imputeGroupValues <- function(df, groups, minFraction = 0.75) {
  if (all(is.na(df))) {
    return(df)
  }

  groups <- as.character(groups)
  minValue <- min(df, na.rm = TRUE)
  groupIndices <- split(seq_along(groups), groups)
  groupTotalCounts <- lengths(groupIndices)
  inputRowNames <- if (is.data.frame(df)) {
    attr(df, "row.names")
  } else {
    rownames(df)
  }

  imputeRow <- function(
    row,
    groupIndices,
    groupTotalCounts,
    minFraction,
    minVal
  ) {
    groupNonNaCounts <- vapply(
      groupIndices,
      function(indices) sum(!is.na(row[indices])),
      integer(1)
    )

    if (
      any(groupNonNaCounts / groupTotalCounts >= minFraction) &&
        any(groupNonNaCounts == 0)
    ) {
      missingGroups <- names(groupNonNaCounts)[groupNonNaCounts == 0]
      for (groupLabel in missingGroups) {
        row[groupIndices[[groupLabel]][1]] <- minVal
      }
    }

    row
  }

  imputed <- apply(
    df,
    1,
    imputeRow,
    groupIndices = groupIndices,
    groupTotalCounts = groupTotalCounts,
    minFraction = minFraction,
    minVal = minValue
  )

  imputedMat <- t(imputed)
  colnames(imputedMat) <- colnames(df)
  if (!is.null(inputRowNames)) {
    rownames(imputedMat) <- rownames(df)
  }

  if (is.data.frame(df)) {
    imputedDf <- as.data.frame(imputedMat)
    if (!is.null(inputRowNames)) {
      attr(imputedDf, "row.names") <- inputRowNames
    }
    return(imputedDf)
  }
  imputedMat
}
