#' Create input options for Proteios matrices
#'
#' Helper to construct an \code{inputOptions} list for
#' \code{inputFormat = "proteios"}. Currently this supports selecting the input
#' delimiter.
#'
#' @param sep Field separator used when reading the data matrix.
#' @return A list suitable to pass as \code{inputOptions}.
#' @export
#' @examples
#' proteiosInputOptions(sep = ";")
proteiosInputOptions <- function(sep = "\t") {
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

#' Create input options for MaxQuant matrices
#'
#' Helper to construct an \code{inputOptions} list for \code{inputFormat =
#' "maxquantpep"} and \code{"maxquantprot"}. Currently this supports selecting
#' the input delimiter.
#'
#' @param sep Field separator used when reading the data matrix.
#' @return A list suitable to pass as \code{inputOptions}.
#' @export
#' @examples
#' maxQuantInputOptions(sep = "\t")
maxQuantInputOptions <- function(sep = "\t") {
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

proteiosToNormalyzer <- function(proteiosFp, sep = "\t") {
  valuesDf <- as.matrix(
    utils::read.table(
      proteiosFp,
      header = FALSE,
      sep = sep,
      stringsAsFactors = FALSE,
      quote = ""
    )
  )
  fullMat <- as.matrix(valuesDf)
  fullMat
}

maxQuantToNormalyzer <- function(maxQuantFp, protLevel, sep = "\t") {
  pepIntensityPattern <- "Intensity\\."

  if (!protLevel) {
    annotCols <- c(
      "Sequence",
      "Mass",
      "Proteins",
      "Leading.razor.protein",
      "PEP",
      "Charges"
    )
    matrixType <- "peptide.txt"
  } else {
    annotCols <- c("Protein.IDs", "Majority.protein.IDs", "Fasta.headers")
    matrixType <- "proteinGroups.txt"
  }

  fullDf <- utils::read.csv(
    maxQuantFp,
    sep = sep,
    stringsAsFactors = FALSE,
    comment.char = "",
    quote = "",
    header = TRUE
  )
  cnames <- colnames(fullDf)
  intensityCols <- cnames[grepl(pepIntensityPattern, cnames)]

  headerNames <- c(annotCols, intensityCols)
  headerNamesTrimmed <- gsub("Intensity.", "", headerNames)

  if (!(all(annotCols %in% colnames(fullDf)))) {
    cli::cli_abort(
      c(
        "Missing expected columns when processing MaxQuant {.val {matrixType}}.",
        i = "Expected columns (spaces instead of dots): {paste(annotCols, collapse = \", \")}.",
        i = "Columns in input data: {paste(colnames(fullDf), collapse = \", \")}."
      ),
      class = "normalyzerde_error",
      call = NULL
    )
  }

  valuesDf <- fullDf[, c(annotCols, intensityCols)]
  rawDf <- rbind(headerNamesTrimmed, valuesDf)
  rawMat <- as.matrix(rawDf)

  return(rawMat)
}
