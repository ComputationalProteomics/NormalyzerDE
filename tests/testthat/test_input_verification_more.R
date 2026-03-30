context("inputVerification (more)")

test_that("loadData supports MaxQuant input formats", {
  pepPath <- system.file(
    "extdata",
    "mq_peptides_100.txt",
    package = "NormalyzerDE"
  )
  protPath <- system.file(
    "extdata",
    "mq_proteinGroups_100.txt",
    package = "NormalyzerDE"
  )
  expect_true(nzchar(pepPath))
  expect_true(nzchar(protPath))

  pepMat <- loadData(pepPath, inputFormat = "maxquantpep")
  expect_true(is.matrix(pepMat))
  expect_true(nrow(pepMat) > 1)
  expect_true(ncol(pepMat) > 1)

  protMat <- loadData(protPath, inputFormat = "maxquantprot")
  expect_true(is.matrix(protMat))
  expect_true(nrow(protMat) > 1)
  expect_true(ncol(protMat) > 1)
})

test_that("loadData supports Proteios input format", {
  proteiosPath <- system.file(
    "extdata",
    "tiny_data_proteios.tsv",
    package = "NormalyzerDE"
  )
  expect_true(nzchar(proteiosPath))

  protMat <- loadData(proteiosPath, inputFormat = "proteios")
  expect_true(is.matrix(protMat))
  expect_true(nrow(protMat) > 0)
  expect_true(ncol(protMat) > 0)
})

test_that("loadData errors for unknown input formats", {
  expect_error(loadData("dummy", inputFormat = "bad"), "Unknown inputFormat")
})

test_that("loadDesign errors when sample/group columns are missing", {
  fp <- withr::local_tempfile(pattern = "design_", fileext = ".tsv")

  utils::write.table(
    data.frame(a = 1, b = 2),
    file = fp,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  expect_error(
    loadDesign(fp, sampleCol = "sample", groupCol = "group"),
    "Both sampleCol value and groupCol value must be present"
  )
})

test_that("verifyValidNumbers errors for below-one values when log transforming", {
  mat <- matrix(c("2", "0.5"), nrow = 1)
  expect_error(
    NormalyzerDE:::verifyValidNumbers(
      mat,
      groups = c("A", "B"),
      noLogTransform = FALSE,
      quiet = TRUE
    ),
    "below-one values"
  )
})

test_that("verifyValidNumbers warns when data looks already log2", {
  mat <- matrix(c("10", "12"), nrow = 1)
  expect_warning(
    NormalyzerDE:::verifyValidNumbers(
      mat,
      groups = c("A", "B"),
      noLogTransform = FALSE,
      quiet = TRUE
    ),
    "log2|noLogTransform=TRUE"
  )
})

test_that("verifyValidNumbers is silent for large linear-scale values", {
  mat <- matrix(c("1000", "2000"), nrow = 1)
  expect_silent(
    NormalyzerDE:::verifyValidNumbers(
      mat,
      groups = c("A", "B"),
      noLogTransform = FALSE,
      quiet = TRUE
    )
  )
})

test_that("verifyValidNumbers is silent when noLogTransform=TRUE", {
  mat <- matrix(c("10", "12"), nrow = 1)
  expect_silent(
    NormalyzerDE:::verifyValidNumbers(
      mat,
      groups = c("A", "B"),
      noLogTransform = TRUE,
      quiet = TRUE
    )
  )
})

test_that("verifyValidNumbers allows signed log2-scale values when noLogTransform=TRUE", {
  mat <- matrix(c("-1", "2"), nrow = 1)
  expect_silent(
    NormalyzerDE:::verifyValidNumbers(
      mat,
      groups = c("A", "B"),
      noLogTransform = TRUE,
      quiet = TRUE
    )
  )
})

test_that("getLowCountSampleFiltered can omit low-count samples", {
  mat <- matrix(c(1, 3, NA, 2, 4, NA), nrow = 2, byrow = TRUE)
  colnames(mat) <- c("s1", "s2", "s3")
  out <- expect_warning(
    NormalyzerDE:::getLowCountSampleFiltered(
      mat,
      groups = c("A", "B", "C"),
      threshold = 2,
      stopIfTooFew = FALSE
    ),
    "does not contain enough"
  )

  expect_equal(colnames(out), c("s1", "s2"))
})

test_that("getLowCountSampleFiltered errors when all samples fail threshold", {
  mat <- matrix(NA_real_, nrow = 2, ncol = 2)
  colnames(mat) <- c("s1", "s2")
  expect_error(
    NormalyzerDE:::getLowCountSampleFiltered(
      mat,
      groups = c("A", "B"),
      threshold = 1,
      stopIfTooFew = TRUE
    ),
    "None of the samples had enough"
  )
})

test_that("verifyDesignMatrix errors for missing/mismatched/duplicate samples", {
  full <- data.frame(
    feature = c("f1", "f2"),
    S1 = c(1, 2),
    S2 = c(3, 4),
    check.names = FALSE
  )

  expect_error(
    NormalyzerDE:::verifyDesignMatrix(
      full,
      data.frame(group = c("A", "B"), stringsAsFactors = FALSE),
      sampleCol = "sample"
    ),
    "Design matrix header must contain sampleCol name"
  )

  expect_error(
    NormalyzerDE:::verifyDesignMatrix(
      full,
      data.frame(
        sample = c("S1", "S3"),
        group = c("A", "B"),
        stringsAsFactors = FALSE
      ),
      sampleCol = "sample"
    ),
    "Not all columns present in design matrix are present in data matrix"
  )

  expect_error(
    NormalyzerDE:::verifyDesignMatrix(
      full,
      data.frame(
        sample = c("S1", "S1"),
        group = c("A", "B"),
        stringsAsFactors = FALSE
      ),
      sampleCol = "sample"
    ),
    "Sample labels must be unique"
  )
})

test_that("preprocessData replaces 0/empty/null and emits messages", {
  mat <- matrix(c("0", "", "null", "1"), nrow = 2, byrow = TRUE)
  msgs <- testthat::capture_messages(
    out <- NormalyzerDE:::preprocessData(mat, quiet = FALSE)
  )

  expect_true(any(grepl("fields with '0'", msgs)))
  expect_true(any(grepl("empty fields were replaced", msgs)))
  expect_true(any(grepl("'null' fields were replaced", msgs)))

  expect_true(all(is.na(out[1, ])))
  expect_equal(out[2, 2], "1")
})

test_that("loadRawDataFromFile errors for missing file and for parse warnings", {
  expect_error(
    NormalyzerDE:::loadRawDataFromFile(file.path(tempdir(), "no_such_file.tsv")),
    "An issue was encountered when attempting to load"
  )

  expect_error(
    NormalyzerDE:::loadRawDataFromFile(NA_character_),
    "Failed to read input file"
  )

  fp <- withr::local_tempfile(pattern = "embedded_nulls_", fileext = ".tsv")
  con <- file(fp, open = "wb")
  bytes <- c(charToRaw("A\tB\n1\t2"), as.raw(0), charToRaw("\n"))
  writeBin(bytes, con)
  close(con)

  expect_error(
    NormalyzerDE:::loadRawDataFromFile(fp),
    "embedded nulls"
  )
})

test_that("verifyMultipleSamplesPresent errors/warns/messages appropriately", {
  mat <- matrix(1, nrow = 1, ncol = 1)

  expect_error(
    NormalyzerDE:::verifyMultipleSamplesPresent(
      mat,
      groups = "A",
      requireReplicates = TRUE,
      quiet = TRUE
    ),
    "At least two samples are required"
  )

  expect_error(
    NormalyzerDE:::verifyMultipleSamplesPresent(
      mat,
      groups = c("A", "A"),
      requireReplicates = TRUE,
      quiet = TRUE
    ),
    "Found less than two distinct sample groups"
  )

  expect_warning(
    NormalyzerDE:::verifyMultipleSamplesPresent(
      mat,
      groups = c("A", "A"),
      requireReplicates = FALSE,
      quiet = FALSE
    ),
    "Found less than two distinct sample groups"
  )

  expect_message(
    NormalyzerDE:::verifyMultipleSamplesPresent(
      mat,
      groups = c("A", "B"),
      requireReplicates = TRUE,
      quiet = FALSE
    ),
    "More than one sample group found"
  )
})
