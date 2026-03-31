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
  expect_error(
    loadData("dummy", inputFormat = "bad"),
    class = "normalyzerde_error"
  )
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
    class = "normalyzerde_error"
  )
})

test_that("loadDesign coerces numeric-like sample IDs and groups for user files", {
  fp <- withr::local_tempfile(pattern = "design_numeric_", fileext = ".tsv")

  utils::write.table(
    data.frame(
      `Sample ID` = c(101, 102),
      `Group ID` = c(1, 1),
      check.names = FALSE
    ),
    file = fp,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  design <- loadDesign(fp, sampleCol = "Sample ID", groupCol = "Group ID")

  expect_type(design[["Sample ID"]], "character")
  expect_equal(design[["Sample ID"]], c("101", "102"))
  expect_true(is.factor(design[["Group ID"]]))
  expect_equal(as.character(design[["Group ID"]]), c("1", "1"))
})

test_that("loadDesign drops blank Excel-export columns made only of tabs", {
  fp <- withr::local_tempfile(pattern = "design_excel_", fileext = ".tsv")
  writeLines(
    c(
      "sample\tgroup\tbatch\t\t",
      "S1\tA\tb1\t\t",
      "S2\tB\tb2\t\t"
    ),
    fp
  )

  design <- loadDesign(fp, sampleCol = "sample", groupCol = "group")

  expect_equal(colnames(design), c("sample", "group", "batch"))
  expect_equal(design$batch, c("b1", "b2"))
})

test_that("setupRawDataObject follows design order and treats non-design columns as annotation", {
  tmpDir <- withr::local_tempdir(pattern = "setup_raw_data_reorder_")
  dataPath <- file.path(tmpDir, "reordered_data.tsv")
  designPath <- file.path(tmpDir, "reordered_design.tsv")

  rawDf <- data.frame(
    feature = c("f1", "f2"),
    note = c("pep1", "pep2"),
    s2 = c(20, 40),
    s1 = c(10, 30),
    extra = c("x", "y"),
    check.names = FALSE
  )
  designDf <- data.frame(
    sample = c("s1", "s2"),
    group = c("A", "B"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  nd_write_table(rawDf, dataPath)
  nd_write_table(designDf, designPath)

  se <- setupRawDataObject(
    dataPath = dataPath,
    designPath = designPath,
    inputFormat = "default"
  )

  assayMat <- SummarizedExperiment::assay(se)
  expect_equal(colnames(assayMat), c("s1", "s2"))
  expect_equal(unname(as.numeric(assayMat[1, ])), c(10, 20))
  expect_true(all(
    c("feature", "note", "extra") %in%
      colnames(SummarizedExperiment::rowData(se))
  ))
})

test_that("setupRawContrastObject follows design order for user-supplied matrices", {
  tmpDir <- withr::local_tempdir(pattern = "setup_raw_contrast_reorder_")
  dataPath <- file.path(tmpDir, "reordered_contrast.tsv")
  designPath <- file.path(tmpDir, "reordered_design.tsv")

  fullDf <- data.frame(
    feature = c("f1", "f2"),
    note = c("pep1", "pep2"),
    s2 = c(20, 40),
    s1 = c(10, 30),
    extra = c("x", "y"),
    check.names = FALSE
  )
  designDf <- data.frame(
    sample = c("s1", "s2"),
    group = c("A", "B"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  nd_write_table(fullDf, dataPath)
  nd_write_table(designDf, designPath)

  se <- setupRawContrastObject(
    dataPath = dataPath,
    designPath = designPath,
    sampleColName = "sample"
  )

  assayMat <- SummarizedExperiment::assay(se)
  expect_equal(colnames(assayMat), c("s1", "s2"))
  expect_equal(unname(as.numeric(assayMat[1, ])), c(10, 20))
  expect_true(all(
    c("feature", "note", "extra") %in%
      colnames(SummarizedExperiment::rowData(se))
  ))
})

test_that("setupRawContrastObject ignores blank Excel-export design columns", {
  tmpDir <- withr::local_tempdir(pattern = "setup_raw_contrast_excel_")
  dataPath <- file.path(tmpDir, "contrast.tsv")
  designPath <- file.path(tmpDir, "design.tsv")

  fullDf <- data.frame(
    feature = c("f1", "f2"),
    s1 = c(10, 30),
    s2 = c(20, 40),
    check.names = FALSE
  )
  nd_write_table(fullDf, dataPath)

  writeLines(
    c(
      "sample\tgroup\tbatch\t\t",
      "s1\tA\tb1\t\t",
      "s2\tB\tb2\t\t"
    ),
    designPath
  )

  se <- setupRawContrastObject(
    dataPath = dataPath,
    designPath = designPath,
    sampleColName = "sample"
  )

  expect_equal(
    colnames(as.data.frame(SummarizedExperiment::colData(se))),
    c("sample", "group", "batch")
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
    class = "normalyzerde_error"
  )
})

test_that("verifyValidNumbers errors for below-one scientific notation", {
  mat <- matrix(c("2", "1e-3"), nrow = 1)
  expect_error(
    NormalyzerDE:::verifyValidNumbers(
      mat,
      groups = c("A", "B"),
      noLogTransform = FALSE,
      quiet = TRUE
    ),
    class = "normalyzerde_error"
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
    class = "normalyzerde_warning"
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

test_that("verifyValidNumbers rejects comma decimals from locale-specific input", {
  mat <- matrix(c("1,23", "4"), nrow = 1)
  expect_error(
    NormalyzerDE:::verifyValidNumbers(
      mat,
      groups = c("A", "B"),
      noLogTransform = TRUE,
      quiet = TRUE
    ),
    class = "normalyzerde_error"
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
    class = "normalyzerde_warning"
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
    class = "normalyzerde_error"
  )
})

test_that("verifyDesignMatrix errors for missing/mismatched/duplicate samples", {
  full <- data.frame(
    feature = c("f1", "f2"),
    S1 = c(1, 2),
    S2 = c(3, 4),
    check.names = FALSE
  )

  nd_expect_error_cases(list(
    "missing sample column" = function() {
      NormalyzerDE:::verifyDesignMatrix(
        full,
        data.frame(group = c("A", "B"), stringsAsFactors = FALSE),
        sampleCol = "sample"
      )
    },
    "mismatched sample ids" = function() {
      NormalyzerDE:::verifyDesignMatrix(
        full,
        data.frame(
          sample = c("S1", "S3"),
          group = c("A", "B"),
          stringsAsFactors = FALSE
        ),
        sampleCol = "sample"
      )
    },
    "duplicate sample ids" = function() {
      NormalyzerDE:::verifyDesignMatrix(
        full,
        data.frame(
          sample = c("S1", "S1"),
          group = c("A", "B"),
          stringsAsFactors = FALSE
        ),
        sampleCol = "sample"
      )
    }
  ))
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
    NormalyzerDE:::loadRawDataFromFile(file.path(
      tempdir(),
      "no_such_file.tsv"
    )),
    class = "normalyzerde_error"
  )

  expect_error(
    NormalyzerDE:::loadRawDataFromFile(NA_character_),
    class = "normalyzerde_error"
  )

  fp <- withr::local_tempfile(pattern = "embedded_nulls_", fileext = ".tsv")
  con <- file(fp, open = "wb")
  bytes <- c(charToRaw("A\tB\n1\t2"), as.raw(0), charToRaw("\n"))
  writeBin(bytes, con)
  close(con)

  expect_error(
    NormalyzerDE:::loadRawDataFromFile(fp),
    class = "normalyzerde_error"
  )
})

test_that("filterOnlyNARows drops fully missing features and preserves alignment", {
  se <- SummarizedExperiment::SummarizedExperiment(
    assay = matrix(
      c(
        1,
        2,
        NA,
        NA,
        3,
        4
      ),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(c("keep1", "drop", "keep2"), c("s1", "s2"))
    ),
    rowData = data.frame(
      feature = c("keep1", "drop", "keep2"),
      label = c("A", "B", "C"),
      stringsAsFactors = FALSE,
      row.names = c("keep1", "drop", "keep2"),
      check.names = FALSE
    ),
    colData = data.frame(
      sample = c("s1", "s2"),
      group = c("A", "B"),
      stringsAsFactors = FALSE
    )
  )

  msgs <- character()
  out <- withCallingHandlers(
    NormalyzerDE:::filterOnlyNARows(se),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )

  expect_equal(nrow(out), 2)
  expect_equal(rownames(SummarizedExperiment::assay(out)), c("keep1", "keep2"))
  expect_equal(
    as.character(SummarizedExperiment::rowData(out)$feature),
    c("keep1", "keep2")
  )
  expect_equal(
    as.character(SummarizedExperiment::rowData(out)$label),
    c("A", "C")
  )
  expect_true(any(grepl("entries with only NA values omitted", msgs)))
})

test_that("verifySummarizedExperiment errors when sample metadata does not match assay columns", {
  se <- SummarizedExperiment::SummarizedExperiment(
    assay = matrix(
      c(1, 2, 3, 4),
      nrow = 2,
      dimnames = list(c("f1", "f2"), c("s1", "s2"))
    ),
    colData = data.frame(
      sample = c("s1", "s3"),
      group = c("A", "B"),
      stringsAsFactors = FALSE
    ),
    rowData = data.frame(feature = c("f1", "f2"))
  )

  expect_error(
    NormalyzerDE:::verifySummarizedExperiment(se, sampleCol = "sample"),
    class = "normalyzerde_error"
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
    class = "normalyzerde_error"
  )

  expect_error(
    NormalyzerDE:::verifyMultipleSamplesPresent(
      mat,
      groups = c("A", "A"),
      requireReplicates = TRUE,
      quiet = TRUE
    ),
    class = "normalyzerde_error"
  )

  expect_warning(
    NormalyzerDE:::verifyMultipleSamplesPresent(
      mat,
      groups = c("A", "A"),
      requireReplicates = FALSE,
      quiet = FALSE
    ),
    class = "normalyzerde_warning"
  )

  expect_message(
    NormalyzerDE:::verifyMultipleSamplesPresent(
      mat,
      groups = c("A", "B"),
      requireReplicates = TRUE,
      quiet = FALSE
    )
  )
})

test_that("getVerifiedNormalyzerObject re-checks groups after omitting low-count samples", {
  mat <- matrix(
    c(
      NA,
      10,
      11,
      NA,
      12,
      13
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("f1", "f2"), c("s1", "s2", "s3"))
  )

  se <- SummarizedExperiment::SummarizedExperiment(
    assay = mat,
    colData = data.frame(
      sample = c("s1", "s2", "s3"),
      group = c("A", "B", "B"),
      stringsAsFactors = FALSE
    ),
    rowData = data.frame(feature = c("f1", "f2"), stringsAsFactors = FALSE)
  )
  S4Vectors::metadata(se) <- list(sample = "sample", group = "group")

  warnings <- character()
  nds <- withCallingHandlers(
    getVerifiedNormalyzerObject(
      jobName = "omit_low_count_groups",
      summarizedExp = se,
      threshold = 1,
      omitSamples = TRUE,
      requireReplicates = FALSE,
      quiet = FALSE,
      noLogTransform = TRUE,
      tinyRunThres = 50
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_equal(as.character(designMatrix(nds)$sample), c("s2", "s3"))
  expect_equal(colnames(filterrawdata(nds)), c("s2", "s3"))
  expect_true(any(grepl(
    "Less than two distinct sample groups found",
    warnings
  )))
  expect_false(any(grepl("Some group conditions have no replicates", warnings)))
})
