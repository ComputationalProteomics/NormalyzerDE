context("diannInputOptions")

test_that("diannInputOptions returns a validated inputOptions list", {
  opts <- diannInputOptions(
    level = "precursor",
    quantityCol = "Precursor.Quantity",
    qCols = "Q.Value",
    qCutoffs = 0.01,
    minPositive = 0.1,
    rt = "RT"
  )

  expect_true(is.list(opts))
  expect_equal(opts$level, "precursor")
  expect_equal(opts$sep, "\t")
  expect_equal(opts$columns$quantity, "Precursor.Quantity")
  expect_true(isTRUE(opts$filters$decoy))
  expect_true(is.list(opts$filters$q))
  expect_equal(opts$filters$q$cols, "Q.Value")
  expect_equal(opts$filters$min_positive, 0.1)
})

test_that("diannInputOptions validates arguments and blocks typos", {
  expect_error(diannInputOptions(level = "bad"), "should be one of")
  expect_error(diannInputOptions(sep = ""), "sep")
  expect_error(diannInputOptions(minPositive = -1), "minPositive")
  expect_error(diannInputOptions(rt = list(column = "RT")), "rt")
  expect_error(diannInputOptions(qCutoff = 0.01), "Unknown argument")
})

test_that("diannInputOptions works end-to-end with DIANN q filtering", {
  tmpDir <- tempdir()
  dataPath <- file.path(tmpDir, "diann_report_qfilter_opts.tsv")
  designPath <- file.path(tmpDir, "diann_design_qfilter_opts.tsv")

  diannReport <- data.frame(
    Run = c("S1", "S2", "S1", "S2"),
    `Precursor.Id` = c("pep1", "pep1", "pep2", "pep2"),
    `Protein.Group` = c("P1", "P1", "P2", "P2"),
    `Precursor.Quantity` = c(100, 200, 300, 400),
    `Q.Value` = c(0.005, 0.005, 0.02, 0.02),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  design <- data.frame(
    sample = c("S1", "S2"),
    group = c("A", "B"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  utils::write.table(
    diannReport,
    file = dataPath,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  utils::write.table(
    design,
    file = designPath,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  opts <- diannInputOptions(
    level = "precursor",
    quantityCol = "Precursor.Quantity",
    qCols = "Q.Value",
    qCutoffs = 0.01,
    rt = FALSE
  )

  se <- setupRawContrastObject(
    dataPath = dataPath,
    designPath = designPath,
    sampleColName = "sample",
    inputFormat = "diann",
    inputOptions = opts
  )

  mat <- SummarizedExperiment::assay(se)
  ann <- as.data.frame(
    SummarizedExperiment::rowData(se),
    stringsAsFactors = FALSE
  )

  expect_true(nrow(mat) == 1)
  expect_true("Precursor.Id" %in% colnames(ann))
  expect_true(ann$Precursor.Id[1] == "pep1")
})
