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

test_that("diannInputOptions validates edge cases for columns/filters/rt", {
  expect_error(diannInputOptions(sampleCol = c("Run", "File.Name")), "sampleCol")
  expect_error(diannInputOptions(extraCols = c("Protein.Group", "")), "extraCols")
  expect_error(diannInputOptions(decoy = NA), "decoy")
  expect_error(diannInputOptions(qEnable = NA), "qEnable")
  expect_error(diannInputOptions(qCols = c("Q.Value", "")), "qCols")
  expect_error(diannInputOptions(qCutoffs = "0.01"), "qCutoffs must be numeric")

  expect_error(diannInputOptions(rt = c(TRUE, FALSE)), "rt must be a single")
  expect_error(diannInputOptions(rt = list(col = "")), "rt\\$col")
  expect_error(diannInputOptions(rt = 1), "rt must be NULL")
})

test_that("diannNormalizeInputOptions validates nested list structure", {
  expect_error(
    NormalyzerDE:::diannNormalizeInputOptions("not_a_list"),
    "inputOptions must be a list"
  )

  expect_error(
    NormalyzerDE:::diannNormalizeInputOptions(list(columns = "oops")),
    "inputOptions\\$columns must be a list"
  )

  expect_error(
    NormalyzerDE:::diannNormalizeInputOptions(list(filters = "oops")),
    "inputOptions\\$filters must be a list"
  )

  expect_error(
    NormalyzerDE:::diannNormalizeInputOptions(list(filters = list(q = 1))),
    "inputOptions\\$filters\\$q must be logical"
  )

  expect_error(
    NormalyzerDE:::diannNormalizeInputOptions(list(rt = 1)),
    "inputOptions\\$rt must be a list"
  )
})

test_that("diannInputOptions works end-to-end with DIANN q filtering", {
  tmpDir <- withr::local_tempdir(pattern = "diann_input_opts_")
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

  nd_write_table(diannReport, dataPath)
  nd_write_table(nd_two_sample_design(), designPath)

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
