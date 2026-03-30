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
  expect_error(diannInputOptions(level = "bad"))
  expect_error(diannInputOptions(sep = ""), class = "normalyzerde_error")
  expect_error(
    diannInputOptions(minPositive = -1),
    class = "normalyzerde_error"
  )
  expect_error(
    diannInputOptions(rt = list(column = "RT")),
    class = "normalyzerde_error"
  )
  expect_error(diannInputOptions(qCutoff = 0.01), class = "normalyzerde_error")
})

test_that("diannInputOptions validates edge cases for columns/filters/rt", {
  expect_error(
    diannInputOptions(sampleCol = c("Run", "File.Name")),
    class = "normalyzerde_error"
  )
  expect_error(
    diannInputOptions(extraCols = c("Protein.Group", "")),
    class = "normalyzerde_error"
  )
  expect_error(diannInputOptions(decoy = NA), class = "normalyzerde_error")
  expect_error(diannInputOptions(qEnable = NA), class = "normalyzerde_error")
  expect_error(
    diannInputOptions(qCols = c("Q.Value", "")),
    class = "normalyzerde_error"
  )
  expect_error(
    diannInputOptions(qCutoffs = "0.01"),
    class = "normalyzerde_error"
  )
  expect_error(
    diannInputOptions(qCutoffs = c(0.01, NA_real_)),
    class = "normalyzerde_error"
  )
  expect_error(
    diannInputOptions(qCutoffs = -0.01),
    class = "normalyzerde_error"
  )
  expect_error(
    diannInputOptions(qCutoffs = 1.01),
    class = "normalyzerde_error"
  )

  expect_error(
    diannInputOptions(rt = c(TRUE, FALSE)),
    class = "normalyzerde_error"
  )
  expect_error(
    diannInputOptions(rt = list(col = "")),
    class = "normalyzerde_error"
  )
  expect_error(diannInputOptions(rt = 1), class = "normalyzerde_error")
})

test_that("diannNormalizeInputOptions validates nested list structure", {
  expect_error(
    NormalyzerDE:::diannNormalizeInputOptions("not_a_list"),
    class = "normalyzerde_error"
  )

  expect_error(
    NormalyzerDE:::diannNormalizeInputOptions(list(columns = "oops")),
    class = "normalyzerde_error"
  )

  expect_error(
    NormalyzerDE:::diannNormalizeInputOptions(list(filters = "oops")),
    class = "normalyzerde_error"
  )

  expect_error(
    NormalyzerDE:::diannNormalizeInputOptions(list(filters = list(q = 1))),
    class = "normalyzerde_error"
  )
  expect_error(
    NormalyzerDE:::diannNormalizeInputOptions(
      list(filters = list(decoy = "yes"))
    ),
    class = "normalyzerde_error"
  )
  expect_error(
    NormalyzerDE:::diannNormalizeInputOptions(
      list(filters = list(q = list(enable = "yes")))
    ),
    class = "normalyzerde_error"
  )
  expect_error(
    NormalyzerDE:::diannNormalizeInputOptions(
      list(filters = list(q = list(cols = c("Q.Value", ""))))
    ),
    class = "normalyzerde_error"
  )
  expect_error(
    NormalyzerDE:::diannNormalizeInputOptions(
      list(filters = list(q = list(cutoffs = c(0.01, NA_real_))))
    ),
    class = "normalyzerde_error"
  )
  expect_error(
    NormalyzerDE:::diannNormalizeInputOptions(
      list(filters = list(q = list(cutoffs = -0.01)))
    ),
    class = "normalyzerde_error"
  )

  expect_error(
    NormalyzerDE:::diannNormalizeInputOptions(list(rt = 1)),
    class = "normalyzerde_error"
  )
})

test_that("diannPrecursor preference keeps explicit auto but upgrades implicit defaults", {
  expect_equal(
    NormalyzerDE:::diannPreferPrecursorInputOptionsForLimpa(
      diannInputOptions()
    )$level,
    "precursor"
  )

  expect_equal(
    NormalyzerDE:::diannPreferPrecursorInputOptionsForLimpa(
      diannInputOptions(level = "auto")
    )$level,
    "auto"
  )

  expect_equal(
    NormalyzerDE:::diannPreferPrecursorInputOptionsForLimpa(
      list(level = "auto")
    )$level,
    "auto"
  )
})

test_that("diannInputOptions works end-to-end with DIANN q filtering", {
  tmpDir <- withr::local_tempdir(pattern = "diann_input_opts_")

  diannReport <- data.frame(
    Run = c("S1", "S2", "S1", "S2"),
    `Precursor.Id` = c("pep1", "pep1", "pep2", "pep2"),
    `Protein.Group` = c("P1", "P1", "P2", "P2"),
    `Precursor.Quantity` = c(100, 200, 300, 400),
    `Q.Value` = c(0.005, 0.005, 0.02, 0.02),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  paths <- nd_write_data_and_design(
    tmp_dir = tmpDir,
    data = diannReport,
    data_name = "diann_report_qfilter_opts.tsv",
    design_name = "diann_design_qfilter_opts.tsv"
  )

  opts <- diannInputOptions(
    level = "precursor",
    quantityCol = "Precursor.Quantity",
    qCols = "Q.Value",
    qCutoffs = 0.01,
    rt = FALSE
  )

  se <- setupRawContrastObject(
    dataPath = paths$dataPath,
    designPath = paths$designPath,
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
