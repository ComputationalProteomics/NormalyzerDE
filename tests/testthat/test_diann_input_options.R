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

test_that("diann helper functions normalize option branches and file headers", {
  headerPath <- withr::local_tempfile(pattern = "diann_header_", fileext = ".tsv")
  writeLines("Run\tQ.Value\tRT", headerPath)
  expect_equal(
    NormalyzerDE:::diannReadHeader(headerPath),
    c("Run", "Q.Value", "RT")
  )

  emptyPath <- withr::local_tempfile(pattern = "diann_empty_", fileext = ".tsv")
  file.create(emptyPath)
  expect_error(
    NormalyzerDE:::diannReadHeader(emptyPath),
    class = "normalyzerde_error"
  )

  expect_equal(
    NormalyzerDE:::diannSelectFirstPresent(c("Missing", "Q.Value"), c("Run", "Q.Value")),
    "Q.Value"
  )
  expect_null(
    NormalyzerDE:::diannSelectFirstPresent(c("Missing", "Absent"), c("Run", "Q.Value"))
  )

  expect_equal(NormalyzerDE:::diannResolveMinPositive(NULL), 0)
  expect_equal(NormalyzerDE:::diannResolveMinPositive("0.25"), 0.25)
  expect_error(
    NormalyzerDE:::diannResolveMinPositive(-1),
    class = "normalyzerde_error"
  )
  expect_error(
    NormalyzerDE:::diannResolveMinPositive("bad"),
    class = "normalyzerde_error"
  )

  expect_null(NormalyzerDE:::diannValidateQCutoffs(NULL))
  expect_equal(NormalyzerDE:::diannValidateQCutoffs(numeric()), numeric())
  expect_equal(NormalyzerDE:::diannValidateQCutoffs(c(0, 1)), c(0, 1))
  expect_error(
    NormalyzerDE:::diannValidateQCutoffs("0.1"),
    class = "normalyzerde_error"
  )
  expect_error(
    NormalyzerDE:::diannValidateQCutoffs(c(0.1, Inf)),
    class = "normalyzerde_error"
  )

  opts_true <- NormalyzerDE:::diannNormalizeInputOptions(
    list(
      columns = list(rt = "Aligned.RT"),
      filters = list(q = FALSE),
      rt = TRUE
    )
  )
  expect_false(opts_true$filterQValue)
  expect_equal(opts_true$rtCol, "Aligned.RT")

  opts_list <- NormalyzerDE:::diannNormalizeInputOptions(
    list(
      columns = list(rt = "Aligned.RT"),
      rt = list()
    )
  )
  expect_equal(opts_list$rtCol, "Aligned.RT")

  opts_list_default <- NormalyzerDE:::diannNormalizeInputOptions(list(rt = list()))
  expect_equal(opts_list_default$rtCol, "RT")

  opts_true_default <- NormalyzerDE:::diannNormalizeInputOptions(list(rt = TRUE))
  expect_equal(opts_true_default$rtCol, "RT")

  opts_char <- NormalyzerDE:::diannNormalizeInputOptions(list(rt = "Observed.RT"))
  expect_equal(opts_char$rtCol, "Observed.RT")

  expect_error(
    NormalyzerDE:::diannNormalizeInputOptions(list(filters = list(q = c(TRUE, FALSE)))),
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

test_that("diann precursor preference and default tracking cover early-return cases", {
  expect_equal(
    NormalyzerDE:::diannPreferPrecursorInputOptionsForLimpa(NULL)$level,
    "precursor"
  )

  marker <- "not_a_list"
  expect_identical(
    NormalyzerDE:::diannPreferPrecursorInputOptionsForLimpa(marker),
    marker
  )

  protein_opts <- list(level = "protein")
  expect_identical(
    NormalyzerDE:::diannPreferPrecursorInputOptionsForLimpa(protein_opts),
    protein_opts
  )

  bad_columns <- list(level = "auto", columns = "oops")
  expect_identical(
    NormalyzerDE:::diannPreferPrecursorInputOptionsForLimpa(bad_columns),
    bad_columns
  )

  feature_cols <- list(level = "auto", columns = list(feature = "Protein.Group"))
  expect_identical(
    NormalyzerDE:::diannPreferPrecursorInputOptionsForLimpa(feature_cols),
    feature_cols
  )

  quantity_cols <- list(
    level = "auto",
    columns = list(quantity = "Precursor.Quantity")
  )
  expect_identical(
    NormalyzerDE:::diannPreferPrecursorInputOptionsForLimpa(quantity_cols),
    quantity_cols
  )

  expect_false(
    NormalyzerDE:::diannWasDefaultedToPrecursorForLimpa(
      list(level = "precursor"),
      list(level = "precursor")
    )
  )
  expect_true(
    NormalyzerDE:::diannWasDefaultedToPrecursorForLimpa(
      list(level = "auto"),
      list(level = "precursor")
    )
  )
  expect_false(
    NormalyzerDE:::diannWasDefaultedToPrecursorForLimpa(
      list(level = "auto"),
      "precursor"
    )
  )

  expect_identical(
    NormalyzerDE:::diannPreferPrecursorInputOptionsForLimpa(
      list(level = "protein", ._normalyzerde_explicit_level = FALSE)
    ),
    list(level = "protein", ._normalyzerde_explicit_level = FALSE)
  )
  expect_identical(
    NormalyzerDE:::diannPreferPrecursorInputOptionsForLimpa(
      list(level = "auto", columns = "oops", ._normalyzerde_explicit_level = FALSE)
    ),
    list(level = "auto", columns = "oops", ._normalyzerde_explicit_level = FALSE)
  )
  expect_identical(
    NormalyzerDE:::diannPreferPrecursorInputOptionsForLimpa(
      list(
        level = "auto",
        columns = list(feature = "Protein.Group"),
        ._normalyzerde_explicit_level = FALSE
      )
    ),
    list(
      level = "auto",
      columns = list(feature = "Protein.Group"),
      ._normalyzerde_explicit_level = FALSE
    )
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
