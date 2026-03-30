context("DIANN input")

test_that("setupRawDataObject reads DIANN pg_matrix and aligns sample names", {
  tmpDir <- withr::local_tempdir(pattern = "diann_pg_matrix_")

  diannMatrix <- data.frame(
    `Protein.Group` = c("P1", "P2"),
    `Protein.Names` = c("Prot1", "Prot2"),
    Genes = c("G1", "G2"),
    `First.Protein.Description` = c("Desc1", "Desc2"),
    `N.Sequences` = c(1, 2),
    `N.Proteotypic.Sequences` = c(1, 2),
    `D:\\path\\S1.raw` = c(100, 0),
    `D:\\path\\S2.raw` = c(200, 300),
    check.names = FALSE
  )

  design <- nd_two_sample_design()

  paths <- nd_write_data_and_design(
    tmp_dir = tmpDir,
    data = diannMatrix,
    design = design,
    data_name = "diann_pg_matrix.tsv",
    design_name = "diann_design.tsv"
  )

  se <- setupRawDataObject(
    dataPath = paths$dataPath,
    designPath = paths$designPath,
    inputFormat = "diann",
    sampleColName = "sample",
    groupColName = "group"
  )

  expect_true(identical(
    colnames(SummarizedExperiment::assay(se)),
    c("S1", "S2")
  ))

  rowDf <- as.data.frame(SummarizedExperiment::rowData(se), optional = TRUE)
  expect_true(all(c("Protein.Group", "Protein.Names") %in% colnames(rowDf)))
})

test_that("setupRawDataObject reads DIANN report.tsv and aggregates duplicates", {
  tmpDir <- withr::local_tempdir(pattern = "diann_report_")

  diannReport <- data.frame(
    `File.Name` = c(
      "D:\\path\\S1.raw",
      "D:\\path\\S1.raw",
      "D:\\path\\S2.raw",
      "D:\\path\\S2.raw"
    ),
    Run = c("S1", "S1", "S2", "S2"),
    `Protein.Group` = c("P1", "P1", "P1", "P2"),
    `Protein.Names` = c("Prot1", "Prot1", "Prot1", "Prot2"),
    Genes = c("G1", "G1", "G1", "G2"),
    `PG.Quantity` = c(100, 150, 200, 50),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  design <- nd_two_sample_design()

  paths <- nd_write_data_and_design(
    tmp_dir = tmpDir,
    data = diannReport,
    design = design,
    data_name = "diann_report.tsv",
    design_name = "diann_design_report.tsv"
  )

  se <- setupRawDataObject(
    dataPath = paths$dataPath,
    designPath = paths$designPath,
    inputFormat = "diann",
    sampleColName = "sample",
    groupColName = "group"
  )

  dataMatrix <- SummarizedExperiment::assay(se)
  rowDf <- as.data.frame(SummarizedExperiment::rowData(se), optional = TRUE)

  p1Index <- which(rowDf$Protein.Group == "P1")
  p2Index <- which(rowDf$Protein.Group == "P2")

  expect_true(length(p1Index) == 1)
  expect_true(length(p2Index) == 1)

  expect_equal(as.numeric(dataMatrix[p1Index, "S1"]), 150)
  expect_equal(as.numeric(dataMatrix[p1Index, "S2"]), 200)
  expect_equal(as.numeric(dataMatrix[p2Index, "S1"]), NA_real_)
  expect_equal(as.numeric(dataMatrix[p2Index, "S2"]), 50)
})

test_that("DIANN report q-value filtering removes failing rows/features", {
  tmpDir <- withr::local_tempdir(pattern = "diann_qfilter_")

  diannReport <- data.frame(
    Run = c("S1", "S2", "S1", "S2"),
    `Precursor.Id` = c("pep1", "pep1", "pep2", "pep2"),
    `Protein.Group` = c("P1", "P1", "P2", "P2"),
    `Precursor.Quantity` = c(100, 200, 300, 400),
    `Q.Value` = c(0.005, 0.005, 0.02, 0.02),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  design <- nd_two_sample_design()

  paths <- nd_write_data_and_design(
    tmp_dir = tmpDir,
    data = diannReport,
    design = design,
    data_name = "diann_report_qfilter.tsv",
    design_name = "diann_design_qfilter.tsv"
  )

  se <- setupRawContrastObject(
    dataPath = paths$dataPath,
    designPath = paths$designPath,
    sampleColName = "sample",
    inputFormat = "diann",
    inputOptions = list(
      level = "precursor",
      filters = list(
        q = list(enable = TRUE, cols = c("Q.Value"), cutoffs = 0.01)
      )
    )
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

test_that("DIANN q-value filtering keeps rows with missing q-values", {
  tmpDir <- withr::local_tempdir(pattern = "diann_qfilter_missing_q_")

  diannReport <- data.frame(
    Run = c("S1", "S2", "S1", "S2"),
    `Precursor.Id` = c("pep1", "pep1", "pep2", "pep2"),
    `Protein.Group` = c("P1", "P1", "P2", "P2"),
    `Precursor.Quantity` = c(100, 200, 300, 400),
    `Q.Value` = c(NA_real_, NA_real_, 0.005, 0.005),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  paths <- nd_write_data_and_design(
    tmp_dir = tmpDir,
    data = diannReport,
    data_name = "diann_report_qfilter_missing_q.tsv",
    design_name = "diann_design_qfilter_missing_q.tsv"
  )

  se <- setupRawContrastObject(
    dataPath = paths$dataPath,
    designPath = paths$designPath,
    sampleColName = "sample",
    inputFormat = "diann",
    inputOptions = diannInputOptions(
      level = "precursor",
      quantityCol = "Precursor.Quantity",
      qCols = "Q.Value",
      qCutoffs = 0.01,
      rt = FALSE
    )
  )

  mat <- SummarizedExperiment::assay(se)
  ann <- as.data.frame(
    SummarizedExperiment::rowData(se),
    stringsAsFactors = FALSE
  )

  expect_equal(nrow(mat), 2)
  expect_equal(ann$Precursor.Id, c("pep1", "pep2"))

  pep1 <- which(ann$Precursor.Id == "pep1")
  pep2 <- which(ann$Precursor.Id == "pep2")

  expect_equal(as.numeric(mat[pep1, ]), c(100, 200))
  expect_equal(as.numeric(mat[pep2, ]), c(300, 400))
})

test_that("DIANN report precursor-level reading adds median RT annotation", {
  tmpDir <- withr::local_tempdir(pattern = "diann_rt_")

  diannReport <- data.frame(
    Run = c("S1", "S1", "S2", "S2"),
    `Precursor.Id` = c("pep1", "pep2", "pep1", "pep2"),
    `Protein.Group` = c("P1", "P2", "P1", "P2"),
    `Precursor.Quantity` = c(100, 200, 300, 400),
    RT = c(10, 20, 30, 40),
    `Q.Value` = c(0.001, 0.001, 0.001, 0.001),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  design <- nd_two_sample_design()

  paths <- nd_write_data_and_design(
    tmp_dir = tmpDir,
    data = diannReport,
    design = design,
    data_name = "diann_report_rt.tsv",
    design_name = "diann_design_rt.tsv"
  )

  se <- setupRawDataObject(
    dataPath = paths$dataPath,
    designPath = paths$designPath,
    inputFormat = "diann",
    sampleColName = "sample",
    groupColName = "group",
    inputOptions = list(
      level = "precursor",
      columns = list(quantity = "Precursor.Quantity"),
      filters = list(
        q = list(enable = TRUE, cols = c("Q.Value"), cutoffs = 0.01)
      ),
      rt = list(col = "RT")
    )
  )

  rowDf <- as.data.frame(SummarizedExperiment::rowData(se), optional = TRUE)
  expect_true("RT" %in% colnames(rowDf))

  pep1 <- rowDf[rowDf$Precursor.Id == "pep1", , drop = FALSE]
  pep2 <- rowDf[rowDf$Precursor.Id == "pep2", , drop = FALSE]
  expect_true(nrow(pep1) == 1)
  expect_true(nrow(pep2) == 1)

  expect_equal(as.numeric(pep1$RT), 20)
  expect_equal(as.numeric(pep2$RT), 30)
})

test_that("DIANN min-positive threshold converts tiny values to NA", {
  tmpDir <- withr::local_tempdir(pattern = "diann_minpos_")

  diannReport <- data.frame(
    Run = c("S1", "S2"),
    `Precursor.Id` = c("pep1", "pep1"),
    `Protein.Group` = c("P1", "P1"),
    `Precursor.Quantity` = c(0.005, 0.02),
    `Q.Value` = c(0.001, 0.001),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  design <- nd_two_sample_design()

  paths <- nd_write_data_and_design(
    tmp_dir = tmpDir,
    data = diannReport,
    design = design,
    data_name = "diann_report_minpos.tsv",
    design_name = "diann_design_minpos.tsv"
  )

  se <- setupRawContrastObject(
    dataPath = paths$dataPath,
    designPath = paths$designPath,
    sampleColName = "sample",
    inputFormat = "diann",
    inputOptions = list(
      level = "precursor",
      columns = list(quantity = "Precursor.Quantity"),
      filters = list(
        min_positive = 0.01,
        q = list(enable = TRUE, cols = c("Q.Value"), cutoffs = 0.01)
      )
    )
  )

  mat <- SummarizedExperiment::assay(se)
  ann <- as.data.frame(
    SummarizedExperiment::rowData(se),
    stringsAsFactors = FALSE
  )

  pep1 <- which(ann$Precursor.Id == "pep1")
  expect_true(length(pep1) == 1)
  expect_true(is.na(as.numeric(mat[pep1, "S1"])))
  expect_equal(as.numeric(mat[pep1, "S2"]), 0.02)
})

test_that("DIANN level option selects precursor vs protein", {
  tmpDir <- withr::local_tempdir(pattern = "diann_level_")

  diannReport <- data.frame(
    Run = c("S1", "S1", "S2", "S2"),
    `Protein.Group` = c("P1", "P1", "P1", "P1"),
    `Protein.Names` = c("Prot1", "Prot1", "Prot1", "Prot1"),
    Genes = c("G1", "G1", "G1", "G1"),
    `Precursor.Id` = c("pep1", "pep2", "pep1", "pep2"),
    `Precursor.Normalised` = c(10, 0, 20, 30),
    `PG.Quantity` = c(100, 150, 200, 50),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  design <- nd_two_sample_design()

  paths <- nd_write_data_and_design(
    tmp_dir = tmpDir,
    data = diannReport,
    design = design,
    data_name = "diann_report_level.tsv",
    design_name = "diann_design_level.tsv"
  )

  seProtein <- setupRawContrastObject(
    dataPath = paths$dataPath,
    designPath = paths$designPath,
    sampleColName = "sample",
    inputFormat = "diann",
    inputOptions = list(level = "protein")
  )

  protMat <- SummarizedExperiment::assay(seProtein)
  protAnn <- as.data.frame(
    SummarizedExperiment::rowData(seProtein),
    stringsAsFactors = FALSE
  )
  expect_true(nrow(protMat) == 1)
  expect_true("Protein.Group" %in% colnames(protAnn))
  expect_equal(as.numeric(protMat[1, "S1"]), 150)
  expect_equal(as.numeric(protMat[1, "S2"]), 200)

  sePrec <- setupRawContrastObject(
    dataPath = paths$dataPath,
    designPath = paths$designPath,
    sampleColName = "sample",
    inputFormat = "diann",
    inputOptions = list(level = "precursor")
  )

  precMat <- SummarizedExperiment::assay(sePrec)
  precAnn <- as.data.frame(
    SummarizedExperiment::rowData(sePrec),
    stringsAsFactors = FALSE
  )
  expect_true(nrow(precMat) == 2)
  expect_true("Precursor.Id" %in% colnames(precAnn))

  pep1 <- which(precAnn$Precursor.Id == "pep1")
  pep2 <- which(precAnn$Precursor.Id == "pep2")
  expect_true(length(pep1) == 1)
  expect_true(length(pep2) == 1)

  expect_equal(as.numeric(precMat[pep1, "S1"]), 10)
  expect_equal(as.numeric(precMat[pep1, "S2"]), 20)
  expect_true(is.na(as.numeric(precMat[pep2, "S1"])))
  expect_equal(as.numeric(precMat[pep2, "S2"]), 30)
})

test_that("setupRawContrastObject reads DIANN report.parquet", {
  testthat::skip_if_not_installed("arrow")

  tmpDir <- withr::local_tempdir(pattern = "diann_parquet_")
  dataPath <- file.path(tmpDir, "diann_report.parquet")
  designPath <- file.path(tmpDir, "diann_design_parquet.tsv")

  diannParquet <- data.frame(
    Run = c("S1", "S2", "S1"),
    `Protein.Group` = c("P1", "P1", "P2"),
    `Protein.Names` = c("Prot1", "Prot1", "Prot2"),
    Genes = c("G1", "G1", "G2"),
    `PG.MaxLFQ` = c(1000, 2000, 3000),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  arrow::write_parquet(diannParquet, dataPath)

  design <- nd_two_sample_design()
  nd_write_table(design, designPath)

  se <- setupRawContrastObject(
    dataPath = dataPath,
    designPath = designPath,
    sampleColName = "sample",
    inputFormat = "diann"
  )

  dataMatrix <- SummarizedExperiment::assay(se)
  expect_true(is.matrix(dataMatrix))
  expect_true(is.numeric(dataMatrix))
  expect_true(identical(colnames(dataMatrix), c("S1", "S2")))
})

test_that("diannChooseReportSpec validates and infers sample/feature/quantity columns", {
  expect_error(
    NormalyzerDE:::diannChooseReportSpec(c("Protein.Group", "PG.Quantity")),
    class = "normalyzerde_error"
  )

  expect_error(
    NormalyzerDE:::diannChooseReportSpec(
      c("Run", "Protein.Group", "PG.Quantity"),
      diannSampleCol = "Missing.Sample"
    ),
    class = "normalyzerde_error"
  )

  expect_error(
    NormalyzerDE:::diannChooseReportSpec(
      c("Run", "Protein.Group", "PG.Quantity"),
      diannFeatureCol = "Missing.Feature"
    ),
    class = "normalyzerde_error"
  )

  expect_error(
    NormalyzerDE:::diannChooseReportSpec(
      c("Run", "Protein.Group", "PG.Quantity"),
      diannQuantityCol = "Missing.Quantity"
    ),
    class = "normalyzerde_error"
  )

  expect_error(
    NormalyzerDE:::diannChooseReportSpec(
      c("Run", "CustomFeature"),
      diannFeatureCol = "CustomFeature"
    ),
    class = "normalyzerde_error"
  )

  expect_error(
    NormalyzerDE:::diannChooseReportSpec(
      c("Run", "Precursor.Id", "Precursor.Quantity"),
      diannLevel = "protein",
      diannQuantityCol = "Precursor.Quantity"
    ),
    class = "normalyzerde_error"
  )

  inferredProtein <- NormalyzerDE:::diannChooseReportSpec(
    c(
      "Run",
      "Protein.Group",
      "Protein.Names",
      "Genes",
      "PG.Quantity",
      "Decoy"
    ),
    diannLevel = "auto"
  )
  expect_equal(inferredProtein$sampleCol, "Run")
  expect_equal(inferredProtein$featureCol, "Protein.Group")
  expect_equal(inferredProtein$quantityCol, "PG.Quantity")
  expect_true(all(
    c("Protein.Group", "Protein.Names") %in% inferredProtein$extraCols
  ))
  expect_equal(inferredProtein$decoyCol, "Decoy")

  inferredPrec <- NormalyzerDE:::diannChooseReportSpec(
    c(
      "File.Name",
      "Precursor.Id",
      "Modified.Sequence",
      "Precursor.Quantity",
      "Q.Value"
    ),
    diannLevel = "auto"
  )
  expect_equal(inferredPrec$sampleCol, "File.Name")
  expect_equal(inferredPrec$featureCol, "Precursor.Id")
  expect_equal(inferredPrec$quantityCol, "Precursor.Quantity")
  expect_true("Precursor.Id" %in% inferredPrec$extraCols)

  inferredFromQuantity <- NormalyzerDE:::diannChooseReportSpec(
    c(
      "Run",
      "Protein.Group",
      "Precursor.Id",
      "Precursor.Quantity"
    ),
    diannLevel = "auto",
    diannQuantityCol = "Precursor.Quantity"
  )
  expect_equal(inferredFromQuantity$featureCol, "Precursor.Id")
  expect_equal(inferredFromQuantity$quantityCol, "Precursor.Quantity")

  expect_error(
    NormalyzerDE:::diannChooseReportSpec(
      c(
        "Run",
        "Protein.Group",
        "Precursor.Id",
        "Precursor.Quantity"
      ),
      diannFeatureCol = "Protein.Group",
      diannQuantityCol = "Precursor.Quantity"
    ),
    class = "normalyzerde_error"
  )

  expect_error(
    NormalyzerDE:::diannChooseReportSpec(
      c(
        "Run",
        "Protein.Group",
        "Precursor.Id",
        "Precursor.Quantity"
      ),
      diannLevel = "protein",
      diannQuantityCol = "Precursor.Quantity"
    ),
    class = "normalyzerde_error"
  )

  expect_error(
    NormalyzerDE:::diannChooseReportSpec(
      c(
        "Run",
        "Protein.Group",
        "Precursor.Id",
        "CustomQuantity"
      ),
      diannLevel = "auto",
      diannQuantityCol = "CustomQuantity"
    ),
    class = "normalyzerde_error"
  )
})

test_that("readDiannToDataFrame warns when requested extra columns are missing", {
  tmpDir <- withr::local_tempdir(pattern = "diann_missing_extra_")

  reportPath <- file.path(tmpDir, "diann_report.tsv")
  diannReport <- data.frame(
    Run = c("S1", "S2"),
    `Protein.Group` = c("P1", "P1"),
    `Protein.Names` = c("Prot1", "Prot1"),
    `PG.Quantity` = c(100, 200),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  nd_write_table(diannReport, reportPath)

  wide <- expect_warning(
    readDiannToDataFrame(
      reportPath,
      inputOptions = list(
        level = "protein",
        columns = list(extra = c("Protein.Names", "MissingExtra"))
      )
    ),
    class = "normalyzerde_warning"
  )

  expect_true("Protein.Group" %in% colnames(wide))
  expect_true("Protein.Names" %in% colnames(wide))
  expect_false("MissingExtra" %in% colnames(wide))
})

test_that("readDiannToDataFrame errors for ambiguous custom quantity columns", {
  tmpDir <- withr::local_tempdir(pattern = "diann_custom_quantity_")

  reportPath <- file.path(tmpDir, "diann_report.tsv")
  diannReport <- data.frame(
    Run = c("S1", "S2"),
    `Protein.Group` = c("P1", "P1"),
    `Precursor.Id` = c("prec1", "prec1"),
    CustomQuantity = c(100, 200),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  nd_write_table(diannReport, reportPath)

  expect_error(
    readDiannToDataFrame(
      reportPath,
      designSampleNames = c("S1", "S2"),
      inputOptions = diannInputOptions(
        level = "auto",
        quantityCol = "CustomQuantity",
        qEnable = FALSE,
        rt = FALSE
      )
    ),
    class = "normalyzerde_error"
  )
})

test_that("readDiannToDataFrame supports empty extraCols and disabling q filtering", {
  tmpDir <- withr::local_tempdir(pattern = "diann_empty_extra_")

  reportPath <- file.path(tmpDir, "diann_report.tsv")
  diannReport <- data.frame(
    Run = c("S1", "S2"),
    `Protein.Group` = c("P1", "P1"),
    `Protein.Names` = c("Prot1", "Prot1"),
    `PG.Quantity` = c(100, 200),
    `PG.Q.Value` = c(0.001, 0.001),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  nd_write_table(diannReport, reportPath)

  wide <- readDiannToDataFrame(
    reportPath,
    inputOptions = list(
      level = "protein",
      columns = list(extra = character()),
      filters = list(q = FALSE)
    )
  )

  expect_true("Protein.Group" %in% colnames(wide))
  expect_false("Protein.Names" %in% colnames(wide))
  expect_equal(wide$S1[1], 100)
  expect_equal(wide$S2[1], 200)
})

test_that("readDiannToDataFrame warns about ambiguous auto inference only when enabled", {
  tmpDir <- withr::local_tempdir(pattern = "diann_auto_ambig_")

  reportPath <- file.path(tmpDir, "diann_report.tsv")
  diannReport <- data.frame(
    Run = c("S1", "S2"),
    `Protein.Group` = c("P1", "P1"),
    `PG.MaxLFQ` = c(1000, 2000),
    `Precursor.Id` = c("pep1", "pep1"),
    `Precursor.Quantity` = c(10000, 20000),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  nd_write_table(diannReport, reportPath)

  expect_no_warning(readDiannToDataFrame(reportPath, inputOptions = NULL))

  withr::local_options(list(NormalyzerDE.warnDiannAutoAmbiguous = TRUE))
  expect_warning(
    readDiannToDataFrame(reportPath, inputOptions = NULL),
    "both protein-level and precursor-level",
    class = "normalyzerde_warning"
  )
})

test_that("normalyzer preQuant=limpa warns for an explicit DIANN auto request", {
  tmpDir <- withr::local_tempdir(pattern = "normalyzer_limpa_diann_auto_ambig_")

  reportPath <- file.path(tmpDir, "diann_report.tsv")
  designPath <- file.path(tmpDir, "design.tsv")

  diannReport <- data.frame(
    Run = c("S1", "S2"),
    `Protein.Group` = c("P1", "P1"),
    `PG.MaxLFQ` = c(1000, 2000),
    `Precursor.Id` = c("pep1", "pep1"),
    `Precursor.Quantity` = c(10000, 20000),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  nd_write_table(diannReport, reportPath)

  design <- data.frame(
    sample = c("S1", "S2"),
    group = c("A", "B"),
    stringsAsFactors = FALSE
  )
  nd_write_table(design, designPath)

  warnings <- character()
  out <- tryCatch(
    withCallingHandlers(
      normalyzer(
        jobName = "diann_auto_ambig",
        designPath = designPath,
        dataPath = reportPath,
        outputDir = tmpDir,
        inputFormat = "diann",
        inputOptions = diannInputOptions(level = "auto", rt = FALSE),
        preQuant = "limpa",
        skipAnalysis = TRUE,
        normalizeRetentionTime = FALSE,
        requireReplicates = FALSE,
        sampleAbundThres = 1,
        quiet = TRUE
      ),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = identity
  )

  expect_true(inherits(out, "error"))
  expect_true(any(grepl("both protein-level and precursor-level", warnings)))
})

test_that("normalyzer preQuant=limpa prefers precursor-level DIANN input by default", {
  testthat::skip_if_not_installed("limpa")

  tmpDir <- withr::local_tempdir(
    pattern = "normalyzer_limpa_diann_precursor_default_"
  )

  reportPath <- file.path(tmpDir, "diann_report.tsv")
  designPath <- file.path(tmpDir, "design.tsv")
  design <- nd_balanced_diann_design()

  nd_write_table(nd_make_ambiguous_diann_report(design = design), reportPath)
  nd_write_table(design, designPath)

  outDefault <- file.path(tmpDir, "default")
  outExplicit <- file.path(tmpDir, "explicit")
  defaultWarnings <- character()
  explicitWarnings <- character()

  withCallingHandlers(
    normalyzer(
      jobName = "diann_limpa_default",
      designPath = designPath,
      dataPath = reportPath,
      outputDir = outDefault,
      inputFormat = "diann",
      preQuant = "limpa",
      skipAnalysis = TRUE,
      normalizeRetentionTime = FALSE,
      requireReplicates = FALSE,
      sampleAbundThres = 1,
      quiet = TRUE
    ),
    warning = function(w) {
      defaultWarnings <<- c(defaultWarnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  withCallingHandlers(
    normalyzer(
      jobName = "diann_limpa_explicit",
      designPath = designPath,
      dataPath = reportPath,
      outputDir = outExplicit,
      inputFormat = "diann",
      inputOptions = diannInputOptions(level = "precursor", rt = FALSE),
      preQuant = "limpa",
      skipAnalysis = TRUE,
      normalizeRetentionTime = FALSE,
      requireReplicates = FALSE,
      sampleAbundThres = 1,
      quiet = TRUE
    ),
    warning = function(w) {
      explicitWarnings <<- c(explicitWarnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  defaultRds <- readRDS(file.path(
    outDefault,
    "diann_limpa_default",
    "diann_limpa_default_limpa_quantified.rds"
  ))
  explicitRds <- readRDS(file.path(
    outExplicit,
    "diann_limpa_explicit",
    "diann_limpa_explicit_limpa_quantified.rds"
  ))

  expect_equal(defaultRds$E, explicitRds$E)
  expect_equal(defaultRds$genes, explicitRds$genes)
  expect_false(any(grepl(
    "both protein-level and precursor-level",
    defaultWarnings
  )))
  expect_false(any(grepl(
    "both protein-level and precursor-level",
    explicitWarnings
  )))
})

test_that("normalyzerDE type=limpa prefers precursor-level DIANN input by default", {
  testthat::skip_if_not_installed("limpa")

  tmpDir <- withr::local_tempdir(
    pattern = "normalyzerde_limpa_diann_precursor_default_"
  )

  reportPath <- file.path(tmpDir, "diann_report.tsv")
  designPath <- file.path(tmpDir, "design.tsv")
  design <- nd_balanced_diann_design()

  nd_write_table(nd_make_ambiguous_diann_report(design = design), reportPath)
  nd_write_table(design, designPath)

  outDefault <- file.path(tmpDir, "default")
  outExplicit <- file.path(tmpDir, "explicit")

  suppressWarnings(normalyzerDE(
    jobName = "diann_limpa_default",
    comparisons = "A-B",
    designPath = designPath,
    dataPath = reportPath,
    outputDir = outDefault,
    inputFormat = "diann",
    type = "limpa",
    limpaProteinIdCol = "Protein.Group",
    logTrans = TRUE,
    quiet = TRUE
  ))

  suppressWarnings(normalyzerDE(
    jobName = "diann_limpa_explicit",
    comparisons = "A-B",
    designPath = designPath,
    dataPath = reportPath,
    outputDir = outExplicit,
    inputFormat = "diann",
    inputOptions = diannInputOptions(level = "precursor", rt = FALSE),
    type = "limpa",
    limpaProteinIdCol = "Protein.Group",
    logTrans = TRUE,
    quiet = TRUE
  ))

  defaultStats <- utils::read.delim(
    file.path(
      outDefault,
      "diann_limpa_default",
      "diann_limpa_default_stats.tsv"
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  explicitStats <- utils::read.delim(
    file.path(
      outExplicit,
      "diann_limpa_explicit",
      "diann_limpa_explicit_stats.tsv"
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  expect_equal(defaultStats, explicitStats)
})

test_that("normalyzerDE type=limpa preserves an explicit DIANN auto request", {
  testthat::skip_if_not_installed("limpa")

  tmpDir <- withr::local_tempdir(
    pattern = "normalyzerde_limpa_diann_explicit_auto_"
  )

  reportPath <- file.path(tmpDir, "diann_report.tsv")
  designPath <- file.path(tmpDir, "design.tsv")

  diannReport <- data.frame(
    Run = c("S1", "S2"),
    `Protein.Group` = c("P1", "P1"),
    `PG.MaxLFQ` = c(1000, 2000),
    `Precursor.Id` = c("pep1", "pep1"),
    `Precursor.Quantity` = c(10000, 20000),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  nd_write_table(diannReport, reportPath)

  design <- data.frame(
    sample = c("S1", "S2"),
    group = c("A", "B"),
    stringsAsFactors = FALSE
  )
  nd_write_table(design, designPath)

  warnings <- character()
  err <- tryCatch(
    withCallingHandlers(
      normalyzerDE(
        jobName = "diann_limpa_auto",
        comparisons = "A-B",
        designPath = designPath,
        dataPath = reportPath,
        outputDir = tmpDir,
        inputFormat = "diann",
        inputOptions = diannInputOptions(level = "auto", rt = FALSE),
        type = "limpa",
        limpaProteinIdCol = "Protein.Group",
        logTrans = TRUE,
        quiet = TRUE
      ),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = identity
  )

  expect_true(inherits(err, "error"))
  expect_true(any(grepl("both protein-level and precursor-level", warnings)))
})

test_that("diannReadReportParquet errors when requested columns are missing", {
  testthat::skip_if_not_installed("arrow")

  tmpDir <- withr::local_tempdir(pattern = "diann_parquet_missing_cols_")

  dataPath <- file.path(tmpDir, "diann_report.parquet")
  arrow::write_parquet(
    data.frame(
      Run = c("S1"),
      `Protein.Group` = c("P1"),
      `PG.Quantity` = c(100),
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    dataPath
  )

  expect_error(
    NormalyzerDE:::diannReadReportParquet(
      dataPath,
      selectCols = c("Run", "Protein.Group", "PG.Quantity", "MissingCol")
    ),
    class = "normalyzerde_error"
  )
})

test_that("readDiannToDataFrame reads parquet DIANN reports and applies filters", {
  testthat::skip_if_not_installed("arrow")

  tmpDir <- withr::local_tempdir(pattern = "diann_parquet_read_")
  reportPath <- file.path(tmpDir, "diann_report.parquet")

  diannReport <- data.frame(
    Run = c("S1", "S2", "S1", "S2"),
    `Protein.Group` = c("P1", "P1", "P2", "P2"),
    `PG.Quantity` = c(100, 200, 300, 400),
    Decoy = c(0, 0, 1, 0),
    `PG.Q.Value` = c(0.005, 0.005, 0.005, 0.02),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  arrow::write_parquet(diannReport, reportPath)

  wideNoQ <- readDiannToDataFrame(
    reportPath,
    inputOptions = list(
      level = "protein",
      columns = list(extra = character()),
      filters = list(q = FALSE)
    )
  )
  expect_true(all(c("S1", "S2") %in% colnames(wideNoQ)))

  wideFiltered <- expect_warning(
    readDiannToDataFrame(
      reportPath,
      inputOptions = list(
        level = "protein",
        columns = list(extra = "MissingExtra"),
        filters = list(q = list(enable = TRUE, cutoffs = 0.01))
      )
    ),
    class = "normalyzerde_warning"
  )

  expect_equal(as.character(wideFiltered$Protein.Group), "P1")
  expect_equal(wideFiltered$S1[1], 100)
  expect_equal(wideFiltered$S2[1], 200)
})

test_that("readDiannToDataFrame errors on duplicate sample names after stripping paths/extensions", {
  tmpDir <- withr::local_tempdir(pattern = "diann_dup_samples_")

  matrixPath <- file.path(tmpDir, "diann_pg_matrix.tsv")
  diannMatrix <- data.frame(
    `Protein.Group` = c("P1", "P2"),
    `S1.raw` = c(1, 2),
    `S1.mzML` = c(3, 4),
    check.names = FALSE
  )
  nd_write_table(diannMatrix, matrixPath)

  expect_error(
    readDiannToDataFrame(matrixPath, designSampleNames = c("S1")),
    class = "normalyzerde_error"
  )
})

test_that("readDiannToDataFrame maps repeated DIA-NN report rows to cleaned design sample names", {
  tmpDir <- withr::local_tempdir(pattern = "diann_report_clean_names_")

  reportPath <- file.path(tmpDir, "diann_report.tsv")
  diannReport <- data.frame(
    `File.Name` = c(
      "/path/S1.raw",
      "/path/S1.raw",
      "/path/S2.raw",
      "/path/S2.raw"
    ),
    `Precursor.Id` = c("pep1", "pep2", "pep1", "pep2"),
    `Precursor.Quantity` = c(100, 150, 200, 250),
    `Q.Value` = c(0.001, 0.001, 0.001, 0.001),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  nd_write_table(diannReport, reportPath)

  out <- readDiannToDataFrame(
    reportPath,
    designSampleNames = c("S1", "S2"),
    inputOptions = diannInputOptions(
      level = "precursor",
      sampleCol = "File.Name",
      quantityCol = "Precursor.Quantity",
      qCols = "Q.Value",
      rt = FALSE
    )
  )

  expect_true(all(c("S1", "S2") %in% colnames(out)))
  expect_equal(out$S1, c(100, 150))
  expect_equal(out$S2, c(200, 250))
})

test_that("readDiannToDataFrame applies min_positive filtering for DIA-NN matrices", {
  tmpDir <- withr::local_tempdir(pattern = "diann_matrix_minpos_")

  matrixPath <- file.path(tmpDir, "diann_pg_matrix.tsv")
  diannMatrix <- data.frame(
    `Protein.Group` = c("P1", "P2"),
    `S1.raw` = c(0.005, 0),
    `S2.raw` = c(0.02, 0.009),
    check.names = FALSE
  )
  nd_write_table(diannMatrix, matrixPath)

  out <- readDiannToDataFrame(
    matrixPath,
    designSampleNames = c("S1", "S2"),
    inputOptions = list(filters = list(min_positive = 0.01))
  )

  expect_equal(colnames(out)[2:3], c("S1", "S2"))
  expect_true(is.na(out$S1[1]))
  expect_true(is.na(out$S1[2]))
  expect_equal(out$S2[1], 0.02)
  expect_true(is.na(out$S2[2]))
})

test_that("readDiannToDataFrame maps File.Name paths to design sample names", {
  tmpDir <- withr::local_tempdir(pattern = "diann_sample_map_")

  reportPath <- file.path(tmpDir, "diann_report_paths.tsv")
  diannReport <- data.frame(
    `File.Name` = c("some/dir/S1.raw", "some/dir/S2.raw"),
    `Protein.Group` = c("P1", "P1"),
    `PG.Quantity` = c(100, 200),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  nd_write_table(diannReport, reportPath)

  wide <- readDiannToDataFrame(
    reportPath,
    designSampleNames = c("S1", "S2"),
    inputOptions = list(level = "protein")
  )

  expect_true(all(c("S1", "S2") %in% colnames(wide)))
  expect_equal(wide$S1[1], 100)
  expect_equal(wide$S2[1], 200)
})

test_that("diannFilterDecoys supports logical, numeric and string decoy columns", {
  reportLogical <- data.frame(
    x = 1:3,
    Decoy = c(TRUE, FALSE, NA),
    stringsAsFactors = FALSE
  )
  outLogical <- NormalyzerDE:::diannFilterDecoys(
    reportLogical,
    decoyCol = "Decoy"
  )
  expect_equal(nrow(outLogical), 2)

  reportNumeric <- data.frame(
    x = 1:3,
    Decoy = c(1, 0, NA),
    stringsAsFactors = FALSE
  )
  outNumeric <- NormalyzerDE:::diannFilterDecoys(
    reportNumeric,
    decoyCol = "Decoy"
  )
  expect_equal(nrow(outNumeric), 2)

  reportString <- data.frame(
    x = 1:7,
    Decoy = c("TRUE", "false", "1", "0", "yes", "no", NA),
    stringsAsFactors = FALSE
  )
  outString <- NormalyzerDE:::diannFilterDecoys(
    reportString,
    decoyCol = "Decoy"
  )
  expect_equal(nrow(outString), 4)
  expect_true(all(outString$x %in% c(2, 4, 6, 7)))
})
