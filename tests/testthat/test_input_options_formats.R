context("input option helpers")

test_that("input option helpers validate sep and block typos", {
  expect_true(is.list(defaultInputOptions()))
  expect_true(is.list(proteiosInputOptions()))
  expect_true(is.list(maxQuantInputOptions()))

  expect_error(defaultInputOptions(sep = ""), "sep")
  expect_error(proteiosInputOptions(sep = NA), "sep")
  expect_error(maxQuantInputOptions(sep = ""), "sep")

  expect_error(defaultInputOptions(s = ","), "Unknown argument")
  expect_error(proteiosInputOptions(s = ","), "Unknown argument")
  expect_error(maxQuantInputOptions(s = ","), "Unknown argument")
})

test_that("setupRawDataObject supports custom delimiter for default input", {
  tmpDir <- tempdir()
  dataPath <- file.path(tmpDir, "default_sep.csv")
  designPath <- file.path(tmpDir, "default_design.tsv")

  writeLines(
    c(
      "id,S1,S2",
      "f1,1,2",
      "f2,3,4"
    ),
    con = dataPath
  )

  design <- data.frame(
    sample = c("S1", "S2"),
    group = c("A", "B"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  utils::write.table(
    design,
    file = designPath,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  se <- setupRawDataObject(
    dataPath = dataPath,
    designPath = designPath,
    inputFormat = "default",
    inputOptions = defaultInputOptions(sep = ","),
    sampleColName = "sample",
    groupColName = "group"
  )

  expect_identical(colnames(SummarizedExperiment::assay(se)), c("S1", "S2"))
  expect_equal(as.numeric(SummarizedExperiment::assay(se)[1, "S1"]), 1)
  expect_equal(as.numeric(SummarizedExperiment::assay(se)[2, "S2"]), 4)
})

test_that("setupRawDataObject supports custom delimiter for Proteios input", {
  tmpDir <- tempdir()
  dataPath <- file.path(tmpDir, "proteios_sep.tsv")
  designPath <- file.path(tmpDir, "proteios_design.tsv")

  writeLines(
    c(
      "id;S1;S2",
      "f1;1;2",
      "f2;3;4"
    ),
    con = dataPath
  )

  design <- data.frame(
    sample = c("S1", "S2"),
    group = c("A", "B"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  utils::write.table(
    design,
    file = designPath,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  se <- setupRawDataObject(
    dataPath = dataPath,
    designPath = designPath,
    inputFormat = "proteios",
    inputOptions = proteiosInputOptions(sep = ";"),
    sampleColName = "sample",
    groupColName = "group"
  )

  expect_identical(colnames(SummarizedExperiment::assay(se)), c("S1", "S2"))
  expect_equal(as.numeric(SummarizedExperiment::assay(se)[1, "S1"]), 1)
  expect_equal(as.numeric(SummarizedExperiment::assay(se)[2, "S2"]), 4)
})

test_that("setupRawDataObject supports custom delimiter for MaxQuant input", {
  tmpDir <- tempdir()
  dataPath <- file.path(tmpDir, "maxquant_pep.csv")
  designPath <- file.path(tmpDir, "maxquant_design.tsv")

  maxQuantPep <- data.frame(
    Sequence = c("PEPTIDE1", "PEPTIDE2"),
    Mass = c(1000, 2000),
    Proteins = c("P1", "P2"),
    Leading.razor.protein = c("P1", "P2"),
    PEP = c(0.01, 0.02),
    Charges = c(2, 2),
    Intensity.S1 = c(10, 30),
    Intensity.S2 = c(20, 40),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  utils::write.table(
    maxQuantPep,
    file = dataPath,
    sep = ",",
    row.names = FALSE,
    quote = FALSE
  )

  design <- data.frame(
    sample = c("S1", "S2"),
    group = c("A", "B"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  utils::write.table(
    design,
    file = designPath,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  se <- setupRawDataObject(
    dataPath = dataPath,
    designPath = designPath,
    inputFormat = "maxquantpep",
    inputOptions = maxQuantInputOptions(sep = ","),
    sampleColName = "sample",
    groupColName = "group"
  )

  expect_identical(colnames(SummarizedExperiment::assay(se)), c("S1", "S2"))
  expect_equal(as.numeric(SummarizedExperiment::assay(se)[1, "S1"]), 10)
  expect_equal(as.numeric(SummarizedExperiment::assay(se)[2, "S2"]), 40)
})
