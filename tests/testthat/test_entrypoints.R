context("entrypoints")

test_that("normalyzer writes outputs when given SummarizedExperiment", {
  data(example_summarized_experiment)

  outDir <- withr::local_tempdir(pattern = "entry_norm_se_")
  jobName <- "entry_norm_se"
  expectedDir <- file.path(outDir, NormalyzerDE:::sanitizeJobName(jobName))

  out <- suppressWarnings(normalyzer(
    jobName = jobName,
    experimentObj = example_summarized_experiment,
    outputDir = outDir,
    normalizeRetentionTime = FALSE,
    skipAnalysis = TRUE,
    quiet = TRUE
  ))

  expect_null(out)
  expect_true(dir.exists(expectedDir))
  expect_true(file.exists(file.path(expectedDir, "submitted_rawdata.txt")))
  expect_gt(
    length(list.files(
      expectedDir,
      pattern = "-normalized\\.txt$",
      full.names = TRUE
    )),
    0
  )
})

test_that("normalyzerDE can compute one-vs-rest without explicit comparisons", {
  mat <- matrix(stats::rnorm(20 * 4, mean = 10, sd = 1), nrow = 20)
  colnames(mat) <- paste0("s", seq_len(ncol(mat)))

  design <- data.frame(
    sample = colnames(mat),
    group = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )
  rownames(design) <- design$sample

  se <- SummarizedExperiment::SummarizedExperiment(
    assay = mat,
    colData = design,
    rowData = data.frame(feature = paste0("f", seq_len(nrow(mat))))
  )

  outDir <- withr::local_tempdir(pattern = "onevsrest_")
  jobName <- "onevsrest"
  expectedDir <- file.path(outDir, NormalyzerDE:::sanitizeJobName(jobName))

  out <- suppressWarnings(normalyzerDE(
    jobName = jobName,
    comparisons = NULL,
    experimentObj = se,
    outputDir = outDir,
    type = "limma",
    oneVsRest = TRUE,
    quiet = TRUE
  ))

  expect_null(out)

  outStatsPath <- file.path(
    expectedDir,
    paste0(basename(expectedDir), "_stats.tsv")
  )
  expect_true(file.exists(outStatsPath))

  restLabel <- NormalyzerDE:::chooseOneVsRestLabel(design$group)
  expectedPCols <- paste0(c("A", "B"), "-", restLabel, "_PValue")

  outDf <- utils::read.table(
    outStatsPath,
    sep = "\t",
    header = TRUE,
    check.names = FALSE
  )
  expect_true(all(expectedPCols %in% colnames(outDf)))
})

test_that("normalyzer supports DIANN report precursors with RT normalization", {
  tmpDir <- withr::local_tempdir(pattern = "diann_normalyzer_")
  dataPath <- file.path(tmpDir, "diann_report_for_normalyzer.tsv")
  designPath <- file.path(tmpDir, "diann_design_for_normalyzer.tsv")

  nFeatures <- 120
  precursors <- paste0("pep", seq_len(nFeatures))
  proteins <- paste0("P", seq_len(nFeatures))

  diannReport <- data.frame(
    Run = rep(c("S1", "S2"), each = nFeatures),
    `Precursor.Id` = rep(precursors, times = 2),
    `Protein.Group` = rep(proteins, times = 2),
    `Precursor.Quantity` = c(
      seq(1000, length.out = nFeatures),
      seq(2000, length.out = nFeatures)
    ),
    RT = c(seq(1, length.out = nFeatures), seq(1.5, length.out = nFeatures)),
    `Q.Value` = rep(0.001, 2 * nFeatures),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  nd_write_table(diannReport, dataPath)
  nd_write_table(nd_two_sample_design(), designPath)

  jobName <- "entry_diann_norm"
  expectedDir <- file.path(tmpDir, NormalyzerDE:::sanitizeJobName(jobName))

  out <- suppressWarnings(normalyzer(
    jobName = jobName,
    designPath = designPath,
    dataPath = dataPath,
    outputDir = tmpDir,
    sampleAbundThres = 5,
    requireReplicates = FALSE,
    inputFormat = "diann",
    inputOptions = list(
      level = "precursor",
      columns = list(quantity = "Precursor.Quantity"),
      filters = list(
        q = list(enable = TRUE, cols = c("Q.Value"), cutoffs = 0.01)
      ),
      rt = list(col = "RT")
    ),
    normalizeRetentionTime = TRUE,
    rtStepSizeMinutes = 1000,
    rtWindowMinCount = 2,
    skipAnalysis = TRUE,
    quiet = TRUE
  ))

  expect_null(out)
  expect_true(file.exists(file.path(expectedDir, "RT-median-normalized.txt")))
  expect_true(file.exists(file.path(expectedDir, "RT-mean-normalized.txt")))
  expect_true(file.exists(file.path(expectedDir, "RT-Loess-normalized.txt")))
})

test_that("normalyzer emits version message and errors when inputs missing", {
  expect_error(
    normalyzer(
      jobName = "missing_inputs",
      quiet = FALSE
    ),
    class = "normalyzerde_error"
  )
})

test_that("normalyzer stays quiet when inputs missing and quiet=TRUE", {
  expect_error(
    normalyzer(
      jobName = "missing_inputs",
      quiet = TRUE
    ),
    class = "normalyzerde_error"
  )
})

test_that("normalyzer runs end-to-end with analysis/plots on tiny example", {
  dataPath <- system.file(package = "NormalyzerDE", "extdata", "tiny_data.tsv")
  designPath <- system.file(
    package = "NormalyzerDE",
    "extdata",
    "tiny_design.tsv"
  )
  expect_true(nzchar(dataPath))
  expect_true(nzchar(designPath))

  tmpDir <- withr::local_tempdir(pattern = "normalyzer_e2e_")
  jobName <- "entry_e2e"
  expectedDir <- file.path(tmpDir, NormalyzerDE:::sanitizeJobName(jobName))

  out <- suppressWarnings(normalyzer(
    jobName = jobName,
    designPath = designPath,
    dataPath = dataPath,
    outputDir = tmpDir,
    normalizeRetentionTime = FALSE,
    skipAnalysis = FALSE,
    quiet = TRUE
  ))

  expect_null(out)
  expect_true(file.exists(file.path(expectedDir, "log2-normalized.txt")))
  expect_true(file.exists(file.path(
    expectedDir,
    paste0("Norm-report-", basename(expectedDir), ".pdf")
  )))
})

test_that("normalyzerDE errors when comparisons are missing and oneVsRest=FALSE", {
  mat <- matrix(stats::rnorm(20 * 4, mean = 10, sd = 1), nrow = 20)
  colnames(mat) <- paste0("s", seq_len(ncol(mat)))

  design <- data.frame(
    sample = colnames(mat),
    group = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )
  rownames(design) <- design$sample

  se <- SummarizedExperiment::SummarizedExperiment(
    assay = mat,
    colData = design,
    rowData = data.frame(feature = paste0("f", seq_len(nrow(mat))))
  )

  expect_error(
    normalyzerDE(
      jobName = "missing_comps",
      comparisons = NULL,
      experimentObj = se,
      quiet = TRUE
    ),
    class = "normalyzerde_error"
  )
})

test_that("normalyzerDE emits version message and errors when inputs missing", {
  expect_error(
    normalyzerDE(jobName = "missing_inputs", quiet = FALSE),
    class = "normalyzerde_error"
  )
})

test_that("normalyzerDE stays quiet when inputs missing and quiet=TRUE", {
  expect_error(
    normalyzerDE(jobName = "missing_inputs", quiet = TRUE),
    class = "normalyzerde_error"
  )
})
