context("entrypoints")

test_that("normalyzer writes outputs when given SummarizedExperiment", {
  data(example_summarized_experiment)

  jobName <- paste0("entry_norm_", sample.int(1e9, 1))
  outDir <- tempdir()
  expectedDir <- file.path(outDir, NormalyzerDE:::sanitizeJobName(jobName))
  if (dir.exists(expectedDir)) {
    unlink(expectedDir, recursive = TRUE, force = TRUE)
  }
  on.exit(unlink(expectedDir, recursive = TRUE, force = TRUE), add = TRUE)

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

  jobName <- paste0("onevsrest_", sample.int(1e9, 1))
  outDir <- tempdir()
  expectedDir <- file.path(outDir, NormalyzerDE:::sanitizeJobName(jobName))
  if (dir.exists(expectedDir)) {
    unlink(expectedDir, recursive = TRUE, force = TRUE)
  }
  on.exit(unlink(expectedDir, recursive = TRUE, force = TRUE), add = TRUE)

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
  tmpDir <- tempdir()
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

  jobName <- paste0("entry_diann_norm_", sample.int(1e9, 1))
  outDir <- tempdir()
  expectedDir <- file.path(outDir, NormalyzerDE:::sanitizeJobName(jobName))
  if (dir.exists(expectedDir)) {
    unlink(expectedDir, recursive = TRUE, force = TRUE)
  }
  on.exit(unlink(expectedDir, recursive = TRUE, force = TRUE), add = TRUE)

  out <- suppressWarnings(normalyzer(
    jobName = jobName,
    designPath = designPath,
    dataPath = dataPath,
    outputDir = outDir,
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
