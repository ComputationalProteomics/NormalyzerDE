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

  design <- nd_make_design(c("A", "A", "B", "B"))
  se <- nd_make_summarized_experiment(
    assay = mat,
    groups = design$group
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

  paths <- nd_write_data_and_design(
    tmp_dir = tmpDir,
    data = diannReport,
    data_name = "diann_report_for_normalyzer.tsv",
    design_name = "diann_design_for_normalyzer.tsv"
  )

  jobName <- "entry_diann_norm"
  expectedDir <- file.path(tmpDir, NormalyzerDE:::sanitizeJobName(jobName))

  out <- suppressWarnings(normalyzer(
    jobName = jobName,
    designPath = paths$designPath,
    dataPath = paths$dataPath,
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
  captured <- nd_capture_conditions(
    normalyzer(
      jobName = "missing_inputs",
      quiet = FALSE
    )
  )

  expect_s3_class(captured$error, "normalyzerde_error")
  expect_gt(length(captured$messages), 0)
  expect_true(any(grepl("NormalyzerDE", captured$messages, fixed = TRUE)))
  expect_true(any(grepl("version", captured$messages, ignore.case = TRUE)))
})

test_that("normalyzer stays quiet when inputs missing and quiet=TRUE", {
  captured <- nd_capture_conditions(
    normalyzer(
      jobName = "missing_inputs",
      quiet = TRUE
    )
  )

  expect_s3_class(captured$error, "normalyzerde_error")
  expect_length(captured$messages, 0)
})

test_that("normalyzer errors on non-empty output directories unless reuseOutputDir=TRUE", {
  data(example_summarized_experiment)

  outDir <- withr::local_tempdir(pattern = "entry_norm_reuse_")
  jobName <- "entry_norm_reuse"
  expectedDir <- file.path(outDir, NormalyzerDE:::sanitizeJobName(jobName))
  dir.create(expectedDir)
  writeLines("stale", file.path(expectedDir, "stale.txt"))

  expect_error(
    normalyzer(
      jobName = jobName,
      experimentObj = example_summarized_experiment,
      outputDir = outDir,
      normalizeRetentionTime = FALSE,
      skipAnalysis = TRUE,
      quiet = TRUE
    ),
    class = "normalyzerde_error"
  )

  out <- suppressWarnings(normalyzer(
    jobName = jobName,
    experimentObj = example_summarized_experiment,
    outputDir = outDir,
    normalizeRetentionTime = FALSE,
    skipAnalysis = TRUE,
    quiet = TRUE,
    reuseOutputDir = TRUE
  ))

  expect_null(out)
  expect_true(file.exists(file.path(expectedDir, "submitted_rawdata.txt")))
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

  se <- nd_make_summarized_experiment(
    assay = mat,
    groups = c("A", "A", "B", "B")
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

test_that("normalyzerDE validates limpaOptions helper input", {
  mat <- matrix(stats::rnorm(20 * 4, mean = 10, sd = 1), nrow = 20)
  colnames(mat) <- paste0("s", seq_len(ncol(mat)))

  se <- nd_make_summarized_experiment(
    assay = mat,
    groups = c("A", "A", "B", "B")
  )

  expect_error(
    normalyzerDE(
      jobName = "bad_limpa_options",
      comparisons = "A-B",
      experimentObj = se,
      type = "limpa",
      quiet = TRUE,
      limpaOptions = list(unknown = TRUE)
    ),
    class = "normalyzerde_error"
  )
})

test_that("normalyzerDE errors on non-empty output directories unless reuseOutputDir=TRUE", {
  mat <- matrix(stats::rnorm(20 * 4, mean = 10, sd = 1), nrow = 20)
  colnames(mat) <- paste0("s", seq_len(ncol(mat)))

  se <- nd_make_summarized_experiment(
    assay = mat,
    groups = c("A", "A", "B", "B")
  )

  outDir <- withr::local_tempdir(pattern = "entry_de_reuse_")
  jobName <- "entry_de_reuse"
  expectedDir <- file.path(outDir, NormalyzerDE:::sanitizeJobName(jobName))
  dir.create(expectedDir)
  writeLines("stale", file.path(expectedDir, "stale.txt"))

  expect_error(
    normalyzerDE(
      jobName = jobName,
      comparisons = "A-B",
      experimentObj = se,
      outputDir = outDir,
      type = "limma",
      quiet = TRUE
    ),
    class = "normalyzerde_error"
  )

  out <- suppressWarnings(normalyzerDE(
    jobName = jobName,
    comparisons = "A-B",
    experimentObj = se,
    outputDir = outDir,
    type = "limma",
    quiet = TRUE,
    reuseOutputDir = TRUE
  ))

  expect_null(out)
  expect_true(file.exists(file.path(
    expectedDir,
    paste0(basename(expectedDir), "_stats.tsv")
  )))
})

test_that("normalyzerDE writes limpa sample weights when requested", {
  testthat::skip_if_not_installed("limpa")

  test_data <- nd_limpa_de_matrix()

  tmpDir <- withr::local_tempdir(pattern = "normalyzerde_limpa_weights_")
  paths <- nd_write_limpa_input_fixture(
    tmp_dir = tmpDir,
    mat = test_data,
    data_name = "limpa_weights_data.tsv",
    design_name = "limpa_weights_design.tsv"
  )

  jobName <- "limpa_weights"
  expectedDir <- file.path(tmpDir, NormalyzerDE:::sanitizeJobName(jobName))

  out <- nd_run_limpa_de_from_paths(
    paths = paths,
    out_dir = tmpDir,
    job_name = jobName,
    limpa_options = limpaOptions(
      deArgs = list(sample.weights = TRUE),
      quantArgs = list(chunk = 10L)
    )
  )

  expect_null(out)

  weightsPath <- file.path(expectedDir, paste0(jobName, "_sample_weights.tsv"))
  expect_true(file.exists(weightsPath))

  weights <- utils::read.delim(
    weightsPath,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  expect_equal(colnames(weights), c("comparison", "sample", "sampleWeight"))
  expect_equal(weights$comparison, rep(".global", ncol(test_data)))
  expect_equal(weights$sample, colnames(test_data))
  expect_false(anyNA(weights$sampleWeight))
})

test_that("normalyzerDE emits version message and errors when inputs missing", {
  captured <- nd_capture_conditions(
    normalyzerDE(jobName = "missing_inputs", quiet = FALSE)
  )

  expect_s3_class(captured$error, "normalyzerde_error")
  expect_gt(length(captured$messages), 0)
  expect_true(any(grepl("NormalyzerDE", captured$messages, fixed = TRUE)))
  expect_true(any(grepl("version", captured$messages, ignore.case = TRUE)))
})

test_that("normalyzerDE stays quiet when inputs missing and quiet=TRUE", {
  captured <- nd_capture_conditions(
    normalyzerDE(jobName = "missing_inputs", quiet = TRUE)
  )

  expect_s3_class(captured$error, "normalyzerde_error")
  expect_length(captured$messages, 0)
})

test_that("normalyzer emits step messages during a successful quiet=FALSE run", {
  dataPath <- system.file(package = "NormalyzerDE", "extdata", "tiny_data.tsv")
  designPath <- system.file(
    package = "NormalyzerDE",
    "extdata",
    "tiny_design.tsv"
  )

  tmpDir <- withr::local_tempdir(pattern = "normalyzer_messages_")
  jobName <- "entry_messages"
  expectedDir <- file.path(tmpDir, NormalyzerDE:::sanitizeJobName(jobName))

  captured <- nd_capture_conditions(
    suppressWarnings(normalyzer(
      jobName = jobName,
      designPath = designPath,
      dataPath = dataPath,
      outputDir = tmpDir,
      normalizeRetentionTime = FALSE,
      skipAnalysis = FALSE,
      quiet = FALSE,
      writeReportAsPngs = TRUE
    ))
  )

  expect_null(captured$error)
  expect_null(captured$value)
  nd_expect_messages(
    captured$messages,
    c(
      "\\[Step 1/5\\].*Load data",
      "\\[Step 2/5\\].*Performing normalizations",
      "\\[Step 3/5\\].*Generating evaluation measures",
      "\\[Step 4/5\\].*written",
      "\\[Step 5/5\\].*generated"
    )
  )

  pngDir <- file.path(expectedDir, "pngs")
  expect_true(dir.exists(pngDir))
  expect_true(file.exists(file.path(pngDir, "1_front_page.png")))
})

test_that("normalyzerDE emits step messages during a successful quiet=FALSE run", {
  dataPath <- system.file(package = "NormalyzerDE", "extdata", "tiny_data.tsv")
  designPath <- system.file(
    package = "NormalyzerDE",
    "extdata",
    "tiny_design.tsv"
  )

  tmpDir <- withr::local_tempdir(pattern = "normalyzerde_messages_")
  jobName <- "entry_de_messages"
  expectedDir <- file.path(tmpDir, NormalyzerDE:::sanitizeJobName(jobName))

  captured <- nd_capture_conditions(
    suppressWarnings(normalyzerDE(
      jobName = jobName,
      comparisons = "4-5",
      designPath = designPath,
      dataPath = dataPath,
      outputDir = tmpDir,
      condCol = "group",
      quiet = FALSE,
      writeReportAsPngs = TRUE
    ))
  )

  expect_null(captured$error)
  expect_null(captured$value)
  nd_expect_messages(
    captured$messages,
    c(
      "\\[Step 1/6\\].*Load data",
      "\\[Step 2/6\\].*technical replicate reduction",
      "\\[Step 3/6\\].*Statistics object prepared",
      "\\[Step 4/6\\].*Contrast calculations done",
      "\\[Step 6/6\\].*Statistics report written"
    )
  )

  pngDir <- file.path(expectedDir, "de_pngs")
  expect_true(dir.exists(pngDir))
  expect_true(file.exists(file.path(pngDir, "1_front.png")))
})
