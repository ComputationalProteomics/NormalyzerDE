context("file inputs")

test_that("normalyzer can read default input files", {
  dataPath <- system.file("extdata", "tiny_data.tsv", package = "NormalyzerDE")
  designPath <- system.file(
    "extdata",
    "tiny_design.tsv",
    package = "NormalyzerDE"
  )
  expect_true(nzchar(dataPath))
  expect_true(nzchar(designPath))

  jobName <- paste0("file_norm_", sample.int(1e9, 1))
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

test_that("normalyzerDE can read default input files", {
  dataPath <- system.file(
    "extdata",
    "tiny_data_log2.tsv",
    package = "NormalyzerDE"
  )
  designPath <- system.file(
    "extdata",
    "tiny_design.tsv",
    package = "NormalyzerDE"
  )
  expect_true(nzchar(dataPath))
  expect_true(nzchar(designPath))

  jobName <- paste0("file_de_", sample.int(1e9, 1))
  outDir <- tempdir()
  expectedDir <- file.path(outDir, NormalyzerDE:::sanitizeJobName(jobName))
  if (dir.exists(expectedDir)) {
    unlink(expectedDir, recursive = TRUE, force = TRUE)
  }
  on.exit(unlink(expectedDir, recursive = TRUE, force = TRUE), add = TRUE)

  out <- suppressWarnings(normalyzerDE(
    jobName = jobName,
    comparisons = "4-5",
    designPath = designPath,
    dataPath = dataPath,
    outputDir = outDir,
    type = "limma",
    logTrans = FALSE,
    quiet = TRUE
  ))

  expect_null(out)

  outStatsPath <- file.path(
    expectedDir,
    paste0(basename(expectedDir), "_stats.tsv")
  )
  expect_true(file.exists(outStatsPath))

  outDf <- utils::read.table(
    outStatsPath,
    sep = "\t",
    header = TRUE,
    check.names = FALSE
  )
  expect_true(all(
    c("4-5_PValue", "4-5_AdjPVal", "4-5_log2FoldChange") %in% colnames(outDf)
  ))
})

test_that("normalyzer can read Excel-exported design files with blank tab columns", {
  dataPath <- system.file("extdata", "tiny_data.tsv", package = "NormalyzerDE")
  designPath <- system.file(
    "extdata",
    "tiny_design.tsv",
    package = "NormalyzerDE"
  )

  tmpDir <- withr::local_tempdir(pattern = "file_norm_excel_design_")
  excelDesignPath <- file.path(tmpDir, "tiny_design_excel.tsv")
  writeLines(paste0(readLines(designPath, warn = FALSE), "\t\t"), excelDesignPath)

  jobName <- "file_norm_excel_design"
  expectedDir <- file.path(tmpDir, NormalyzerDE:::sanitizeJobName(jobName))

  out <- suppressWarnings(normalyzer(
    jobName = jobName,
    designPath = excelDesignPath,
    dataPath = dataPath,
    outputDir = tmpDir,
    normalizeRetentionTime = FALSE,
    skipAnalysis = TRUE,
    quiet = TRUE
  ))

  expect_null(out)
  expect_true(file.exists(file.path(expectedDir, "submitted_rawdata.txt")))
})

test_that("normalyzerDE can read Excel-exported design files with blank tab columns", {
  dataPath <- system.file(
    "extdata",
    "tiny_data_log2.tsv",
    package = "NormalyzerDE"
  )
  designPath <- system.file(
    "extdata",
    "tiny_design.tsv",
    package = "NormalyzerDE"
  )

  tmpDir <- withr::local_tempdir(pattern = "file_de_excel_design_")
  excelDesignPath <- file.path(tmpDir, "tiny_design_excel.tsv")
  writeLines(paste0(readLines(designPath, warn = FALSE), "\t\t"), excelDesignPath)

  jobName <- "file_de_excel_design"
  expectedDir <- file.path(tmpDir, NormalyzerDE:::sanitizeJobName(jobName))

  out <- suppressWarnings(normalyzerDE(
    jobName = jobName,
    comparisons = "4-5",
    designPath = excelDesignPath,
    dataPath = dataPath,
    outputDir = tmpDir,
    type = "limma",
    logTrans = FALSE,
    quiet = TRUE
  ))

  expect_null(out)
  expect_true(file.exists(file.path(
    expectedDir,
    paste0(basename(expectedDir), "_stats.tsv")
  )))
})
