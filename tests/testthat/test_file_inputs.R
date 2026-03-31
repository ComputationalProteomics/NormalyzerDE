context("file inputs")

nd_example_file_paths <- function(data_name) {
  list(
    dataPath = system.file("extdata", data_name, package = "NormalyzerDE"),
    designPath = system.file(
      "extdata",
      "tiny_design.tsv",
      package = "NormalyzerDE"
    )
  )
}

nd_local_job_paths <- function(pattern, job_name) {
  outDir <- withr::local_tempdir(
    pattern = pattern,
    .local_envir = parent.frame()
  )
  list(
    outDir = outDir,
    expectedDir = file.path(outDir, NormalyzerDE:::sanitizeJobName(job_name))
  )
}

nd_write_excel_design_export <- function(design_path, tmp_dir) {
  excel_design_path <- file.path(tmp_dir, "tiny_design_excel.tsv")
  writeLines(
    paste0(readLines(design_path, warn = FALSE), "\t\t"),
    excel_design_path
  )
  excel_design_path
}

test_that("normalyzer can read default input files", {
  paths <- nd_example_file_paths("tiny_data.tsv")
  expect_true(all(nzchar(unlist(paths, use.names = FALSE))))
  jobName <- paste0("file_norm_", sample.int(1e9, 1))
  jobPaths <- nd_local_job_paths("file_norm_", jobName)

  out <- suppressWarnings(normalyzer(
    jobName = jobName,
    designPath = paths$designPath,
    dataPath = paths$dataPath,
    outputDir = jobPaths$outDir,
    normalizeRetentionTime = FALSE,
    skipAnalysis = TRUE,
    quiet = TRUE
  ))

  expect_null(out)
  expect_true(dir.exists(jobPaths$expectedDir))
  expect_true(file.exists(file.path(
    jobPaths$expectedDir,
    "submitted_rawdata.txt"
  )))
  expect_gt(
    length(list.files(
      jobPaths$expectedDir,
      pattern = "-normalized\\.txt$",
      full.names = TRUE
    )),
    0
  )
})

test_that("normalyzerDE can read default input files", {
  paths <- nd_example_file_paths("tiny_data_log2.tsv")
  expect_true(all(nzchar(unlist(paths, use.names = FALSE))))
  jobName <- paste0("file_de_", sample.int(1e9, 1))
  jobPaths <- nd_local_job_paths("file_de_", jobName)

  out <- suppressWarnings(normalyzerDE(
    jobName = jobName,
    comparisons = "4-5",
    designPath = paths$designPath,
    dataPath = paths$dataPath,
    outputDir = jobPaths$outDir,
    type = "limma",
    logTrans = FALSE,
    quiet = TRUE
  ))

  expect_null(out)

  outStatsPath <- file.path(
    jobPaths$expectedDir,
    paste0(basename(jobPaths$expectedDir), "_stats.tsv")
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
  paths <- nd_example_file_paths("tiny_data.tsv")
  expect_true(all(nzchar(unlist(paths, use.names = FALSE))))
  jobName <- "file_norm_excel_design"
  jobPaths <- nd_local_job_paths("file_norm_excel_design_", jobName)
  excelDesignPath <- nd_write_excel_design_export(
    design_path = paths$designPath,
    tmp_dir = jobPaths$outDir
  )

  out <- suppressWarnings(normalyzer(
    jobName = jobName,
    designPath = excelDesignPath,
    dataPath = paths$dataPath,
    outputDir = jobPaths$outDir,
    normalizeRetentionTime = FALSE,
    skipAnalysis = TRUE,
    quiet = TRUE
  ))

  expect_null(out)
  expect_true(file.exists(file.path(
    jobPaths$expectedDir,
    "submitted_rawdata.txt"
  )))
})

test_that("normalyzerDE can read Excel-exported design files with blank tab columns", {
  paths <- nd_example_file_paths("tiny_data_log2.tsv")
  expect_true(all(nzchar(unlist(paths, use.names = FALSE))))
  jobName <- "file_de_excel_design"
  jobPaths <- nd_local_job_paths("file_de_excel_design_", jobName)
  excelDesignPath <- nd_write_excel_design_export(
    design_path = paths$designPath,
    tmp_dir = jobPaths$outDir
  )

  out <- suppressWarnings(normalyzerDE(
    jobName = jobName,
    comparisons = "4-5",
    designPath = excelDesignPath,
    dataPath = paths$dataPath,
    outputDir = jobPaths$outDir,
    type = "limma",
    logTrans = FALSE,
    quiet = TRUE
  ))

  expect_null(out)
  expect_true(file.exists(file.path(
    jobPaths$expectedDir,
    paste0(basename(jobPaths$expectedDir), "_stats.tsv")
  )))
})
