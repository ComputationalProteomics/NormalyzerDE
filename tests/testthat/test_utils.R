context("utils.R")

test_that("getIndexList", {
  out <- getIndexList(c(1, 2, 3))
  expect_true(all.equal(out, list("1" = 1, "2" = 2, "3" = 3)))

  out <- getIndexList(c(1, 2, 1, 2, 3))
  expect_true(all.equal(out, list("1" = c(1, 3), "2" = c(2, 4), "3" = 5)))
})

test_that("getRowNAFilterContrast", {
  data("example_stat_data")
  data("example_design")
  levels <- example_design$group
  out <- getRowNAFilterContrast(head(example_stat_data), levels, minCount = 3)
  expect_true(
    all.equal(
      out,
      c(
        "1" = TRUE,
        "2" = TRUE,
        "3" = FALSE,
        "4" = FALSE,
        "5" = TRUE,
        "6" = TRUE
      )
    )
  )
})

test_that("getReplicateSortedData_constant", {
  rawMat <- matrix(
    c(
      1,
      1,
      1,
      2,
      2,
      2,
      3,
      3,
      3,
      4,
      4,
      4
    ),
    nrow = 3
  )
  colnames(rawMat) <- c("A", "A", "B", "B")

  groups <- c("A", "A", "B", "B")
  sortedMat <- getReplicateSortedData(rawMat, groups)

  expect_true(
    all.equal(
      rawMat,
      sortedMat
    )
  )
})

test_that("getReplicateSortedData_reordering", {
  rawMat <- matrix(
    c(
      3,
      3,
      3,
      1,
      1,
      1,
      4,
      4,
      4,
      2,
      2,
      2
    ),
    nrow = 3
  )
  colnames(rawMat) <- c("B", "A", "B", "A")

  expectedMat <- matrix(
    c(
      1,
      1,
      1,
      2,
      2,
      2,
      3,
      3,
      3,
      4,
      4,
      4
    ),
    nrow = 3
  )
  colnames(expectedMat) <- c("A", "A", "B", "B")

  groups <- c("B", "A", "B", "A")
  sortedMat <- getReplicateSortedData(rawMat, groups)

  expect_true(
    all.equal(
      sortedMat,
      expectedMat
    )
  )
})


test_that("filterLowRep", {
  test_data <- data.frame(
    c(NA, 1, 1, NA),
    c(1, NA, 1, NA),
    c(3, 3, 3, 0),
    c(5, 3, NA, 0),
    c(NA, 5, NA, 0),
    c(7, 9, NA, 0)
  )
  colnames(test_data) <- c("a1", "a2", "a3", "b1", "b2", "b3")
  groups <- c(rep("A", 3), rep("B", 3))

  expected_out_data <- data.frame(
    "a1" = c(NA, 1),
    "a2" = c(1, NA),
    "a3" = c(3, 3),
    "b1" = c(5, 3),
    "b2" = c(NA, 5),
    "b3" = c(7, 9)
  )

  out <- filterLowRep(test_data, groups, leastRep = 2)

  expect_true(
    all.equal(
      expected_out_data,
      out
    )
  )

  out2 <- filterLowRep(test_data, groups, leastRep = 0)

  expect_true(
    all.equal(
      test_data,
      out2
    )
  )
})


test_that("imputeGroupValues", {
  test_data <- data.frame(
    c(NA, 1, 1, NA),
    c(1, NA, 1, NA),
    c(3, 3, 3, NA),
    c(5, 3, NA, NA),
    c(NA, 5, NA, NA),
    c(NA, 5, NA, NA),
    c(NA, 5, NA, NA),
    c(7, 9, NA, 1)
  )
  colnames(test_data) <- c("a1", "a2", "a3", "b1", "b2", "b3", "c1", "c2")
  groups <- c(rep("A", 3), rep("B", 3), rep("C", 2))

  expected_out_data <- data.frame(
    "a1" = c(NA, 1, 1, NA),
    "a2" = c(1, NA, 1, NA),
    "a3" = c(3, 3, 3, NA),
    "b1" = c(5, 3, 1, NA),
    "b2" = c(NA, 5, NA, NA),
    "b3" = c(NA, 5, NA, NA),
    "c1" = c(NA, 5, 1, NA),
    "c2" = c(7, 9, NA, 1)
  )

  expected_out_data2 <- data.frame(
    "a1" = c(NA, 1, 1, 1),
    "a2" = c(1, NA, 1, NA),
    "a3" = c(3, 3, 3, NA),
    "b1" = c(5, 3, 1, 1),
    "b2" = c(NA, 5, NA, NA),
    "b3" = c(NA, 5, NA, NA),
    "c1" = c(NA, 5, 1, NA),
    "c2" = c(7, 9, NA, 1)
  )

  out <- imputeGroupValues(test_data, groups, minFraction = 1)

  expect_true(
    all.equal(
      expected_out_data,
      out
    )
  )

  out2 <- imputeGroupValues(test_data, groups, minFraction = 0.25)

  expect_true(
    all.equal(
      expected_out_data2,
      out2
    )
  )
})

test_that("imputeGroupValues_ungrouped_samples", {
  test_data <- data.frame(
    "s1" = c(NA),
    "s2" = c(2),
    "s3" = c(NA),
    "s4" = c(3)
  )
  groups <- c("B", "A", "B", "A")

  out <- imputeGroupValues(test_data, groups, minFraction = 1)

  expect_true(
    all.equal(
      out,
      data.frame("s1" = c(2), "s2" = c(2), "s3" = c(NA_real_), "s4" = c(3))
    )
  )
})

test_that("imputeGroupValues ignores unused factor levels", {
  test_data <- data.frame(
    "s1" = c(NA),
    "s2" = c(2),
    "s3" = c(NA),
    "s4" = c(3)
  )
  groups <- factor(c("B", "A", "B", "A"), levels = c("A", "B", "C"))

  out <- imputeGroupValues(test_data, groups, minFraction = 1)

  expect_true(
    all.equal(
      out,
      data.frame("s1" = c(2), "s2" = c(2), "s3" = c(NA_real_), "s4" = c(3))
    )
  )
})

test_that("setupJobDir sanitizes jobName", {
  parentDir <- tempfile()
  dir.create(parentDir)
  on.exit(unlink(parentDir, recursive = TRUE), add = TRUE)

  jobDir <- expect_warning(
    setupJobDir("../my job", parentDir),
    class = "normalyzerde_warning"
  )

  expect_true(startsWith(normalizePath(jobDir), normalizePath(parentDir)))
  expect_false(grepl("[/\\\\]", basename(jobDir)))
  expect_false(grepl("\\.\\.", basename(jobDir)))
})

test_that("sanitizeJobName rejects invalid names before and after sanitization", {
  expect_error(sanitizeJobName(character()), class = "normalyzerde_error")
  expect_error(sanitizeJobName(NA_character_), class = "normalyzerde_error")
  expect_error(sanitizeJobName("..."), class = "normalyzerde_error")
  expect_equal(sanitizeJobName("analysis_run"), "analysis_run")
})

test_that("createDirectory errors for files and impossible nested paths", {
  parentDir <- withr::local_tempdir(pattern = "create_dir_")

  existingFile <- file.path(parentDir, "existing.txt")
  writeLines("x", existingFile)
  expect_error(createDirectory(existingFile), class = "normalyzerde_error")

  blockingFile <- file.path(parentDir, "blocking")
  writeLines("x", blockingFile)
  expect_error(
    createDirectory(file.path(blockingFile, "child")),
    class = "normalyzerde_error"
  )

  existingDir <- file.path(parentDir, "already_there")
  dir.create(existingDir)
  expect_no_error(createDirectory(existingDir))
})

test_that("elapsedSecondsBetweenSystimes returns a difftime", {
  start <- as.POSIXct("2024-01-01 00:00:00", tz = "UTC")
  end <- as.POSIXct("2024-01-01 00:00:42", tz = "UTC")

  elapsed <- elapsedSecondsBetweenSystimes(start, end)

  expect_s3_class(elapsed, "difftime")
  expect_equal(as.numeric(elapsed, units = "secs"), 42)
})

test_that("setupTestData returns named columns and bounded RT values", {
  set.seed(1)
  out <- setupTestData(
    nSamples = 3,
    nFeatures = 4,
    rtMin = 10,
    rtMax = 20,
    mean = 5,
    sd = 0
  )

  expect_equal(dim(out), c(4L, 5L))
  expect_equal(colnames(out), c("feature", "RT", "S1", "S2", "S3"))
  expect_equal(out$feature, paste0("feature_", seq_len(4)))
  expect_true(all(out$RT >= 10 & out$RT <= 20))
  expect_true(all(as.matrix(out[, c("S1", "S2", "S3")]) == 5))
})

test_that("setupJobDir reuses an existing empty directory", {
  parentDir <- tempfile()
  dir.create(parentDir)
  on.exit(unlink(parentDir, recursive = TRUE), add = TRUE)

  jobName <- "my_job"
  existingDir <- file.path(parentDir, jobName)
  dir.create(existingDir)

  jobDir <- expect_silent(setupJobDir(jobName, parentDir))

  expect_true(dir.exists(jobDir))
  expect_identical(normalizePath(jobDir), normalizePath(existingDir))
})

test_that("setupJobDir errors when reusing a non-empty directory by default", {
  parentDir <- tempfile()
  dir.create(parentDir)
  on.exit(unlink(parentDir, recursive = TRUE), add = TRUE)

  existingDir <- file.path(parentDir, "my_job")
  dir.create(existingDir)
  writeLines("old output", file.path(existingDir, "stale.txt"))

  expect_error(
    setupJobDir("my_job", parentDir),
    class = "normalyzerde_error"
  )
})

test_that("setupJobDir can reuse a non-empty directory when explicitly requested", {
  parentDir <- tempfile()
  dir.create(parentDir)
  on.exit(unlink(parentDir, recursive = TRUE), add = TRUE)

  existingDir <- file.path(parentDir, "my_job")
  dir.create(existingDir)
  writeLines("old output", file.path(existingDir, "stale.txt"))

  jobDir <- expect_warning(
    setupJobDir("my_job", parentDir, reuseOutputDir = TRUE),
    class = "normalyzerde_warning"
  )

  expect_identical(normalizePath(jobDir), normalizePath(existingDir))
})

test_that("setupJobDir defaults to the working directory when outputDir is NULL", {
  workDir <- withr::local_tempdir(pattern = "setup_job_dir_")

  withr::with_dir(workDir, {
    jobDir <- expect_no_error(setupJobDir("default_dir_job", NULL))
    expect_identical(
      normalizePath(jobDir),
      normalizePath(file.path(workDir, "default_dir_job"))
    )
  })
})
