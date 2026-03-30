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

test_that("setupJobDir reuses existing directory by default", {
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

test_that("setupJobDir warns when reusing a non-empty directory", {
  parentDir <- tempfile()
  dir.create(parentDir)
  on.exit(unlink(parentDir, recursive = TRUE), add = TRUE)

  existingDir <- file.path(parentDir, "my_job")
  dir.create(existingDir)
  writeLines("old output", file.path(existingDir, "stale.txt"))

  jobDir <- expect_warning(
    setupJobDir("my_job", parentDir),
    class = "normalyzerde_warning"
  )

  expect_identical(normalizePath(jobDir), normalizePath(existingDir))
})
