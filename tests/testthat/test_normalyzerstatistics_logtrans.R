test_that("NormalyzerStatistics logTrans replaces non-finite values with NA", {
  mat <- matrix(
    c(
      1,
      0,
      2,
      4,
      -1,
      NA
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("feat1", "feat2"), c("s1", "s2", "s3"))
  )

  se <- SummarizedExperiment::SummarizedExperiment(
    assay = mat,
    rowData = data.frame(id = c("feat1", "feat2")),
    colData = data.frame(sample = c("s1", "s2", "s3"))
  )

  expect_warning(
    nst <- NormalyzerStatistics(se, logTrans = TRUE),
    "Non-finite values produced by log2 transform"
  )

  dm <- dataMat(nst)
  expect_equal(dm["feat1", "s1"], log2(1))
  expect_true(is.na(dm["feat1", "s2"]))
  expect_equal(dm["feat1", "s3"], log2(2))

  expect_equal(dm["feat2", "s1"], log2(4))
  expect_true(is.na(dm["feat2", "s2"]))
  expect_true(is.na(dm["feat2", "s3"]))
})

test_that("NormalyzerStatistics logTrans does not warn for existing NA values", {
  mat <- matrix(
    c(1, NA, 2),
    nrow = 1,
    dimnames = list("feat1", c("s1", "s2", "s3"))
  )
  se <- SummarizedExperiment::SummarizedExperiment(
    assay = mat,
    rowData = data.frame(id = "feat1"),
    colData = data.frame(sample = c("s1", "s2", "s3"))
  )

  nst <- expect_no_warning(NormalyzerStatistics(se, logTrans = TRUE))
  dm <- dataMat(nst)
  expect_equal(dm["feat1", "s1"], log2(1))
  expect_true(is.na(dm["feat1", "s2"]))
  expect_equal(dm["feat1", "s3"], log2(2))
})

test_that("NormalyzerStatistics keeps 0 when logTrans is FALSE", {
  mat <- matrix(
    c(0, 1),
    nrow = 1,
    dimnames = list("feat1", c("s1", "s2"))
  )
  se <- SummarizedExperiment::SummarizedExperiment(
    assay = mat,
    rowData = data.frame(id = "feat1"),
    colData = data.frame(sample = c("s1", "s2"))
  )

  nst <- NormalyzerStatistics(se, logTrans = FALSE)
  dm <- dataMat(nst)
  expect_equal(dm["feat1", "s1"], 0)
  expect_equal(dm["feat1", "s2"], 1)
})
