context("NormalyzerDataset run modes")

test_that("checkSingleReplicateRun detects non-replicated groups", {
  mat <- matrix(stats::rnorm(60 * 3, mean = 10, sd = 1), nrow = 60)
  colnames(mat) <- c("s1", "s2", "s3")

  design <- data.frame(
    sample = colnames(mat),
    group = c("A", "A", "B"),
    stringsAsFactors = FALSE
  )
  rownames(design) <- design$sample

  se <- SummarizedExperiment::SummarizedExperiment(
    assay = mat,
    colData = design,
    rowData = data.frame(feature = paste0("f", seq_len(nrow(mat))))
  )
  S4Vectors::metadata(se) <- list(sample = "sample", group = "group")

  nds <- getVerifiedNormalyzerObject(
    jobName = "single_rep",
    summarizedExp = se,
    threshold = 0,
    omitSamples = FALSE,
    requireReplicates = FALSE,
    quiet = TRUE,
    noLogTransform = TRUE,
    tinyRunThres = 50
  )

  expect_true(singleReplicateRun(nds))
})

test_that("checkSingleReplicateRun detects singleton sample groups", {
  mat <- matrix(stats::rnorm(60 * 3, mean = 10, sd = 1), nrow = 60)
  colnames(mat) <- c("s1", "s2", "s3")

  design <- data.frame(
    sample = colnames(mat),
    group = c("A", "A", "A"),
    stringsAsFactors = FALSE
  )
  rownames(design) <- design$sample

  se <- SummarizedExperiment::SummarizedExperiment(
    assay = mat,
    colData = design,
    rowData = data.frame(feature = paste0("f", seq_len(nrow(mat))))
  )
  S4Vectors::metadata(se) <- list(sample = "sample", group = "group")

  nds <- getVerifiedNormalyzerObject(
    jobName = "singleton_group",
    summarizedExp = se,
    threshold = 0,
    omitSamples = FALSE,
    requireReplicates = FALSE,
    quiet = TRUE,
    noLogTransform = TRUE,
    tinyRunThres = 50
  )

  expect_true(singleReplicateRun(nds))
})

test_that("NormalyzerDataset stores sample names and reports tiny runs with RT data", {
  raw <- matrix(c(1, 3, 2, 4), nrow = 2)
  colnames(raw) <- c("s1", "s2")

  design <- data.frame(
    sample = c("s1", "s2"),
    group = c("A", "A"),
    stringsAsFactors = FALSE
  )

  annot <- as.matrix(data.frame(
    feature = c("f1", "f2"),
    RT = c(10, 20),
    check.names = FALSE
  ))

  nds <- expect_message(
    NormalyzerDataset(
      jobName = "tiny_run",
      designMatrix = design,
      rawData = raw,
      annotationData = annot,
      sampleNameCol = "sample",
      groupNameCol = "group",
      tinyRunThres = 50,
      quiet = FALSE
    )
  )

  expect_true(isTinyRun(nds))
  expect_equal(rawData(nds), raw)
  expect_equal(sampleNames(nds), c("s1", "s2"))
  expect_equal(retentionTimes(nds), c(10, 20))
  expect_true(singleReplicateRun(nds))
})

test_that("NormalyzerDataset helpers report replicate and RT edge cases", {
  raw <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2)
  colnames(raw) <- c("s1", "s2", "s3")

  design <- data.frame(
    sample = c("s1", "s2", "s3"),
    group = c("A", "A", "B"),
    stringsAsFactors = FALSE
  )

  annot_no_rt <- as.matrix(data.frame(feature = c("f1", "f2"), check.names = FALSE))
  annot_multi_rt <- as.matrix(data.frame(
    feature = c("f1", "f2"),
    RT = c(10, 20),
    `Observed RT` = c(11, 21),
    check.names = FALSE
  ))

  nds <- NormalyzerDataset(
    jobName = "nonreplicated",
    designMatrix = design,
    rawData = raw,
    annotationData = annot_no_rt,
    sampleNameCol = "sample",
    groupNameCol = "group",
    tinyRunThres = 0,
    quiet = TRUE
  )

  expect_message(detectSingleReplicate(nds, quiet = FALSE))
  expect_message(getRTColumn(annot_no_rt, quiet = FALSE))
  expect_null(getRTColumn(annot_no_rt, quiet = TRUE))
  expect_error(getRTColumn(annot_multi_rt, quiet = TRUE), class = "normalyzerde_error")
})
