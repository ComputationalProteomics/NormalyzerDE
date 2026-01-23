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
    tinyRunThres = 50
  )

  expect_true(singleReplicateRun(nds))
})
