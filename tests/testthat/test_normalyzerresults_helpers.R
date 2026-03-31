context("NormalyzerResults helpers")

nd_make_results_dataset <- function(include_rt = FALSE, tiny_run_thres = 50) {
  raw <- matrix(
    c(
      10,
      11,
      12,
      13,
      11,
      12,
      13,
      14,
      12,
      13,
      14,
      15,
      13,
      14,
      15,
      16,
      14,
      15,
      16,
      17
    ),
    nrow = 5,
    byrow = TRUE
  )
  colnames(raw) <- paste0("s", seq_len(ncol(raw)))

  design <- data.frame(
    sample = colnames(raw),
    group = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  annot_df <- data.frame(
    feature = paste0("f", seq_len(nrow(raw))),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (include_rt) {
    annot_df$RT <- seq_len(nrow(raw))
  }

  NormalyzerDataset(
    jobName = if (include_rt) "results_rt" else "results_basic",
    designMatrix = design,
    rawData = raw,
    annotationData = as.matrix(annot_df),
    sampleNameCol = "sample",
    groupNameCol = "group",
    tinyRunThres = tiny_run_thres,
    quiet = TRUE
  )
}

test_that("NormalyzerResults constructor and nds setter round-trip", {
  nds_a <- nd_make_results_dataset(include_rt = FALSE)
  nds_b <- nd_make_results_dataset(include_rt = TRUE)

  nr <- NormalyzerResults(nds_a)
  nr <- `nds<-`(nr, nds_b)

  expect_identical(nds(nr), nds_b)
})

test_that("performNormalizations skips VSN for tiny runs", {
  nr <- NormalyzerResults(nd_make_results_dataset(include_rt = FALSE))

  out <- expect_no_error(
    suppressMessages(suppressWarnings(
      performNormalizations(nr, rtNorm = FALSE, quiet = FALSE)
    ))
  )

  expect_true("log2" %in% names(normalizations(out)))
  expect_false("VSN" %in% names(normalizations(out)))
})

test_that("performNormalizations skips VSN on no-log workflows without RT normalization", {
  nr <- NormalyzerResults(nd_make_results_dataset(
    include_rt = FALSE,
    tiny_run_thres = 1
  ))

  out <- expect_no_error(
    suppressMessages(suppressWarnings(
      performNormalizations(
        nr,
        rtNorm = FALSE,
        rtWindowMinCount = 1,
        noLogTransform = TRUE,
        quiet = FALSE
      )
    ))
  )

  expect_true(all(
    c("GI", "median", "mean", "Quantile", "CycLoess", "RLR") %in%
      names(normalizations(out))
  ))
  expect_false("VSN" %in% names(normalizations(out)))
})

test_that("performNormalizations skips RT-VSN for log-scale RT workflows", {
  nr <- NormalyzerResults(nd_make_results_dataset(
    include_rt = TRUE,
    tiny_run_thres = 1
  ))

  out <- expect_no_error(
    suppressMessages(suppressWarnings(
      performNormalizations(
        nr,
        rtNorm = TRUE,
        rtWindowMinCount = 1,
        noLogTransform = TRUE,
        quiet = FALSE
      )
    ))
  )

  expect_true(all(
    c("RT-median", "RT-mean", "RT-Loess") %in% names(normalizations(out))
  ))
  expect_false("RT-VSN" %in% names(normalizations(out)))
})
