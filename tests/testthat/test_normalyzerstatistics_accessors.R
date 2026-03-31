context("NormalyzerStatistics accessors and helpers")

test_that("NormalyzerStatistics accessors and setters round-trip", {
  mat <- matrix(
    c(1, 2, 3, 4),
    nrow = 2,
    dimnames = list(c("f1", "f2"), c("s1", "s2"))
  )

  design <- data.frame(
    sample = c("s1", "s2"),
    group = c("A", "B"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  se <- SummarizedExperiment::SummarizedExperiment(
    assay = mat,
    colData = design,
    rowData = data.frame(id = rownames(mat), check.names = FALSE)
  )

  nst <- NormalyzerStatistics(se, logTrans = FALSE)

  expect_true(is.matrix(annotMat(nst)))
  expect_true(is.matrix(dataMat(nst)))
  expect_true(is.data.frame(designDf(nst)))
  expect_true(is.matrix(filteredDataMat(nst)))
  expect_true(is.logical(filteringContrast(nst)))

  expect_true(is.list(backendData(nst)))
  backendData(nst) <- list(example = 1)
  expect_equal(backendData(nst)$example, 1)

  expect_true(is.character(comparisons(nst)))
  nst <- `comparisons<-`(nst, c("A-B"))
  expect_equal(comparisons(nst), "A-B")

  expect_true(is.character(condCol(nst)))
  nst <- `condCol<-`(nst, c("A", "B"))
  expect_equal(condCol(nst), c("A", "B"))

  nst <- `contrastSplitter<-`(nst, "-")
  expect_equal(contrastSplitter(nst), "-")

  outMat <- dataMat(nst) * 2
  nst <- `dataMat<-`(nst, outMat)
  expect_equal(dataMat(nst), outMat)

  newDesign <- design
  newDesign$group <- c("A", "A")
  nst <- `designDf<-`(nst, newDesign)
  expect_equal(designDf(nst)$group, c("A", "A"))

  nst <- `pairwiseCompsP<-`(nst, list("A-B" = c(0.1, 0.2)))
  nst <- `pairwiseCompsFdr<-`(nst, list("A-B" = c(0.2, 0.4)))
  nst <- `pairwiseCompsAve<-`(nst, list("A-B" = c(1, 2)))
  nst <- `pairwiseCompsFold<-`(nst, list("A-B" = c(-1, 1)))

  expect_equal(pairwiseCompsP(nst)[["A-B"]], c(0.1, 0.2))
  expect_equal(pairwiseCompsFdr(nst)[["A-B"]], c(0.2, 0.4))
  expect_equal(pairwiseCompsAve(nst)[["A-B"]], c(1, 2))
  expect_equal(pairwiseCompsFold(nst)[["A-B"]], c(-1, 1))
})

test_that("chooseOneVsRestLabel selects an unused label", {
  expect_equal(chooseOneVsRestLabel(c("A", "B")), "rest")

  used <- c("rest", "others", "all_other", "all_others")
  expect_equal(chooseOneVsRestLabel(used), "rest1")
  expect_equal(chooseOneVsRestLabel(c(used, "rest1")), "rest2")

  expect_error(
    chooseOneVsRestLabel("A", candidates = character()),
    class = "normalyzerde_error"
  )
})

test_that("batchCol accessor preserves labels", {
  mat <- matrix(
    c(1, 2, 3, 4),
    nrow = 2,
    dimnames = list(c("f1", "f2"), c("s1", "s2"))
  )

  design <- data.frame(
    sample = c("s1", "s2"),
    group = c("A", "B"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  se <- SummarizedExperiment::SummarizedExperiment(
    assay = mat,
    colData = design,
    rowData = data.frame(id = rownames(mat), check.names = FALSE)
  )

  nst <- NormalyzerStatistics(se, logTrans = FALSE)
  nst <- `batchCol<-`(nst, c("batch1", "batch2"))

  expect_type(batchCol(nst), "character")
  expect_equal(batchCol(nst), c("batch1", "batch2"))
})

test_that("sanitizeLimmaDesign and calculateLimmaContrast work with coefMap", {
  set.seed(1)
  dataMat <- matrix(stats::rnorm(20), nrow = 5)
  colnames(dataMat) <- paste0("s", seq_len(ncol(dataMat)))

  designDf <- data.frame(
    group = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )

  limmaDesignRaw <- stats::model.matrix(~ 0 + Variable, data = transform(
    designDf,
    Variable = as.factor(group)
  ))

  limmaPrepared <- sanitizeLimmaDesign(limmaDesignRaw)
  limmaDesign <- limmaPrepared$design
  coefMap <- limmaPrepared$coefMap

  limmaFit <- limma::lmFit(dataMat, limmaDesign)

  out <- calculateLimmaContrast(
    dataMat,
    limmaDesign,
    limmaFit,
    levels = c("A", "B"),
    useIntensityTrend = FALSE,
    coefMap = coefMap
  )

  expect_true(is.list(out))
  expect_true(all(c("P", "FDR", "Ave", "Fold") %in% names(out)))
  expect_equal(length(out$P), nrow(dataMat))

  expect_error(
    calculateLimmaContrast(
      dataMat,
      limmaDesign,
      limmaFit,
      levels = c("A", "B"),
      useIntensityTrend = FALSE,
      coefMap = coefMap[1]
    ),
    class = "normalyzerde_error"
  )

  expect_error(
    calculateLimmaContrast(
      dataMat,
      limmaDesign,
      limmaFit,
      levels = c("A", "B"),
      useIntensityTrend = FALSE,
      coefMap = coefMap[c("VariableA")]
    ),
    class = "normalyzerde_error"
  )
})

test_that("contrast helper validators handle model and contrast edge cases", {
  design <- data.frame(
    group = c("A", "A", "B", "B"),
    batch = c("x", "y", "x", "y"),
    stringsAsFactors = FALSE
  )

  model_no_batch <- setupModelFromDesign(design, condCol = "group")
  mm_no_batch <- stats::model.matrix(model_no_batch)
  expect_true(all(c("VariableA", "VariableB") %in% colnames(mm_no_batch)))

  model_with_batch <- setupModelFromDesign(
    design,
    condCol = "group",
    batchCol = "batch",
    type = "limma"
  )
  mm_with_batch <- stats::model.matrix(model_with_batch)
  expect_true(all(c("VariableA", "VariableB") %in% colnames(mm_with_batch)))
  expect_true(any(grepl("^Batch", colnames(mm_with_batch))))

  expect_error(
    setupModelFromDesign(
      design,
      condCol = "group",
      batchCol = "batch",
      type = "welch"
    ),
    class = "normalyzerde_error"
  )

  expect_equal(parseContrastLevels("A-B", "-"), c("A", "B"))
  expect_error(
    parseContrastLevels("A-B-C", "-"),
    class = "normalyzerde_error"
  )

  expect_no_error(assertContrastLevelsPresent(c("A", "B"), c("A", "B", "C")))
  expect_error(
    assertContrastLevelsPresent(c("X", "B"), c("A", "B", "C")),
    class = "normalyzerde_error"
  )
  expect_error(
    assertContrastLevelsPresent(c("A", "Y"), c("A", "B", "C")),
    class = "normalyzerde_error"
  )

  expect_no_error(verifyContrasts(c("A", "B"), c("A-B")))
  expect_error(
    verifyContrasts(c("A", "B"), c("A")),
    class = "normalyzerde_error"
  )
  expect_error(
    verifyContrasts(c("A", "B"), c("A-C")),
    class = "normalyzerde_error"
  )
})

test_that(".getContrastSplitter falls back to defaults when needed", {
  expect_equal(.getContrastSplitter(list(), default = ":"), ":")

  mat <- matrix(
    c(1, 2, 3, 4),
    nrow = 2,
    dimnames = list(c("f1", "f2"), c("s1", "s2"))
  )
  design <- data.frame(
    sample = c("s1", "s2"),
    group = c("A", "B"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  se <- SummarizedExperiment::SummarizedExperiment(
    assay = mat,
    colData = design,
    rowData = data.frame(id = rownames(mat), check.names = FALSE)
  )

  nst <- NormalyzerStatistics(se, logTrans = FALSE)
  nst <- `contrastSplitter<-`(nst, "")

  expect_equal(.getContrastSplitter(nst, default = ":"), ":")
})
