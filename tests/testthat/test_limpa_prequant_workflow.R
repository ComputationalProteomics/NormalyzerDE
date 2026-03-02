context("limpa pre-quant workflow")

test_that("normalyzer preQuant='limpa' writes quantified RDS (by-row)", {
  testthat::skip_if_not_installed("limpa")

  set.seed(1)
  mat <- matrix(stats::rnorm(20 * 4, mean = 10, sd = 1), nrow = 20)
  colnames(mat) <- paste0("s", seq_len(ncol(mat)))
  mat[sample.int(length(mat), 8)] <- NA_real_

  design <- data.frame(
    sample = colnames(mat),
    group = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  rownames(design) <- design$sample

  se <- SummarizedExperiment::SummarizedExperiment(
    assay = mat,
    colData = design,
    rowData = data.frame(`Protein.Group` = paste0("P", seq_len(nrow(mat))))
  )

  outDir <- withr::local_tempdir(pattern = "prequant_byrow_")
  jobName <- "prequant_byrow"
  expectedDir <- file.path(outDir, NormalyzerDE:::sanitizeJobName(jobName))

  out <- suppressWarnings(normalyzer(
    jobName = jobName,
    experimentObj = se,
    outputDir = outDir,
    preQuant = "limpa",
    limpaByRow = TRUE,
    limpaQuantArgs = list(chunk = 10L, verbose = FALSE),
    noLogTransform = TRUE,
    normalizeRetentionTime = FALSE,
    skipAnalysis = TRUE,
    quiet = TRUE,
    sampleAbundThres = 1,
    requireReplicates = FALSE
  ))

  expect_null(out)
  expect_true(file.exists(file.path(expectedDir, "log2-normalized.txt")))

  rdsPath <- file.path(
    expectedDir,
    paste0(basename(expectedDir), "_limpa_quantified.rds")
  )
  expect_true(file.exists(rdsPath))

  y <- readRDS(rdsPath)
  expect_true(inherits(y, "EList"))
  keepRows <- rowSums(!is.na(mat)) > 0
  expect_equal(nrow(y$E), sum(keepRows))
  expect_equal(ncol(y$E), ncol(mat))
  expect_false(anyNA(y$E))
})

test_that("normalyzer preQuant='limpa' can summarize peptides to proteins", {
  testthat::skip_if_not_installed("limpa")

  set.seed(1)
  n_proteins <- 5
  peptides_per_protein <- 2
  n_samples <- 4

  protein_ids <- paste0("P", seq_len(n_proteins))
  peptide_protein <- rep(protein_ids, each = peptides_per_protein)

  mat <- matrix(
    stats::rnorm(length(peptide_protein) * n_samples, mean = 10, sd = 0.5),
    nrow = length(peptide_protein)
  )
  colnames(mat) <- paste0("s", seq_len(ncol(mat)))
  mat[sample.int(length(mat), 6)] <- NA_real_

  design <- data.frame(
    sample = colnames(mat),
    group = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  rownames(design) <- design$sample

  rowAnno <- data.frame(
    `Protein.Group` = peptide_protein,
    `Protein.Names` = rep(
      paste0("Prot_", protein_ids),
      each = peptides_per_protein
    ),
    check.names = FALSE
  )

  se <- SummarizedExperiment::SummarizedExperiment(
    assay = mat,
    colData = design,
    rowData = rowAnno
  )

  outDir <- withr::local_tempdir(pattern = "prequant_protein_")
  jobName <- "prequant_protein"
  expectedDir <- file.path(outDir, NormalyzerDE:::sanitizeJobName(jobName))

  out <- suppressWarnings(normalyzer(
    jobName = jobName,
    experimentObj = se,
    outputDir = outDir,
    preQuant = "limpa",
    limpaProteinIdCol = "Protein.Group",
    limpaQuantArgs = list(chunk = 10L, verbose = FALSE),
    noLogTransform = TRUE,
    normalizeRetentionTime = FALSE,
    skipAnalysis = TRUE,
    quiet = TRUE,
    sampleAbundThres = 1,
    requireReplicates = FALSE
  ))

  expect_null(out)

  rdsPath <- file.path(
    expectedDir,
    paste0(basename(expectedDir), "_limpa_quantified.rds")
  )
  y <- readRDS(rdsPath)
  expect_true(inherits(y, "EList"))
  expect_equal(nrow(y$E), n_proteins)
  expect_true("Protein.Group" %in% colnames(y$genes))
})

test_that("calculateContrasts can reuse quantified EList from limpaQuantifiedRds", {
  testthat::skip_if_not_installed("limpa")

  raw <- matrix(
    c(
      10,
      NA,
      NA,
      10,
      NA,
      NA,
      10,
      10,
      10,
      11,
      11,
      11,
      9,
      9,
      NA,
      9,
      9,
      NA,
      NA,
      NA,
      NA,
      8,
      8,
      8
    ),
    nrow = 4,
    byrow = TRUE
  )
  colnames(raw) <- paste0("s", seq_len(ncol(raw)))
  rownames(raw) <- paste0("f", seq_len(nrow(raw)))

  design <- data.frame(
    sample = colnames(raw),
    group = c("A", "A", "A", "B", "B", "B"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  rownames(design) <- design$sample

  yQuant <- limpa::dpcQuantByRow(
    y = raw,
    chunk = 10L,
    verbose = FALSE
  )
  colnames(yQuant$other$n.observations) <- colnames(yQuant$E)
  colnames(yQuant$other$standard.error) <- colnames(yQuant$E)

  rdsPath <- withr::local_tempfile(pattern = "limpa_quant_", fileext = ".rds")
  saveRDS(yQuant, file = rdsPath)

  se <- SummarizedExperiment::SummarizedExperiment(
    assay = yQuant$E,
    colData = design,
    rowData = data.frame(feature = rownames(yQuant$E))
  )

  nst <- NormalyzerStatistics(se, logTrans = FALSE)
  out <- calculateContrasts(
    nst,
    comparisons = "A-B",
    condCol = "group",
    type = "limpa",
    leastRepCount = 2,
    limpaQuantifiedRds = rdsPath,
    limpaKeep = "elist",
    limpaPostQuantNorm = "median"
  )

  pvals <- pairwiseCompsP(out)[["A-B"]]
  expect_length(pvals, nrow(yQuant$E))
  expect_true(is.na(pvals[1]))
  expect_true(any(!is.na(pvals[-1])))

  backend <- backendData(out)[["limpa"]]
  expect_type(backend, "list")
  expect_equal(backend$quantifiedRds, rdsPath)

  yUsed <- backend$elists[[".global"]]
  expect_true(inherits(yUsed, "EList"))

  linearMedians <- matrixStats::colMedians(2**yUsed$E)
  expect_lt(max(linearMedians) - min(linearMedians), 1e-6)
})

test_that("normalyzerDE auto-detects quantified RDS next to Normalyzer output matrix", {
  testthat::skip_if_not_installed("limpa")

  set.seed(1)
  nFeatures <- 30
  raw <- matrix(
    stats::rnorm(nFeatures * 4, mean = 10, sd = 1),
    nrow = nFeatures
  )
  raw[1, ] <- c(10, NA, 15, NA)
  raw[2, ] <- c(11, NA, 11, 11)
  colnames(raw) <- paste0("s", seq_len(ncol(raw)))

  design <- data.frame(
    sample = colnames(raw),
    group = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  rownames(design) <- design$sample

  se <- SummarizedExperiment::SummarizedExperiment(
    assay = raw,
    colData = design,
    rowData = data.frame(
      feature = c("sparse", paste0("f", seq_len(nFeatures - 1))),
      check.names = FALSE
    )
  )

  outDir <- withr::local_tempdir(pattern = "prequant_autodetect_")
  jobNameNorm <- "prequant_autodetect"
  expectedDirNorm <- file.path(
    outDir,
    NormalyzerDE:::sanitizeJobName(jobNameNorm)
  )

  out <- suppressWarnings(normalyzer(
    jobName = jobNameNorm,
    experimentObj = se,
    outputDir = outDir,
    preQuant = "limpa",
    limpaByRow = TRUE,
    limpaQuantArgs = list(chunk = 10L, verbose = FALSE),
    noLogTransform = TRUE,
    normalizeRetentionTime = FALSE,
    skipAnalysis = TRUE,
    quiet = TRUE,
    sampleAbundThres = 1,
    requireReplicates = FALSE
  ))
  expect_null(out)

  log2Path <- file.path(expectedDirNorm, "log2-normalized.txt")
  expect_true(file.exists(log2Path))

  extraRds <- file.path(expectedDirNorm, "unrelated_limpa_quantified.rds")
  file.create(extraRds)

  designPath <- withr::local_tempfile(
    pattern = "design_autodetect_",
    fileext = ".tsv"
  )
  nd_write_table(design, designPath)

  jobNameDE <- "de_autodetect"
  expectedDirDE <- file.path(
    outDir,
    NormalyzerDE:::sanitizeJobName(jobNameDE)
  )

  outDE <- suppressWarnings(normalyzerDE(
    jobName = jobNameDE,
    comparisons = "A-B",
    designPath = designPath,
    dataPath = log2Path,
    outputDir = outDir,
    type = "limpa",
    leastRepCount = 2,
    limpaPostQuantNorm = "median",
    quiet = TRUE
  ))
  expect_null(outDE)

  outStatsPath <- file.path(
    expectedDirDE,
    paste0(basename(expectedDirDE), "_stats.tsv")
  )
  expect_true(file.exists(outStatsPath))

  outDf <- utils::read.table(
    outStatsPath,
    sep = "\t",
    header = TRUE,
    check.names = FALSE
  )
  pCol <- "A-B_PValue"
  expect_true(pCol %in% colnames(outDf))

  sparseRow <- outDf$feature == "sparse"
  expect_equal(sum(sparseRow), 1)
  expect_true(is.na(outDf[sparseRow, pCol]))
  expect_false(is.na(outDf[outDf$feature == "f1", pCol]))
})

test_that("autoDetectLimpaQuantifiedRds warns for single non-canonical cache file", {
  tmpDir <- withr::local_tempdir(pattern = "limpa_autodetect_warn_")
  file.create(file.path(tmpDir, "unrelated_limpa_quantified.rds"))

  out <- expect_warning(
    NormalyzerDE:::autoDetectLimpaQuantifiedRds(
      file.path(tmpDir, "log2-normalized.txt")
    ),
    "does not match the expected canonical filename"
  )
  expect_null(out)
})

test_that("normalyzerDE stops when multiple quantified RDS files are found", {
  tmpDir <- withr::local_tempdir(pattern = "multirds_")

  file.create(file.path(tmpDir, "a_limpa_quantified.rds"))
  file.create(file.path(tmpDir, "b_limpa_quantified.rds"))

  expect_error(
    normalyzerDE(
      jobName = "multirds_de",
      comparisons = "A-B",
      designPath = file.path(tmpDir, "design.tsv"),
      dataPath = file.path(tmpDir, "log2-normalized.txt"),
      type = "limpa",
      quiet = TRUE
    ),
    "Multiple '\\*\\_limpa\\_quantified\\.rds'"
  )
})

test_that("normalyzerDE prevents same-method double normalization for limpa", {
  expect_error(
    normalyzerDE(
      jobName = "dblnorm",
      comparisons = "A-B",
      designPath = "design.tsv",
      dataPath = file.path("some_dir", "median-normalized.txt"),
      type = "limpa",
      limpaPostQuantNorm = "median",
      quiet = TRUE
    ),
    "already be normalized"
  )
})

test_that("normalyzerDE warns for potential double normalization for limpa", {
  outDir <- withr::local_tempdir(pattern = "dblnorm_warn_")

  expect_warning(
    expect_error(
      normalyzerDE(
        jobName = "dblwarn",
        comparisons = "A-B",
        designPath = file.path(outDir, "design_missing.tsv"),
        dataPath = file.path(outDir, "median-normalized.tsv"),
        outputDir = outDir,
        type = "limpa",
        limpaPostQuantNorm = "GI",
        quiet = TRUE
      )
    ),
    "double-normalization"
  )
})
