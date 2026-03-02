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

  jobName <- paste0("prequant_byrow_", sample.int(1e9, 1))
  outDir <- tempdir()
  expectedDir <- file.path(outDir, NormalyzerDE:::sanitizeJobName(jobName))
  if (dir.exists(expectedDir)) {
    unlink(expectedDir, recursive = TRUE, force = TRUE)
  }
  on.exit(unlink(expectedDir, recursive = TRUE, force = TRUE), add = TRUE)

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
    `Protein.Names` = rep(paste0("Prot_", protein_ids), each = peptides_per_protein),
    check.names = FALSE
  )

  se <- SummarizedExperiment::SummarizedExperiment(
    assay = mat,
    colData = design,
    rowData = rowAnno
  )

  jobName <- paste0("prequant_protein_", sample.int(1e9, 1))
  outDir <- tempdir()
  expectedDir <- file.path(outDir, NormalyzerDE:::sanitizeJobName(jobName))
  if (dir.exists(expectedDir)) {
    unlink(expectedDir, recursive = TRUE, force = TRUE)
  }
  on.exit(unlink(expectedDir, recursive = TRUE, force = TRUE), add = TRUE)

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
      10, NA, NA, 10, NA, NA, # sparse row -> should be filtered for leastRepCount=2
      10, 10, 10, 11, 11, 11, # well observed
      9, 9, NA, 9, 9, NA, # some missing but >=2 observed per group
      NA, NA, NA, 8, 8, 8 # observed only in B
    ),
    nrow = 4,
    byrow = TRUE
  )
  colnames(raw) <- paste0("s", seq_len(ncol(raw)))

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
  rowIds <- as.character(seq_len(nrow(yQuant$E)))
  rownames(yQuant$E) <- rowIds
  rownames(yQuant$other$n.observations) <- rowIds
  rownames(yQuant$other$standard.error) <- rowIds
  colnames(yQuant$other$n.observations) <- colnames(yQuant$E)
  colnames(yQuant$other$standard.error) <- colnames(yQuant$E)

  rdsPath <- file.path(tempdir(), paste0("limpa_quant_", sample.int(1e9, 1), ".rds"))
  saveRDS(yQuant, file = rdsPath)
  on.exit(unlink(rdsPath), add = TRUE)

  se <- SummarizedExperiment::SummarizedExperiment(
    assay = yQuant$E,
    colData = design,
    rowData = data.frame(feature = paste0("f", seq_len(nrow(yQuant$E))))
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

