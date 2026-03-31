context("limpa pre-quant workflow")

test_that("normalyzer preQuant='limpa' writes quantified RDS from limpaOptions helper (by-row)", {
  testthat::skip_if_not_installed("limpa")

  fixture <- nd_make_limpa_prequant_fixture()

  outDir <- withr::local_tempdir(pattern = "prequant_byrow_")
  jobName <- "prequant_byrow"
  expectedDir <- file.path(outDir, NormalyzerDE:::sanitizeJobName(jobName))

  out <- nd_run_prequant_limpa(
    se = fixture$se,
    out_dir = outDir,
    job_name = jobName,
    limpa_options = limpaOptions(
      byRow = TRUE,
      quantArgs = list(chunk = 10L, verbose = FALSE)
    )
  )

  expect_null(out)
  expect_true(file.exists(file.path(expectedDir, "log2-normalized.txt")))

  rdsPath <- file.path(
    expectedDir,
    paste0(basename(expectedDir), "_limpa_quantified.rds")
  )
  expect_true(file.exists(rdsPath))

  y <- readRDS(rdsPath)
  expect_true(inherits(y, "EList"))
  keepRows <- rowSums(!is.na(fixture$mat)) > 0
  expect_equal(nrow(y$E), sum(keepRows))
  expect_equal(ncol(y$E), ncol(fixture$mat))
  expect_false(anyNA(y$E))
})

test_that("normalyzer preQuant='limpa' can summarize peptides to proteins", {
  testthat::skip_if_not_installed("limpa")

  n_proteins <- 5
  peptides_per_protein <- 2

  protein_ids <- paste0("P", seq_len(n_proteins))
  peptide_protein <- rep(protein_ids, each = peptides_per_protein)

  rowAnno <- data.frame(
    `Protein.Group` = peptide_protein,
    `Protein.Names` = rep(
      paste0("Prot_", protein_ids),
      each = peptides_per_protein
    ),
    check.names = FALSE
  )

  fixture <- nd_make_limpa_prequant_fixture(
    n_features = length(peptide_protein),
    row_data = rowAnno,
    missing_count = 6L,
    sd = 0.5
  )

  outDir <- withr::local_tempdir(pattern = "prequant_protein_")
  jobName <- "prequant_protein"
  expectedDir <- file.path(outDir, NormalyzerDE:::sanitizeJobName(jobName))

  out <- nd_run_prequant_limpa(
    se = fixture$se,
    out_dir = outDir,
    job_name = jobName,
    limpa_options = limpaOptions(
      proteinIdCol = "Protein.Group",
      quantArgs = list(chunk = 10L, verbose = FALSE)
    )
  )

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

test_that("normalyzer preQuant='limpa' emits step messages during a successful quiet=FALSE run", {
  testthat::skip_if_not_installed("limpa")

  fixture <- nd_make_limpa_prequant_fixture()

  outDir <- withr::local_tempdir(pattern = "prequant_messages_")
  jobName <- "prequant_messages"

  captured <- nd_capture_conditions(
    nd_run_prequant_limpa(
      se = fixture$se,
      out_dir = outDir,
      job_name = jobName,
      quiet = FALSE,
      limpa_options = limpaOptions(
        byRow = TRUE,
        quantArgs = list(chunk = 10L, verbose = FALSE)
      )
    )
  )

  expect_null(captured$error)
  expect_null(captured$value)
  nd_expect_messages(
    captured$messages,
    c(
      "\\[Step 1/6\\].*Load data",
      "\\[Step 2/6\\].*Running limpa pre-quantification",
      "\\[Step 2/6\\].*limpa quantification completed",
      "\\[Step 3/6\\].*Performing normalizations",
      "\\[Step 4/6\\].*Skipped evaluation measures",
      "\\[Step 6/6\\].*Skipped plot generation"
    )
  )
})

test_that("normalyzer preQuant='limpa' reports when unique protein IDs skip summarization", {
  testthat::skip_if_not_installed("limpa")

  fixture <- nd_make_limpa_prequant_fixture(
    n_features = 12L,
    missing_count = 6L
  )

  outDir <- withr::local_tempdir(pattern = "prequant_unique_protein_ids_")
  jobName <- "prequant_unique_protein_ids"

  captured <- nd_capture_conditions(
    nd_run_prequant_limpa(
      se = fixture$se,
      out_dir = outDir,
      job_name = jobName,
      quiet = FALSE,
      limpa_options = limpaOptions(
        proteinIdCol = "Protein.Group",
        quantArgs = list(chunk = 10L, verbose = FALSE)
      )
    )
  )

  expect_null(captured$error)
  expect_null(captured$value)
  nd_expect_messages(
    captured$messages,
    "contains no duplicated identifiers"
  )

  rdsPath <- file.path(
    outDir,
    jobName,
    paste0(jobName, "_limpa_quantified.rds")
  )
  y <- readRDS(rdsPath)
  expect_equal(nrow(y$E), sum(rowSums(!is.na(fixture$mat)) > 0))
})

test_that("normalyzer preQuant='limpa' informs when supplied dpc overrides dpcMethod", {
  testthat::skip_if_not_installed("limpa")

  fixture <- nd_make_limpa_prequant_fixture(
    n_features = 12L,
    missing_count = 6L
  )

  outDir <- withr::local_tempdir(pattern = "prequant_dpc_override_")

  captured <- nd_capture_conditions(
    nd_run_prequant_limpa(
      se = fixture$se,
      out_dir = outDir,
      job_name = "prequant_dpc_override",
      quiet = FALSE,
      limpa_options = limpaOptions(
        byRow = TRUE,
        dpc = c(0, 1),
        dpcMethod = "dpc",
        quantArgs = list(chunk = 10L, verbose = TRUE)
      )
    )
  )

  expect_null(captured$error)
  expect_null(captured$value)
  nd_expect_messages(
    captured$messages,
    c(
      "ignoring .*dpcMethod",
      "Saved quantified EList"
    )
  )
})

test_that("calculateContrasts can reuse quantified EList from limpaOptions(quantifiedRds)", {
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

  design <- nd_make_design(c("A", "A", "A", "B", "B", "B"))

  yQuant <- limpa::dpcQuantByRow(
    y = raw,
    chunk = 10L,
    verbose = FALSE
  )
  colnames(yQuant$other$n.observations) <- colnames(yQuant$E)
  colnames(yQuant$other$standard.error) <- colnames(yQuant$E)

  rdsPath <- withr::local_tempfile(pattern = "limpa_quant_", fileext = ".rds")
  saveRDS(yQuant, file = rdsPath)

  se <- nd_make_summarized_experiment(
    assay = yQuant$E,
    groups = design$group,
    row_data = data.frame(feature = rownames(yQuant$E))
  )

  nst <- NormalyzerStatistics(se, logTrans = FALSE)
  out <- calculateContrasts(
    nst,
    comparisons = "A-B",
    condCol = "group",
    type = "limpa",
    leastRepCount = 2,
    limpaOptions = limpaOptions(
      quantifiedRds = rdsPath,
      keep = "elist",
      postQuantNorm = "median"
    )
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

test_that("calculateContrasts errors when quantifiedRds rows cannot be matched safely", {
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
      8,
      7,
      7,
      7,
      7,
      7,
      7,
      6,
      6,
      6,
      6,
      6,
      6
    ),
    nrow = 6,
    byrow = TRUE,
    dimnames = list(paste0("orig", seq_len(6)), paste0("s", seq_len(6)))
  )

  yQuant <- limpa::dpcQuantByRow(
    y = raw,
    chunk = 10L,
    verbose = FALSE
  )
  colnames(yQuant$other$n.observations) <- colnames(yQuant$E)
  colnames(yQuant$other$standard.error) <- colnames(yQuant$E)

  completed <- yQuant$E[6:1, , drop = FALSE]
  rownames(completed) <- paste0("different", seq_len(nrow(completed)))

  se <- nd_make_summarized_experiment(
    assay = completed,
    groups = c("A", "A", "A", "B", "B", "B"),
    row_data = data.frame(
      feature = paste0("feat", seq_len(nrow(completed))),
      row.names = rownames(completed),
      check.names = FALSE
    )
  )

  rdsPath <- withr::local_tempfile(pattern = "limpa_quant_", fileext = ".rds")
  saveRDS(yQuant, file = rdsPath)

  nst <- NormalyzerStatistics(se, logTrans = FALSE)
  expect_error(
    calculateContrasts(
      nst,
      comparisons = "A-B",
      condCol = "group",
      type = "limpa",
      leastRepCount = 1,
      limpaOptions = limpaOptions(quantifiedRds = rdsPath, keep = "elist")
    ),
    class = "normalyzerde_error"
  )
})

test_that("normalyzerDE ignores nearby quantified RDS files unless requested explicitly", {
  testthat::skip_if_not_installed("limpa")

  test_data <- matrix(
    c(
      10,
      11,
      10,
      13,
      12,
      11,
      NA,
      NA,
      NA,
      9,
      9,
      10,
      5,
      5,
      5,
      NA,
      NA,
      NA,
      7,
      8,
      7,
      7,
      7,
      7
    ),
    nrow = 4,
    byrow = TRUE
  )
  colnames(test_data) <- paste0("s", seq_len(ncol(test_data)))

  design <- nd_make_design(c("A", "A", "A", "B", "B", "B"))
  data_df <- data.frame(
    feature = paste0("f", seq_len(nrow(test_data))),
    as.data.frame(test_data, check.names = FALSE),
    check.names = FALSE
  )

  outDir <- withr::local_tempdir(pattern = "prequant_no_autodetect_")
  paths <- nd_write_data_and_design(
    tmp_dir = outDir,
    data = data_df,
    design = design,
    data_name = "limpa_data.tsv",
    design_name = "limpa_design.tsv"
  )

  file.create(file.path(outDir, "a_limpa_quantified.rds"))
  file.create(file.path(outDir, "b_limpa_quantified.rds"))

  jobNameImplicit <- "de_no_quantified_rds"
  expectedDirImplicit <- file.path(
    outDir,
    NormalyzerDE:::sanitizeJobName(jobNameImplicit)
  )

  outImplicit <- suppressWarnings(normalyzerDE(
    jobName = jobNameImplicit,
    comparisons = "A-B",
    designPath = paths$designPath,
    dataPath = paths$dataPath,
    outputDir = outDir,
    type = "limpa",
    logTrans = FALSE,
    leastRepCount = 1,
    limpaOptions = limpaOptions(
      byRow = TRUE,
      quantArgs = list(chunk = 10L, verbose = FALSE)
    ),
    quiet = TRUE
  ))
  expect_null(outImplicit)

  implicitStatsPath <- file.path(
    expectedDirImplicit,
    paste0(basename(expectedDirImplicit), "_stats.tsv")
  )
  expect_true(file.exists(implicitStatsPath))

  implicitDf <- utils::read.table(
    implicitStatsPath,
    sep = "\t",
    header = TRUE,
    check.names = FALSE
  )
  expect_true("A-B_PValue" %in% colnames(implicitDf))
})

test_that("normalyzerDE prevents same-method double normalization for limpa", {
  expect_error(
    normalyzerDE(
      jobName = "dblnorm",
      comparisons = "A-B",
      designPath = "design.tsv",
      dataPath = file.path("some_dir", "median-normalized.txt"),
      type = "limpa",
      limpaOptions = limpaOptions(postQuantNorm = "median"),
      quiet = TRUE
    ),
    class = "normalyzerde_error"
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
        limpaOptions = limpaOptions(postQuantNorm = "GI"),
        quiet = TRUE
      ),
      class = "normalyzerde_error"
    ),
    class = "normalyzerde_warning"
  )
})
