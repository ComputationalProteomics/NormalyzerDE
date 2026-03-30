context("limpa backend")

test_that("resolveLimpaQuantByRowFn supports old and new limpa APIs", {
  new_ns <- new.env(parent = emptyenv())
  new_ns$dpcQuantByRow <- function(...) "new"

  old_ns <- new.env(parent = emptyenv())
  old_ns$dpcImpute <- function(...) "old"

  both_ns <- new.env(parent = emptyenv())
  both_ns$dpcQuantByRow <- function(...) "new"
  both_ns$dpcImpute <- function(...) "old"

  expect_identical(
    NormalyzerDE:::resolveLimpaQuantByRowFn(new_ns),
    new_ns$dpcQuantByRow
  )
  expect_identical(
    NormalyzerDE:::resolveLimpaQuantByRowFn(old_ns),
    old_ns$dpcImpute
  )
  expect_identical(
    NormalyzerDE:::resolveLimpaQuantByRowFn(both_ns),
    both_ns$dpcQuantByRow
  )
  expect_error(
    NormalyzerDE:::resolveLimpaQuantByRowFn(new.env(parent = emptyenv())),
    class = "normalyzerde_error"
  )
})

test_that("calculateContrasts supports type='limpa'", {
  test_data <- matrix(
    c(
      10,
      11,
      NA,
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
      7,
      NA,
      NA,
      NA,
      NA,
      NA,
      NA
    ),
    nrow = 5,
    byrow = TRUE
  )
  colnames(test_data) <- paste0("s", seq_len(ncol(test_data)))

  se <- nd_make_summarized_experiment(
    assay = test_data,
    groups = c("A", "A", "A", "B", "B", "B")
  )

  nst <- NormalyzerStatistics(se, logTrans = FALSE)

  if (!requireNamespace("limpa", quietly = TRUE)) {
    expect_error(
      calculateContrasts(
        nst,
        comparisons = "A-B",
        condCol = "group",
        type = "limpa"
      ),
      class = "normalyzerde_error"
    )
    return()
  }

  out <- calculateContrasts(
    nst,
    comparisons = "A-B",
    condCol = "group",
    type = "limpa",
    leastRepCount = 1,
    limpaOptions = limpaOptions(
      quantArgs = list(dpc.slope = 0.7, chunk = 10L),
      deArgs = list(prior.n = 5)
    )
  )

  pvals <- pairwiseCompsP(out)[["A-B"]]
  expect_length(pvals, nrow(test_data))
  expect_true(!all(is.na(pvals)))

  fdrs <- pairwiseCompsFdr(out)[["A-B"]]
  expect_length(fdrs, nrow(test_data))
})

test_that("limpa backend warns when input does not look log2-transformed", {
  testthat::skip_if_not_installed("limpa")

  test_data <- matrix(
    c(
      100,
      110,
      NA,
      130,
      120,
      110,
      NA,
      NA,
      NA,
      90,
      90,
      100,
      50,
      50,
      50,
      NA,
      NA,
      NA,
      70,
      80,
      70,
      70,
      70,
      70,
      NA,
      NA,
      NA,
      NA,
      NA,
      NA
    ),
    nrow = 5,
    byrow = TRUE
  )
  colnames(test_data) <- paste0("s", seq_len(ncol(test_data)))

  se <- nd_make_summarized_experiment(
    assay = test_data,
    groups = c("A", "A", "A", "B", "B", "B")
  )
  nst <- NormalyzerStatistics(se, logTrans = FALSE)

  expect_warning(
    expect_error(
      calculateContrasts(
        nst,
        comparisons = "A-C",
        condCol = "group",
        type = "limpa",
        limpaOptions = limpaOptions(proteinIdCol = NULL)
      ),
      class = "normalyzerde_error"
    ),
    class = "normalyzerde_warning"
  )
})

test_that("limpa backend warns when input looks protein-level", {
  testthat::skip_if_not_installed("limpa")

  test_data <- matrix(
    c(
      10,
      11,
      NA,
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

  se <- nd_make_summarized_experiment(
    assay = test_data,
    groups = c("A", "A", "A", "B", "B", "B"),
    row_data = data.frame(
      `Protein.Group` = paste0("P", seq_len(nrow(test_data)))
    )
  )

  out <- testthat::expect_warning(
    calculateContrasts(
      NormalyzerStatistics(se, logTrans = FALSE),
      comparisons = "A-B",
      condCol = "group",
      type = "limpa",
      leastRepCount = 1,
      limpaOptions = limpaOptions(quantArgs = list(chunk = 10L))
    ),
    class = "normalyzerde_warning"
  )

  expect_equal(nrow(dataMat(out)), nrow(test_data))
  expect_equal(ncol(dataMat(out)), ncol(test_data))
})

test_that("limpa backend supports one-vs-rest contrasts", {
  if (!requireNamespace("limpa", quietly = TRUE)) {
    skip("limpa not installed")
  }

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

  se <- nd_make_summarized_experiment(
    assay = test_data,
    groups = c("A", "A", "B", "B", "C", "C")
  )

  nst <- NormalyzerStatistics(se, logTrans = FALSE)
  out <- calculateContrasts(
    nst,
    condCol = "group",
    type = "limpa",
    oneVsRest = TRUE
  )

  expect_true(all(
    c("A-rest", "B-rest", "C-rest") %in% names(pairwiseCompsP(out))
  ))
})

test_that("limpa backend can summarize peptides to proteins via dpcQuant", {
  if (!requireNamespace("limpa", quietly = TRUE)) {
    skip("limpa not installed")
  }

  set.seed(1)
  n_proteins <- 6
  peptides_per_protein <- 2
  n_samples <- 6

  protein_ids <- paste0("P", seq_len(n_proteins))
  peptide_protein <- rep(protein_ids, each = peptides_per_protein)

  group <- c("A", "A", "A", "B", "B", "B")
  protein_expr <- matrix(
    rnorm(n_proteins * n_samples, mean = 10, sd = 1),
    nrow = n_proteins
  )
  protein_expr[, group == "B"] <- protein_expr[, group == "B"] +
    seq_len(n_proteins) / n_proteins

  test_data <- matrix(
    NA_real_,
    nrow = length(peptide_protein),
    ncol = n_samples
  )
  for (ii in seq_along(peptide_protein)) {
    protein_index <- match(peptide_protein[ii], protein_ids)
    test_data[ii, ] <- protein_expr[protein_index, ] +
      rnorm(n_samples, sd = 0.1)
  }
  missing_mask <- matrix(runif(length(test_data)) < 0.1, nrow = nrow(test_data))
  test_data[missing_mask] <- NA_real_

  colnames(test_data) <- paste0("s", seq_len(n_samples))

  row_anno <- data.frame(
    `Protein.Group` = peptide_protein,
    `Protein.Names` = rep(
      paste0("Prot_", protein_ids),
      each = peptides_per_protein
    ),
    Proteotypic = rep(c(0, 1), times = n_proteins),
    `Precursor.Charge` = rep(c(2, 3), times = n_proteins),
    check.names = FALSE
  )

  se <- nd_make_summarized_experiment(
    assay = test_data,
    groups = group,
    row_data = row_anno
  )

  nst <- NormalyzerStatistics(se, logTrans = FALSE)
  out <- calculateContrasts(
    nst,
    comparisons = "A-B",
    condCol = "group",
    type = "limpa",
    limpaOptions = limpaOptions(
      proteinIdCol = "Protein.Group",
      quantArgs = list(dpc.slope = 0.7, chunk = 10L)
    )
  )

  expect_equal(nrow(dataMat(out)), n_proteins)
  expect_true("Protein.Group" %in% colnames(annotMat(out)))
  expect_true("Protein.Names" %in% colnames(annotMat(out)))
  expect_false("Proteotypic" %in% colnames(annotMat(out)))
  expect_false("Precursor.Charge" %in% colnames(annotMat(out)))
  expect_length(pairwiseCompsP(out)[["A-B"]], n_proteins)
})

test_that("limpaOptions(byRow = TRUE) keeps duplicate protein IDs as separate rows", {
  testthat::skip_if_not_installed("limpa")

  set.seed(1)
  n_proteins <- 6
  peptides_per_protein <- 2
  n_samples <- 6

  protein_ids <- paste0("P", seq_len(n_proteins))
  peptide_protein <- rep(protein_ids, each = peptides_per_protein)

  group <- c("A", "A", "A", "B", "B", "B")
  protein_expr <- matrix(
    rnorm(n_proteins * n_samples, mean = 10, sd = 1),
    nrow = n_proteins
  )
  protein_expr[, group == "B"] <- protein_expr[, group == "B"] +
    seq_len(n_proteins) / n_proteins

  test_data <- matrix(
    NA_real_,
    nrow = length(peptide_protein),
    ncol = n_samples
  )
  for (ii in seq_along(peptide_protein)) {
    protein_index <- match(peptide_protein[ii], protein_ids)
    test_data[ii, ] <- protein_expr[protein_index, ] +
      rnorm(n_samples, sd = 0.1)
  }
  missing_mask <- matrix(runif(length(test_data)) < 0.1, nrow = nrow(test_data))
  test_data[missing_mask] <- NA_real_

  colnames(test_data) <- paste0("s", seq_len(n_samples))

  row_anno <- data.frame(`Protein.Group` = peptide_protein)

  se <- nd_make_summarized_experiment(
    assay = test_data,
    groups = group,
    row_data = row_anno
  )

  nst <- NormalyzerStatistics(se, logTrans = FALSE)
  out <- calculateContrasts(
    nst,
    comparisons = "A-B",
    condCol = "group",
    type = "limpa",
    limpaOptions = limpaOptions(
      byRow = TRUE,
      keep = "elist",
      quantArgs = list(dpc.slope = 0.7, chunk = 10L)
    )
  )

  backend <- backendData(out)[["limpa"]]
  expect_type(backend, "list")
  expect_false(isTRUE(backend$usedDpcQuant))

  expect_equal(nrow(dataMat(out)), nrow(test_data))
})

test_that("limpa backend can auto-estimate DPC via limpaOptions(dpcMethod = ...)", {
  testthat::skip_if_not_installed("limpa")

  test_data <- matrix(
    c(
      10,
      11,
      NA,
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
      7,
      NA,
      NA,
      NA,
      NA,
      NA,
      NA
    ),
    nrow = 5,
    byrow = TRUE
  )
  colnames(test_data) <- paste0("s", seq_len(ncol(test_data)))

  design <- data.frame(
    sample = colnames(test_data),
    group = c("A", "A", "A", "B", "B", "B")
  )
  rownames(design) <- design$sample

  se <- SummarizedExperiment::SummarizedExperiment(
    assay = test_data,
    colData = design,
    rowData = data.frame(feature = paste0("f", seq_len(nrow(test_data))))
  )

  nst <- NormalyzerStatistics(se, logTrans = FALSE)
  out <- calculateContrasts(
    nst,
    comparisons = "A-B",
    condCol = "group",
    type = "limpa",
    leastRepCount = 1,
    limpaOptions = limpaOptions(
      quantArgs = list(chunk = 10L),
      keep = "fit",
      dpcMethod = "dpc"
    )
  )

  backend <- backendData(out)[["limpa"]]
  expect_type(backend, "list")
  expect_equal(backend$dpcMethod, "dpc")
  expect_true(is.numeric(backend$dpc))
  expect_length(backend$dpc, 2)

  expect_true(".global" %in% names(backend$fits))
  expect_true(inherits(backend$fits[[".global"]], "MArrayLM"))
  expect_true(is.matrix(backend$designs[[".global"]]))

  if (exists("dpcON", envir = asNamespace("limpa"), inherits = FALSE)) {
    out_on <- calculateContrasts(
      nst,
      comparisons = "A-B",
      condCol = "group",
      type = "limpa",
      leastRepCount = 1,
      limpaOptions = limpaOptions(
        quantArgs = list(chunk = 10L),
        keep = "fit",
        dpcMethod = "dpcON",
        dpcArgs = list(robust = TRUE)
      )
    )

    backend_on <- backendData(out_on)[["limpa"]]
    expect_type(backend_on, "list")
    expect_equal(backend_on$dpcMethod, "dpcON")
    expect_true(is.numeric(backend_on$dpc))
    expect_length(backend_on$dpc, 2)
    expect_true(isTRUE(backend_on$dpcArgs$robust))

    expect_true(".global" %in% names(backend_on$fits))
    expect_true(inherits(backend_on$fits[[".global"]], "MArrayLM"))
    expect_true(is.matrix(backend_on$designs[[".global"]]))
  }
})

test_that("limpa sample weights are accessible without keeping full fits", {
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

  design <- data.frame(
    sample = colnames(test_data),
    group = c("A", "A", "A", "B", "B", "B")
  )
  rownames(design) <- design$sample

  se <- SummarizedExperiment::SummarizedExperiment(
    assay = test_data,
    colData = design,
    rowData = data.frame(feature = paste0("f", seq_len(nrow(test_data))))
  )

  out <- calculateContrasts(
    NormalyzerStatistics(se, logTrans = FALSE),
    comparisons = "A-B",
    condCol = "group",
    type = "limpa",
    leastRepCount = 1,
    limpaOptions = limpaOptions(
      deArgs = list(sample.weights = TRUE),
      quantArgs = list(chunk = 10L)
    )
  )

  backend <- backendData(out)[["limpa"]]
  expect_type(backend, "list")
  expect_true(is.list(backend$sampleWeights))
  expect_true(".global" %in% names(backend$sampleWeights))
  expect_length(backend$fits, 0)

  weights <- getLimpaSampleWeights(out)
  expect_s3_class(weights, "data.frame")
  expect_equal(colnames(weights), c("comparison", "sample", "sampleWeight"))
  expect_equal(weights$comparison, rep(".global", ncol(test_data)))
  expect_equal(weights$sample, colnames(test_data))
  expect_equal(
    weights$sampleWeight,
    as.numeric(backend$sampleWeights[[".global"]])
  )
  expect_equal(names(backend$sampleWeights[[".global"]]), colnames(test_data))
})

test_that("limpa backend validates limpaOptions(dpcArgs = ...) by method", {
  testthat::skip_if_not_installed("limpa")

  test_data <- matrix(
    c(
      10,
      11,
      NA,
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

  design <- data.frame(
    sample = colnames(test_data),
    group = c("A", "A", "A", "B", "B", "B")
  )
  rownames(design) <- design$sample

  se <- SummarizedExperiment::SummarizedExperiment(
    assay = test_data,
    colData = design,
    rowData = data.frame(feature = paste0("f", seq_len(nrow(test_data))))
  )

  expect_error(
    calculateContrasts(
      NormalyzerStatistics(se, logTrans = FALSE),
      comparisons = "A-B",
      condCol = "group",
      type = "limpa",
      limpaOptions = limpaOptions(
        dpcMethod = "dpc",
        dpcArgs = list(verbose = FALSE)
      )
    ),
    class = "normalyzerde_error"
  )

  expect_error(
    calculateContrasts(
      NormalyzerStatistics(se, logTrans = FALSE),
      comparisons = "A-B",
      condCol = "group",
      type = "limpa",
      limpaOptions = limpaOptions(
        dpcMethod = "dpcCN",
        dpcArgs = list(maxit = 10)
      )
    ),
    class = "normalyzerde_error"
  )

  if (exists("dpcON", envir = asNamespace("limpa"), inherits = FALSE)) {
    expect_error(
      calculateContrasts(
        NormalyzerStatistics(se, logTrans = FALSE),
        comparisons = "A-B",
        condCol = "group",
        type = "limpa",
        limpaOptions = limpaOptions(
          dpcMethod = "dpcON",
          dpcArgs = list(maxit = 10)
        )
      ),
      class = "normalyzerde_error"
    )
  }
})

test_that("limpa backend accepts and records limpaOptions(quantArgs = ...)", {
  testthat::skip_if_not_installed("limpa")

  test_data <- matrix(
    c(
      10,
      11,
      NA,
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

  design <- data.frame(
    sample = colnames(test_data),
    group = c("A", "A", "A", "B", "B", "B")
  )
  rownames(design) <- design$sample

  se <- SummarizedExperiment::SummarizedExperiment(
    assay = test_data,
    colData = design,
    rowData = data.frame(feature = paste0("f", seq_len(nrow(test_data))))
  )

  out <- calculateContrasts(
    NormalyzerStatistics(se, logTrans = FALSE),
    comparisons = "A-B",
    condCol = "group",
    type = "limpa",
    leastRepCount = 1,
    limpaOptions = limpaOptions(
      keep = "fit",
      quantArgs = list(
        sd.quantile.for.logFC = 0.8,
        dpc.slope = 0.1,
        chunk = 1L
      )
    )
  )

  backend <- backendData(out)[["limpa"]]
  expect_type(backend, "list")
  expect_type(backend$quantArgs, "list")
  expect_equal(backend$quantArgs$sd.quantile.for.logFC, 0.8)
  expect_equal(backend$quantArgs$dpc.slope, 0.1)
  expect_equal(backend$quantArgs$chunk, 1L)
  expect_false("dpc" %in% names(backend$quantArgs))
  expect_false("protein.id" %in% names(backend$quantArgs))
  expect_identical(backend$quantArgs$verbose, FALSE)
})

test_that("limpa backend can apply quantile normalization after dpcQuant", {
  testthat::skip_if_not_installed("limpa")

  test_data <- matrix(
    c(
      10,
      11,
      NA,
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

  design <- data.frame(
    sample = colnames(test_data),
    group = c("A", "A", "A", "B", "B", "B")
  )
  rownames(design) <- design$sample

  se <- SummarizedExperiment::SummarizedExperiment(
    assay = test_data,
    colData = design,
    rowData = data.frame(feature = paste0("f", seq_len(nrow(test_data))))
  )

  out <- calculateContrasts(
    NormalyzerStatistics(se, logTrans = FALSE),
    comparisons = "A-B",
    condCol = "group",
    type = "limpa",
    leastRepCount = 1,
    limpaOptions = limpaOptions(
      keep = "elist",
      postQuantNorm = "quantile",
      quantArgs = list(chunk = 10L)
    )
  )

  backend <- backendData(out)[["limpa"]]
  expect_equal(backend$postQuantNorm, "quantile")

  y <- backend$elists[[".global"]]
  expect_true(inherits(y, "EList"))

  sortedE <- apply(y$E, 2, sort)
  ref <- sortedE[, 1]
  for (ii in seq_len(ncol(sortedE))) {
    expect_equal(sortedE[, ii], ref, tolerance = 1e-10)
  }
})

test_that("limpa backend keeps fits keyed by one-vs-rest comparisons", {
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

  design <- data.frame(
    sample = colnames(test_data),
    group = c("A", "A", "B", "B", "C", "C")
  )
  rownames(design) <- design$sample

  se <- SummarizedExperiment::SummarizedExperiment(
    assay = test_data,
    colData = design,
    rowData = data.frame(feature = paste0("f", seq_len(nrow(test_data))))
  )

  nst <- NormalyzerStatistics(se, logTrans = FALSE)
  out <- calculateContrasts(
    nst,
    condCol = "group",
    type = "limpa",
    oneVsRest = TRUE,
    limpaOptions = limpaOptions(keep = "fit")
  )

  backend <- backendData(out)[["limpa"]]
  expect_type(backend, "list")
  expect_true(all(c("A-rest", "B-rest", "C-rest") %in% names(backend$fits)))
})

test_that("limpa backend can keep the quantified EList when summarizing peptides", {
  testthat::skip_if_not_installed("limpa")

  set.seed(1)
  n_proteins <- 6
  peptides_per_protein <- 2
  n_samples <- 6

  protein_ids <- paste0("P", seq_len(n_proteins))
  peptide_protein <- rep(protein_ids, each = peptides_per_protein)

  group <- c("A", "A", "A", "B", "B", "B")
  protein_expr <- matrix(
    rnorm(n_proteins * n_samples, mean = 10, sd = 1),
    nrow = n_proteins
  )
  protein_expr[, group == "B"] <- protein_expr[, group == "B"] +
    seq_len(n_proteins) / n_proteins

  test_data <- matrix(
    NA_real_,
    nrow = length(peptide_protein),
    ncol = n_samples
  )
  for (ii in seq_along(peptide_protein)) {
    protein_index <- match(peptide_protein[ii], protein_ids)
    test_data[ii, ] <- protein_expr[protein_index, ] +
      rnorm(n_samples, sd = 0.1)
  }
  missing_mask <- matrix(runif(length(test_data)) < 0.1, nrow = nrow(test_data))
  test_data[missing_mask] <- NA_real_

  colnames(test_data) <- paste0("s", seq_len(n_samples))

  design <- data.frame(sample = colnames(test_data), group = group)
  rownames(design) <- design$sample

  row_anno <- data.frame(`Protein.Group` = peptide_protein)

  se <- SummarizedExperiment::SummarizedExperiment(
    assay = test_data,
    colData = design,
    rowData = row_anno
  )

  nst <- NormalyzerStatistics(se, logTrans = FALSE)
  out <- calculateContrasts(
    nst,
    comparisons = "A-B",
    condCol = "group",
    type = "limpa",
    limpaOptions = limpaOptions(
      proteinIdCol = "Protein.Group",
      keep = "elist",
      quantArgs = list(chunk = 10L)
    )
  )

  backend <- backendData(out)[["limpa"]]
  expect_type(backend, "list")
  expect_true(isTRUE(backend$usedDpcQuant))
  expect_equal(backend$proteinIdCol, "Protein.Group")
  expect_true(inherits(backend$quantifiedEList, "EList"))
  expect_true(inherits(backend$elists[[".global"]], "EList"))
  expect_equal(nrow(dataMat(out)), n_proteins)
})
