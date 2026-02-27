context("limpa backend")

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

  if (!requireNamespace("limpa", quietly = TRUE)) {
    expect_error(
      calculateContrasts(
        nst,
        comparisons = "A-B",
        condCol = "group",
        type = "limpa"
      ),
      "requires.*limpa"
    )
    return()
  }

  out <- calculateContrasts(
    nst,
    comparisons = "A-B",
    condCol = "group",
    type = "limpa",
    leastRepCount = 1,
    limpaQuantArgs = list(dpc.slope = 0.7, chunk = 10L),
    limpaDEArgs = list(prior.n = 5)
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

  expect_warning(
    expect_error(
      calculateContrasts(
        nst,
        comparisons = "A-C",
        condCol = "group",
        type = "limpa",
        limpaProteinIdCol = NULL
      ),
      "issues in your contrast"
    ),
    "log2"
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

  design <- data.frame(
    sample = colnames(test_data),
    group = c("A", "A", "A", "B", "B", "B")
  )
  rownames(design) <- design$sample

  se <- SummarizedExperiment::SummarizedExperiment(
    assay = test_data,
    colData = design,
    rowData = data.frame(`Protein.Group` = paste0("P", seq_len(nrow(test_data))))
  )

  out <- testthat::expect_warning(
    calculateContrasts(
      NormalyzerStatistics(se, logTrans = FALSE),
      comparisons = "A-B",
      condCol = "group",
      type = "limpa",
      leastRepCount = 1,
      limpaQuantArgs = list(chunk = 10L)
    ),
    "no peptide/precursor-to-protein summarization"
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

  design <- data.frame(sample = colnames(test_data), group = group)
  rownames(design) <- design$sample

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
    limpaProteinIdCol = "Protein.Group",
    limpaQuantArgs = list(dpc.slope = 0.7, chunk = 10L)
  )

  expect_equal(nrow(dataMat(out)), n_proteins)
  expect_true("Protein.Group" %in% colnames(annotMat(out)))
  expect_true("Protein.Names" %in% colnames(annotMat(out)))
  expect_false("Proteotypic" %in% colnames(annotMat(out)))
  expect_false("Precursor.Charge" %in% colnames(annotMat(out)))
  expect_length(pairwiseCompsP(out)[["A-B"]], n_proteins)
})

test_that("limpaByRow keeps duplicate protein IDs as separate rows", {
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
    limpaByRow = TRUE,
    limpaKeep = "elist",
    limpaQuantArgs = list(dpc.slope = 0.7, chunk = 10L)
  )

  backend <- backendData(out)[["limpa"]]
  expect_type(backend, "list")
  expect_false(isTRUE(backend$usedDpcQuant))

  expect_equal(nrow(dataMat(out)), nrow(test_data))
})

test_that("limpa backend can auto-estimate DPC via limpaDpcMethod", {
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
    limpaQuantArgs = list(chunk = 10L),
    limpaKeep = "fit",
    limpaDpcMethod = "dpc"
  )

  backend <- backendData(out)[["limpa"]]
  expect_type(backend, "list")
  expect_equal(backend$dpcMethod, "dpc")
  expect_true(is.numeric(backend$dpc))
  expect_length(backend$dpc, 2)

  expect_true(".global" %in% names(backend$fits))
  expect_true(inherits(backend$fits[[".global"]], "MArrayLM"))
  expect_true(is.matrix(backend$designs[[".global"]]))
})

test_that("limpa backend validates limpaDpcArgs by method", {
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
      limpaDpcMethod = "dpc",
      limpaDpcArgs = list(verbose = FALSE)
    ),
    "limpaDpcArgs"
  )

  expect_error(
    calculateContrasts(
      NormalyzerStatistics(se, logTrans = FALSE),
      comparisons = "A-B",
      condCol = "group",
      type = "limpa",
      limpaDpcMethod = "dpcCN",
      limpaDpcArgs = list(maxit = 10)
    ),
    "limpaDpcArgs"
  )
})

test_that("limpa backend accepts and records limpaQuantArgs", {
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
    limpaKeep = "fit",
    limpaQuantArgs = list(
      sd.quantile.for.logFC = 0.8,
      dpc.slope = 0.1,
      chunk = 1L
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
    limpaKeep = "fit"
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
    limpaProteinIdCol = "Protein.Group",
    limpaKeep = "elist",
    limpaQuantArgs = list(chunk = 10L)
  )

  backend <- backendData(out)[["limpa"]]
  expect_type(backend, "list")
  expect_true(isTRUE(backend$usedDpcQuant))
  expect_equal(backend$proteinIdCol, "Protein.Group")
  expect_true(inherits(backend$quantifiedEList, "EList"))
  expect_true(inherits(backend$elists[[".global"]], "EList"))
  expect_equal(nrow(dataMat(out)), n_proteins)
})
