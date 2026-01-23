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
    limpaDpcSlope = 0.7,
    limpaChunk = 10L,
    limpaDEArgs = list(prior.n = 5)
  )

  pvals <- pairwiseCompsP(out)[["A-B"]]
  expect_length(pvals, nrow(test_data))
  expect_true(!all(is.na(pvals)))

  fdrs <- pairwiseCompsFdr(out)[["A-B"]]
  expect_length(fdrs, nrow(test_data))
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
    limpaDpcSlope = 0.7,
    limpaChunk = 10L
  )

  expect_equal(nrow(dataMat(out)), n_proteins)
  expect_true("Protein.Group" %in% colnames(annotMat(out)))
  expect_length(pairwiseCompsP(out)[["A-B"]], n_proteins)
})
