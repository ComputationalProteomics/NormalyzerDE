context("limpa quantifiedRds validation")

nd_completed_limpa_matrix <- function() {
  matrix(
    c(
      10, 10, 10, 11, 11, 11,
      9, 9, 9, 10, 10, 10,
      8, 8, 8, 9, 9, 9,
      7, 7, 7, 8, 8, 8
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(paste0("f", seq_len(4)), paste0("s", seq_len(6)))
  )
}

nd_cached_genes_df <- function(E) {
  genes <- data.frame(
    feature = paste0("g", seq_len(nrow(E))),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  row_ids <- rownames(E)
  if (
    !is.null(row_ids) &&
      !anyNA(row_ids) &&
      anyDuplicated(row_ids) == 0 &&
      all(nzchar(row_ids))
  ) {
    rownames(genes) <- row_ids
  }

  genes
}

nd_cached_elist <- function(
  E = nd_completed_limpa_matrix(),
  n_observations = matrix(1L, nrow(E), ncol(E), dimnames = dimnames(E)),
  standard_error = matrix(0.1, nrow(E), ncol(E), dimnames = dimnames(E)),
  genes = nd_cached_genes_df(E)
) {
  structure(
    list(
      E = E,
      other = list(
        n.observations = n_observations,
        standard.error = standard_error
      ),
      genes = genes
    ),
    class = "EList"
  )
}

nd_calculate_with_cached_elist <- function(
  cached,
  assay = nd_completed_limpa_matrix()
) {
  testthat::skip_if_not_installed("limpa")

  rds_path <- withr::local_tempfile(
    pattern = "limpa_quantified_",
    fileext = ".rds"
  )
  saveRDS(cached, file = rds_path)

  se <- nd_make_summarized_experiment(
    assay = assay,
    groups = c("A", "A", "A", "B", "B", "B")
  )

  calculateContrasts(
    NormalyzerStatistics(se, logTrans = FALSE),
    comparisons = "A-B",
    condCol = "group",
    type = "limpa",
    leastRepCount = 1,
    limpaOptions = limpaOptions(quantifiedRds = rds_path, keep = "elist")
  )
}

test_that("calculateContrasts validates quantifiedRds EList structure", {
  base_matrix <- nd_completed_limpa_matrix()

  expect_error(
    nd_calculate_with_cached_elist(list(E = base_matrix)),
    regexp = "must contain a limma .*EList",
    class = "normalyzerde_error"
  )

  expect_error(
    nd_calculate_with_cached_elist(
      structure(
        list(
          E = 1,
          other = list(
            n.observations = matrix(1L, nrow = 1, ncol = 1),
            standard.error = matrix(0.1, nrow = 1, ncol = 1)
          ),
          genes = data.frame(
            feature = "g1",
            stringsAsFactors = FALSE,
            check.names = FALSE
          )
        ),
        class = "EList"
      )
    ),
    regexp = "must contain a matrix element .*E",
    class = "normalyzerde_error"
  )

  no_colnames <- base_matrix
  colnames(no_colnames) <- NULL
  expect_error(
    nd_calculate_with_cached_elist(
      nd_cached_elist(E = no_colnames)
    ),
    regexp = "must contain sample column names",
    class = "normalyzerde_error"
  )

  duplicate_rows <- base_matrix
  rownames(duplicate_rows) <- c("dup", "dup", "ok3", "ok4")
  expect_error(
    nd_calculate_with_cached_elist(
      nd_cached_elist(E = duplicate_rows)
    ),
    regexp = "must have unique, non-empty row names",
    class = "normalyzerde_error"
  )
})

test_that("calculateContrasts validates quantifiedRds observation metadata", {
  base_matrix <- nd_completed_limpa_matrix()

  expect_error(
    nd_calculate_with_cached_elist(
      nd_cached_elist(n_observations = NULL)
    ),
    regexp = "must contain .*other\\$n\\.observations",
    class = "normalyzerde_error"
  )

  expect_error(
    nd_calculate_with_cached_elist(
      nd_cached_elist(
        n_observations = matrix(
          1L,
          nrow = nrow(base_matrix) - 1,
          ncol = ncol(base_matrix),
          dimnames = list(rownames(base_matrix)[-1], colnames(base_matrix))
        )
      )
    ),
    regexp = "other\\$n\\.observations.*same dimensions",
    class = "normalyzerde_error"
  )

  expect_error(
    nd_calculate_with_cached_elist(
      nd_cached_elist(standard_error = NULL)
    ),
    regexp = "must contain .*other\\$standard\\.error",
    class = "normalyzerde_error"
  )

  expect_error(
    nd_calculate_with_cached_elist(
      nd_cached_elist(
        standard_error = matrix(
          0.1,
          nrow = nrow(base_matrix),
          ncol = ncol(base_matrix) - 1,
          dimnames = list(rownames(base_matrix), colnames(base_matrix)[-1])
        )
      )
    ),
    regexp = "other\\$standard\\.error.*same dimensions",
    class = "normalyzerde_error"
  )
})

test_that("calculateContrasts validates quantifiedRds alignment with the data matrix", {
  base_matrix <- nd_completed_limpa_matrix()

  assay_with_na <- base_matrix
  assay_with_na[1, 1] <- NA_real_
  expect_error(
    nd_calculate_with_cached_elist(
      nd_cached_elist(),
      assay = assay_with_na
    ),
    regexp = "completed expression matrix without NA",
    class = "normalyzerde_error"
  )

  missing_sample <- base_matrix[, -1, drop = FALSE]
  expect_error(
    nd_calculate_with_cached_elist(
      nd_cached_elist(E = missing_sample)
    ),
    regexp = "missing sample columns required by the data matrix",
    class = "normalyzerde_error"
  )

  bad_genes <- data.frame(
    feature = "g1",
    row.names = "unmatched",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  expect_error(
    nd_calculate_with_cached_elist(
      nd_cached_elist(genes = bad_genes)
    ),
    regexp = "genes rows could not be aligned",
    class = "normalyzerde_error"
  )
})

test_that("calculateContrasts errors when quantifiedRds file is missing", {
  testthat::skip_if_not_installed("limpa")

  missing_rds <- file.path(
    withr::local_tempdir(pattern = "missing_quantified_rds_"),
    "missing_quantified.rds"
  )

  se <- nd_make_summarized_experiment(
    assay = nd_completed_limpa_matrix(),
    groups = c("A", "A", "A", "B", "B", "B")
  )

  expect_error(
    calculateContrasts(
      NormalyzerStatistics(se, logTrans = FALSE),
      comparisons = "A-B",
      condCol = "group",
      type = "limpa",
      leastRepCount = 1,
      limpaOptions = limpaOptions(quantifiedRds = missing_rds)
    ),
    regexp = "file does not exist",
    class = "normalyzerde_error"
  )
})
