context("DIANN parquet design alignment")

test_that("DIANN parquet setup respects design subsets and order", {
  testthat::skip_if_not_installed("arrow")

  tmp_dir <- withr::local_tempdir(pattern = "diann_subset_")
  data_path <- nd_write_diann_parquet_fixture(tmp_dir)

  subset_design <- nd_balanced_diann_design()[c(6, 3, 1), , drop = FALSE]
  se <- nd_load_diann_contrast_object(data_path, subset_design)

  expect_identical(
    colnames(SummarizedExperiment::assay(se)),
    subset_design$sample
  )
  expect_identical(
    as.character(SummarizedExperiment::colData(se)$sample),
    subset_design$sample
  )
  expect_false("S2" %in% colnames(SummarizedExperiment::assay(se)))
})

test_that("DIANN parquet retains samples whose rows are fully removed by q filtering", {
  testthat::skip_if_not_installed("arrow")

  tmp_dir <- withr::local_tempdir(pattern = "diann_qdrop_sample_")
  data_path <- file.path(tmp_dir, "diann_report_qdrop.parquet")
  report <- nd_make_diann_precursor_report()
  report$`Q.Value`[report$Run == "S2"] <- 0.02
  nd_write_parquet(report, data_path)

  se <- nd_load_diann_contrast_object(data_path, nd_balanced_diann_design())
  assay <- SummarizedExperiment::assay(se)

  expect_identical(colnames(assay), nd_balanced_diann_design()$sample)
  expect_true(all(is.na(assay[, "S2"])))
  expect_false(all(is.na(assay[, "S1"])))
})

test_that("limma route is invariant to DIANN parquet design row order", {
  testthat::skip_if_not_installed("arrow")

  tmp_dir <- withr::local_tempdir(pattern = "diann_limma_order_")
  data_path <- nd_write_diann_parquet_fixture(tmp_dir)

  design <- nd_balanced_diann_design()
  reversed_design <- design[rev(seq_len(nrow(design))), , drop = FALSE]

  out_default <- nd_run_diann_contrast(data_path, design, type = "limma")
  out_reversed <- nd_run_diann_contrast(
    data_path,
    reversed_design,
    type = "limma"
  )

  expect_identical(
    colnames(dataMat(out_default)),
    design$sample
  )
  expect_identical(
    colnames(dataMat(out_reversed)),
    reversed_design$sample
  )

  expect_equal(
    nd_extract_contrast_table(out_default, "A-B", "Precursor.Id"),
    nd_extract_contrast_table(out_reversed, "A-B", "Precursor.Id")
  )
})

test_that("limpa route is invariant to DIANN parquet design row order", {
  testthat::skip_if_not_installed("arrow")
  testthat::skip_if_not_installed("limpa")

  tmp_dir <- withr::local_tempdir(pattern = "diann_limpa_order_")
  data_path <- nd_write_diann_parquet_fixture(tmp_dir)

  design <- nd_balanced_diann_design()
  reversed_design <- design[rev(seq_len(nrow(design))), , drop = FALSE]

  out_default <- nd_run_diann_contrast(data_path, design, type = "limpa")
  out_reversed <- nd_run_diann_contrast(
    data_path,
    reversed_design,
    type = "limpa"
  )

  expect_identical(
    colnames(dataMat(out_default)),
    design$sample
  )
  expect_identical(
    colnames(dataMat(out_reversed)),
    reversed_design$sample
  )

  expect_equal(
    nd_extract_contrast_table(out_default, "A-B", "Protein.Group"),
    nd_extract_contrast_table(out_reversed, "A-B", "Protein.Group")
  )
  expect_equal(
    nd_extract_limpa_elist_matrix(out_default),
    nd_extract_limpa_elist_matrix(out_reversed)
  )
})
