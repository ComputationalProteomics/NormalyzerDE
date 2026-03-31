context("calculateStatistics.R")

data("example_data_only_values")


# test_that("reduceTechnicalReplicates", {
#
#     tech_rep <- c("a", "a", "b", "b", "c", "c", "d", "d")
#     test_data <- data.frame(
#         c(1,1,1),
#         c(1,2,1),
#         c(3,3,3),
#         c(5,3,3),
#         c(5,5,4),
#         c(5,5,5),
#         c(7,7,7),
#         c(7,9,7))
#     colnames(test_data) <- c("a1", "a2", "b1", "b2", "c1", "c2", "d1", "d2")
#
#     expect_out_data <- as.matrix(data.frame(
#         "a"=c(1,1.5,1),
#         "b"=c(4,3,3),
#         "c"=c(5,5,4.5),
#         "d"=c(7,8,7)))
#
#     out <- reduceTechnicalReplicates(test_data, tech_rep)
#
#     expect_that(
#         all.equal(
#             expect_out_data,
#             out
#         ),
#         is_true()
#     )
# })
#
# test_that("reduceDesignTechRep", {
#
#     test_df <- data.frame(
#         sample=c("a1", "a2", "a3", "b1", "b2", "c1", "c2", "d1"),
#         group=c(rep("A", 5), rep("B", 3)),
#         techrep=c("a", "a", "a", "b", "b", "c", "c", "d")
#     )
#
#     expected_out_df <- data.frame(
#         sample=c("a1", "b1", "c1", "d1"),
#         group=c(rep("A", 2), rep("B", 2)),
#         techrep=c("a", "b", "c", "d")
#     )
#
#     out <- reduceDesignTechRep(test_df, test_df$techrep)
#
#     expect_that(
#         all.equal(
#             expected_out_df,
#             out
#         ),
#         is_true()
#     )
# })

test_that("reduceTechnicalReplicates", {
  test_data <- data.frame(
    c(1, 1, 1),
    c(1, 1.5, NA),
    c(4, 2, 3),
    c(3, 2.5, 3),
    c(5, 3.5, 3),
    c(4, 3, 4),
    c(6, 7, 5),
    c(7, 9, NA)
  )
  colnames(test_data) <- c("a1", "a2", "a3", "b1", "b2", "c1", "c2", "d1")

  test_df <- data.frame(
    sample = c("a1", "a2", "a3", "b1", "b2", "c1", "c2", "d1"),
    group = c(rep("A", 5), rep("B", 3)),
    techrep = c("a", "a", "a", "b", "b", "c", "c", "d")
  )

  se <- SummarizedExperiment::SummarizedExperiment(
    assay = as.matrix(test_data),
    colData = test_df,
    rowData = data.frame(annot = paste0("Pep", seq_len(3)))
  )

  expect_out_data <- as.matrix(data.frame(
    "a1.a2.a3" = c(2, 1.5, 2),
    "b1.b2" = c(4, 3, 3),
    "c1.c2" = c(5, 5, 4.5),
    "d1" = c(7, 9, NA)
  ))

  expect_out_design <- data.frame(
    sample = c("a1.a2.a3", "b1.b2", "c1.c2", "d1"),
    group = c(rep("A", 2), rep("B", 2)),
    techrep = c("a", "b", "c", "d")
  )
  expect_out_design$sample <- as.character(expect_out_design$sample)
  rownames(expect_out_design) <- expect_out_design$sample

  expect_out_annot <- data.frame(
    annot = paste0("Pep", seq_len(3))
  )

  out_se <- reduceTechnicalReplicates(se, "techrep", "sample")

  expect_true(
    all.equal(
      expect_out_data,
      SummarizedExperiment::assay(out_se)
    ),
    "Reduced data check"
  )

  expect_true(
    all.equal(
      expect_out_design,
      data.frame(SummarizedExperiment::colData(out_se))
    ),
    "Reduced design check"
  )

  expect_true(
    all.equal(
      expect_out_annot,
      data.frame(SummarizedExperiment::rowData(out_se))
    ),
    "Annotation remains identical"
  )
})


# Statistics
test_that("calculateWelch_data_test", {
  small_df <- log2(head(example_data_only_values))

  # small_df content
  #
  # s_500amol_1 s_500amol_2 s_500amol_3 s_2500amol_1 s_2500amol_2 s_2500amol_3
  # [1,]    26.54878    26.52997    26.63352     26.89629     26.79226     26.74493
  # [2,]    26.62507    26.39669    26.75476     26.56381     26.36598     26.17868
  # [3,]    28.42273    28.45091    28.58193           NA     28.56563     28.66707
  # [4,]    25.53961    25.90364    25.38250     23.43634     24.55913     24.19715
  # [5,]    27.36890    27.29489    27.37707     26.60043     26.34719     26.22433
  # [6,]    28.90688    28.83352    28.85170     28.12027     28.31122     28.08823

  expected_p <- c(0.01484, 0.21877, 0.17375, 0.0271, 0.01001, 0.00592)
  expected_fdr <- c(0.02968, 0.21877, 0.20851, 0.04065, 0.02968, 0.02968)
  expected_fold <- c(-0.2404, 0.22268, -0.13116, 1.54437, 0.95631, 0.69079)
  expected_ave <- c(26.70042, 26.5699, 28.52404, 24.94016, 26.82584, 28.36693)

  header <- c(
    rep(1, 3),
    rep(2, 3),
    rep(3, 3),
    rep(4, 3),
    rep(5, 3),
    rep(6, 3),
    rep(7, 3),
    rep(8, 3),
    rep(9, 3)
  )

  out <- calculateWelch(small_df, header, c(4, 5))

  expect_true(all.equal(expected_p, round(out[["P"]], 5)))
  expect_true(all.equal(expected_fdr, round(out[["FDR"]], 5)))
  expect_true(all.equal(expected_fold, round(out[["Fold"]], 5)))
  expect_true(all.equal(expected_ave, round(out[["Ave"]], 5)))
})

test_that("calculateWelch_limited", {
  small_df <- data.frame(
    "a1" = c(1, 1, 1),
    "a2" = c(1.1, 1, 2),
    "a3" = c(0.9, 1, 3),
    "b1" = c(2, 1, 4),
    "b2" = c(2.1, 1, 5),
    "b3" = c(1.9, 1, 6)
  )

  expected_p <- c(0.00026, NA, 0.02131)
  expected_fdr <- c(0.00051, NA, 0.02131)
  expected_fold <- c(-1, 0, -3)
  expected_ave <- c(1.5, 1.0, 3.5)

  out <- calculateWelch(small_df, c(1, 1, 1, 2, 2, 2), c(1, 2))

  expect_true(all.equal(expected_p, round(out[["P"]], 5)))
  expect_true(all.equal(expected_fdr, round(out[["FDR"]], 5)))
  expect_true(all.equal(expected_fold, round(out[["Fold"]], 5)))
  expect_true(all.equal(expected_ave, round(out[["Ave"]], 5)))
})

test_that("calculateLimmaContrast_data_test", {
  small_df <- log2(head(example_data_only_values))[, seq(10, 15)]

  # small_df content
  #
  # s_500amol_1 s_500amol_2 s_500amol_3 s_2500amol_1 s_2500amol_2 s_2500amol_3
  # [1,]    26.54878    26.52997    26.63352     26.89629     26.79226     26.74493
  # [2,]    26.62507    26.39669    26.75476     26.56381     26.36598     26.17868
  # [3,]    28.42273    28.45091    28.58193           NA     28.56563     28.66707
  # [4,]    25.53961    25.90364    25.38250     23.43634     24.55913     24.19715
  # [5,]    27.36890    27.29489    27.37707     26.60043     26.34719     26.22433
  # [6,]    28.90688    28.83352    28.85170     28.12027     28.31122     28.08823

  header <- c(rep(4, 3), rep(5, 3))
  levels <- c(4, 5)

  Variable <- as.factor(header)
  model <- ~ 0 + Variable
  limmaDesign <- stats::model.matrix(model)
  limmaFit <- limma::lmFit(small_df, limmaDesign)

  out <- calculateLimmaContrast(
    small_df,
    limmaDesign,
    limmaFit,
    levels,
    useIntensityTrend = FALSE
  )

  expected_p <- c(0.0169, 0.14704, 0.2116, 0.00154, 8e-05, 0.00012)
  expected_fdr <- c(0.02535, 0.17645, 0.2116, 0.00308, 0.00037, 0.00037)
  expected_ave <- c(26.69096, 26.48083, 28.53765, 24.8364, 26.8688, 28.51864)
  expected_fold <- c(-0.2404, 0.22268, -0.13116, 1.54437, 0.95631, 0.69079)

  # These P-values vary slightly with versions
  # Commented out not to break tests with version updates
  #expect_true(all.equal(expected_p, round(out[["P"]], 5)))
  #expect_true(all.equal(expected_fdr, round(out[["FDR"]], 5)))
  expect_true(all.equal(expected_fold, round(out[["Fold"]], 5)))
  expect_true(all.equal(expected_ave, round(out[["Ave"]], 5)))
})

test_that("calculateLimmaContrast_data_batch_test", {
  small_df <- log2(head(example_data_only_values))[, seq(10, 15)]

  # small_df content
  #
  # s_500amol_1 s_500amol_2 s_500amol_3 s_2500amol_1 s_2500amol_2 s_2500amol_3
  # [1,]    26.54878    26.52997    26.63352     26.89629     26.79226     26.74493
  # [2,]    26.62507    26.39669    26.75476     26.56381     26.36598     26.17868
  # [3,]    28.42273    28.45091    28.58193           NA     28.56563     28.66707
  # [4,]    25.53961    25.90364    25.38250     23.43634     24.55913     24.19715
  # [5,]    27.36890    27.29489    27.37707     26.60043     26.34719     26.22433
  # [6,]    28.90688    28.83352    28.85170     28.12027     28.31122     28.08823

  header <- c(rep(4, 3), rep(5, 3))
  batch <- c(rep(c(1, 2), 3))
  levels <- c(4, 5)

  Variable <- as.factor(header)
  Batch <- as.factor(batch)
  model <- ~ 0 + Variable + Batch
  limmaDesign <- stats::model.matrix(model)
  limmaFit <- limma::lmFit(small_df, limmaDesign)

  out <- calculateLimmaContrast(
    small_df,
    limmaDesign,
    limmaFit,
    levels,
    useIntensityTrend = FALSE
  )

  expected_p <- c(0.04056, 0.28436, 0.28509, 0.00571, 0.00053, 0.00047)
  expected_fdr <- c(0.06083, 0.28509, 0.28509, 0.01142, 0.00158, 0.00158)
  expected_ave <- c(26.69096, 26.48083, 28.53765, 24.8364, 26.8688, 28.51864)
  expected_fold <- c(-0.24587, 0.17468, -0.12881, 1.49441, 0.95415, 0.64867)

  # These P-values vary slightly with versions
  # Commented out not to break tests with version updates
  #expect_true(all.equal(expected_p, round(out[["P"]], 5)))
  #expect_true(all.equal(expected_fdr, round(out[["FDR"]], 5)))
  expect_true(all.equal(expected_fold, round(out[["Fold"]], 5)))
  expect_true(all.equal(expected_ave, round(out[["Ave"]], 5)))
})

test_that("calculateContrasts_impute", {
  test_data <- matrix(
    c(
      10,
      12,
      11,
      9,
      NA,
      NA,
      NA,
      NA,
      5,
      NA,
      4,
      NA,
      NA,
      NA,
      NA,
      NA,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0
    ),
    nrow = 3,
    byrow = TRUE
  )
  colnames(test_data) <- c("A1", "A2", "A3", "A4", "B1", "B2", "B3", "B4")

  test_df <- data.frame(
    sample = colnames(test_data),
    group = c(rep("A", 4), rep("B", 4))
  )
  rownames(test_df) <- test_df$sample

  se <- SummarizedExperiment::SummarizedExperiment(
    assay = test_data,
    colData = test_df,
    rowData = data.frame(annot = paste0("Pep", seq_len(3)))
  )

  nst_no_impute <- NormalyzerStatistics(se)
  out_no_impute <- suppressWarnings(
    calculateContrasts(
      nst_no_impute,
      comparisons = c("A-B"),
      condCol = "group",
      leastRepCount = 0,
      impute = FALSE
    )
  )
  fold_no_impute <- pairwiseCompsFold(out_no_impute)[["A-B"]]

  expect_true(is.na(fold_no_impute[1]))
  expect_true(is.na(fold_no_impute[2]))
  expect_true(all.equal(0, round(fold_no_impute[3], 5)))

  nst_impute1 <- NormalyzerStatistics(se)
  out_impute1 <- suppressWarnings(
    calculateContrasts(
      nst_impute1,
      comparisons = c("A-B"),
      condCol = "group",
      leastRepCount = 0,
      impute = TRUE,
      imputeMinFraction = 1
    )
  )
  fold_impute1 <- pairwiseCompsFold(out_impute1)[["A-B"]]

  expect_true(all.equal(10.5, round(fold_impute1[1], 5)))
  expect_true(is.na(fold_impute1[2]))
  expect_true(all.equal(0, round(fold_impute1[3], 5)))

  nst_impute05 <- NormalyzerStatistics(se)
  out_impute05 <- suppressWarnings(
    calculateContrasts(
      nst_impute05,
      comparisons = c("A-B"),
      condCol = "group",
      leastRepCount = 0,
      impute = TRUE,
      imputeMinFraction = 0.5
    )
  )
  fold_impute05 <- pairwiseCompsFold(out_impute05)[["A-B"]]

  expect_true(all.equal(10.5, round(fold_impute05[1], 5)))
  expect_true(all.equal(4.5, round(fold_impute05[2], 5)))
  expect_true(all.equal(0, round(fold_impute05[3], 5)))
})

test_that("calculateContrasts_subsetByComparison_impute_scope", {
  test_data <- matrix(
    c(
      10,
      11,
      NA,
      NA,
      NA,
      NA
    ),
    nrow = 1,
    byrow = TRUE
  )
  colnames(test_data) <- c("A1", "A2", "B1", "B2", "C1", "C2")

  test_df <- data.frame(
    sample = colnames(test_data),
    group = c(rep("A", 2), rep("B", 2), rep("C", 2))
  )
  rownames(test_df) <- test_df$sample

  se <- SummarizedExperiment::SummarizedExperiment(
    assay = test_data,
    colData = test_df,
    rowData = data.frame(annot = "Pep1")
  )

  nst_global <- NormalyzerStatistics(se)
  out_global <- calculateContrasts(
    nst_global,
    comparisons = c("A-B", "B-C"),
    condCol = "group",
    type = "welch",
    leastRepCount = 0,
    impute = TRUE,
    imputeMinFraction = 1,
    subsetByComparison = FALSE
  )
  globalFolds <- pairwiseCompsFold(out_global)

  expect_true(all.equal(0.5, round(globalFolds[["A-B"]][1], 5)))
  expect_true(all.equal(0, round(globalFolds[["B-C"]][1], 5)))

  nst_subset <- NormalyzerStatistics(se)
  out_subset <- calculateContrasts(
    nst_subset,
    comparisons = c("A-B", "B-C"),
    condCol = "group",
    type = "welch",
    leastRepCount = 0,
    impute = TRUE,
    imputeMinFraction = 1,
    subsetByComparison = TRUE
  )
  subsetFolds <- pairwiseCompsFold(out_subset)

  expect_true(all.equal(0.5, round(subsetFolds[["A-B"]][1], 5)))
  expect_true(is.na(subsetFolds[["B-C"]][1]))
})

test_that("generateAnnotatedMatrix keeps comparison-specific averages", {
  test_data <- matrix(
    c(
      10,
      11,
      NA,
      NA,
      NA,
      NA
    ),
    nrow = 1,
    byrow = TRUE
  )
  colnames(test_data) <- c("A1", "A2", "B1", "B2", "C1", "C2")

  test_df <- data.frame(
    sample = colnames(test_data),
    group = c(rep("A", 2), rep("B", 2), rep("C", 2))
  )
  rownames(test_df) <- test_df$sample

  se <- SummarizedExperiment::SummarizedExperiment(
    assay = test_data,
    colData = test_df,
    rowData = data.frame(annot = "Pep1")
  )

  nst_multi <- NormalyzerStatistics(se)
  out_multi <- calculateContrasts(
    nst_multi,
    comparisons = c("A-B", "B-C"),
    condCol = "group",
    type = "welch",
    leastRepCount = 0,
    impute = TRUE,
    imputeMinFraction = 1,
    subsetByComparison = TRUE
  )
  annot_multi <- generateAnnotatedMatrix(out_multi)

  expect_false("featureAvg" %in% colnames(annot_multi))
  expect_true(all(
    c("A-B_featureAvg", "B-C_featureAvg") %in% colnames(annot_multi)
  ))
  expect_equal(
    annot_multi[["A-B_featureAvg"]],
    pairwiseCompsAve(out_multi)[["A-B"]]
  )
  expect_equal(
    annot_multi[["B-C_featureAvg"]],
    pairwiseCompsAve(out_multi)[["B-C"]]
  )

  nst_single <- NormalyzerStatistics(se)
  out_single <- calculateContrasts(
    nst_single,
    comparisons = "A-B",
    condCol = "group",
    type = "welch",
    leastRepCount = 0,
    impute = TRUE,
    imputeMinFraction = 1,
    subsetByComparison = TRUE
  )
  annot_single <- generateAnnotatedMatrix(out_single)

  expect_true("featureAvg" %in% colnames(annot_single))
  expect_equal(
    annot_single[["featureAvg"]],
    pairwiseCompsAve(out_single)[["A-B"]]
  )
})

test_that("generateAnnotatedMatrix validates and applies custom comparison labels", {
  test_data <- matrix(
    c(
      10,
      11,
      NA,
      NA,
      NA,
      NA
    ),
    nrow = 1,
    byrow = TRUE
  )
  colnames(test_data) <- c("A1", "A2", "B1", "B2", "C1", "C2")

  test_df <- data.frame(
    sample = colnames(test_data),
    group = c(rep("A", 2), rep("B", 2), rep("C", 2))
  )
  rownames(test_df) <- test_df$sample

  se <- SummarizedExperiment::SummarizedExperiment(
    assay = test_data,
    colData = test_df,
    rowData = data.frame(annot = "Pep1")
  )

  nst <- NormalyzerStatistics(se)
  out <- calculateContrasts(
    nst,
    comparisons = c("A-B", "B-C"),
    condCol = "group",
    type = "welch",
    leastRepCount = 0,
    impute = TRUE,
    imputeMinFraction = 1,
    subsetByComparison = TRUE
  )

  expect_error(
    generateAnnotatedMatrix(out, compLabels = "only_one"),
    class = "normalyzerde_error"
  )

  annot <- generateAnnotatedMatrix(
    out,
    prefixSep = ".",
    compLabels = c("first", "second")
  )

  expect_true(all(c("first.PValue", "second.PValue") %in% colnames(annot)))
  expect_true(all(
    c("first.featureAvg", "second.featureAvg") %in% colnames(annot)
  ))
})

nd_make_one_vs_rest_se <- function(test_data, groups, batch = NULL) {
  design <- data.frame(
    sample = colnames(test_data),
    group = groups,
    stringsAsFactors = FALSE
  )

  if (!is.null(batch)) {
    design$batch <- batch
  }

  rownames(design) <- design$sample

  SummarizedExperiment::SummarizedExperiment(
    assay = test_data,
    colData = design,
    rowData = data.frame(annot = paste0("Pep", seq_len(nrow(test_data))))
  )
}

test_that("calculateContrasts_oneVsRest_welch", {
  test_data <- matrix(
    c(
      1,
      1,
      2,
      2,
      3,
      3,
      0,
      1,
      0,
      1,
      1,
      2
    ),
    nrow = 2,
    byrow = TRUE
  )
  colnames(test_data) <- c("A1", "A2", "B1", "B2", "C1", "C2")

  se <- nd_make_one_vs_rest_se(
    test_data = test_data,
    groups = c(rep("A", 2), rep("B", 2), rep("C", 2))
  )

  nst <- NormalyzerStatistics(se)
  out <- calculateContrasts(
    nst,
    condCol = "group",
    type = "welch",
    oneVsRest = TRUE
  )

  folds <- pairwiseCompsFold(out)
  expect_true(all(c("A-rest", "B-rest", "C-rest") %in% names(folds)))

  expect_true(all.equal(-1.5, round(folds[["A-rest"]][1], 5)))
  expect_true(all.equal(0, round(folds[["B-rest"]][1], 5)))
  expect_true(all.equal(1.5, round(folds[["C-rest"]][1], 5)))

  expect_true(all.equal(-0.5, round(folds[["A-rest"]][2], 5)))
  expect_true(all.equal(-0.5, round(folds[["B-rest"]][2], 5)))
  expect_true(all.equal(1, round(folds[["C-rest"]][2], 5)))
})

test_that("calculateContrasts_oneVsRest_welch_rejects_batch", {
  test_data <- matrix(
    c(
      1,
      1,
      2,
      2
    ),
    nrow = 1,
    byrow = TRUE
  )
  colnames(test_data) <- c("A1", "A2", "B1", "B2")

  se <- nd_make_one_vs_rest_se(
    test_data = test_data,
    groups = c(rep("A", 2), rep("B", 2)),
    batch = c("b1", "b2", "b1", "b2")
  )

  nst <- NormalyzerStatistics(se)
  expect_error(
    calculateContrasts(
      nst,
      condCol = "group",
      batchCol = "batch",
      type = "welch",
      oneVsRest = TRUE
    ),
    class = "normalyzerde_error"
  )
})

test_that("calculateContrasts_oneVsRestGroups", {
  test_data <- matrix(
    c(
      1,
      1,
      2,
      2,
      3,
      3
    ),
    nrow = 1,
    byrow = TRUE
  )
  colnames(test_data) <- c("A1", "A2", "B1", "B2", "C1", "C2")

  se <- nd_make_one_vs_rest_se(
    test_data = test_data,
    groups = c(rep("A", 2), rep("B", 2), rep("C", 2))
  )

  nst <- NormalyzerStatistics(se)
  out <- calculateContrasts(
    nst,
    condCol = "group",
    type = "welch",
    oneVsRest = TRUE,
    oneVsRestGroups = c("B", "C")
  )

  folds <- pairwiseCompsFold(out)
  expect_true(all(c("B-rest", "C-rest") %in% names(folds)))
  expect_true(!("A-rest" %in% names(folds)))
})

test_that("calculateContrasts errors when oneVsRestGroups is empty", {
  test_data <- matrix(
    c(
      1,
      1,
      2,
      2
    ),
    nrow = 1,
    byrow = TRUE
  )
  colnames(test_data) <- c("A1", "A2", "B1", "B2")

  se <- nd_make_one_vs_rest_se(
    test_data = test_data,
    groups = c(rep("A", 2), rep("B", 2))
  )

  nst <- NormalyzerStatistics(se)
  expect_error(
    calculateContrasts(
      nst,
      condCol = "group",
      type = "welch",
      oneVsRest = TRUE,
      oneVsRestGroups = character()
    ),
    class = "normalyzerde_error"
  )
})

test_that("calculateContrasts errors for unknown statistics type", {
  test_data <- matrix(
    c(
      1,
      1,
      2,
      2
    ),
    nrow = 1,
    byrow = TRUE
  )
  colnames(test_data) <- c("A1", "A2", "B1", "B2")

  test_df <- data.frame(
    sample = colnames(test_data),
    group = c(rep("A", 2), rep("B", 2)),
    stringsAsFactors = FALSE
  )
  rownames(test_df) <- test_df$sample

  se <- SummarizedExperiment::SummarizedExperiment(
    assay = test_data,
    colData = test_df,
    rowData = data.frame(annot = "Pep1")
  )

  nst <- NormalyzerStatistics(se)
  expect_error(
    calculateContrasts(
      nst,
      comparisons = "A-B",
      condCol = "group",
      type = "mystery"
    ),
    class = "normalyzerde_error"
  )
})

test_that("calculateContrasts one-vs-rest errors when filtering removes all rows", {
  test_data <- matrix(
    c(
      NA,
      NA,
      2,
      2
    ),
    nrow = 1,
    byrow = TRUE
  )
  colnames(test_data) <- c("A1", "A2", "B1", "B2")

  test_df <- data.frame(
    sample = colnames(test_data),
    group = c(rep("A", 2), rep("B", 2)),
    stringsAsFactors = FALSE
  )
  rownames(test_df) <- test_df$sample

  se <- SummarizedExperiment::SummarizedExperiment(
    assay = test_data,
    colData = test_df,
    rowData = data.frame(annot = "Pep1")
  )

  nst <- NormalyzerStatistics(se)
  expect_error(
    calculateContrasts(
      nst,
      condCol = "group",
      type = "welch",
      leastRepCount = 1,
      oneVsRest = TRUE,
      oneVsRestGroups = "A"
    ),
    class = "normalyzerde_error"
  )
})

test_that("calculateContrasts subsetByComparison errors when filtering removes all rows", {
  test_data <- matrix(
    c(
      NA,
      NA,
      2,
      2
    ),
    nrow = 1,
    byrow = TRUE
  )
  colnames(test_data) <- c("A1", "A2", "B1", "B2")

  test_df <- data.frame(
    sample = colnames(test_data),
    group = c(rep("A", 2), rep("B", 2)),
    stringsAsFactors = FALSE
  )
  rownames(test_df) <- test_df$sample

  se <- SummarizedExperiment::SummarizedExperiment(
    assay = test_data,
    colData = test_df,
    rowData = data.frame(annot = "Pep1")
  )

  nst <- NormalyzerStatistics(se)
  expect_error(
    calculateContrasts(
      nst,
      comparisons = "A-B",
      condCol = "group",
      type = "welch",
      leastRepCount = 1,
      subsetByComparison = TRUE
    ),
    class = "normalyzerde_error"
  )
})

test_that("calculateContrasts_oneVsRest_impute_limma", {
  test_data <- matrix(
    c(
      NA,
      NA,
      1,
      1,
      2,
      2
    ),
    nrow = 1,
    byrow = TRUE
  )
  colnames(test_data) <- c("A1", "A2", "B1", "B2", "C1", "C2")

  test_df <- data.frame(
    sample = colnames(test_data),
    group = c(rep("A", 2), rep("B", 2), rep("C", 2))
  )
  rownames(test_df) <- test_df$sample

  se <- SummarizedExperiment::SummarizedExperiment(
    assay = test_data,
    colData = test_df,
    rowData = data.frame(annot = "Pep1")
  )

  nst <- NormalyzerStatistics(se)
  out <- suppressWarnings(
    calculateContrasts(
      nst,
      condCol = "group",
      type = "limma",
      leastRepCount = 0,
      impute = TRUE,
      imputeMinFraction = 1,
      oneVsRest = TRUE,
      oneVsRestGroups = c("A")
    )
  )

  folds <- pairwiseCompsFold(out)
  expect_true(all.equal(-0.5, round(folds[["A-rest"]][1], 5)))
})

test_that("calculateContrasts_oneVsRest_limma_supports_batch", {
  set.seed(1)
  mat <- matrix(stats::rnorm(20 * 6, mean = 10, sd = 1), nrow = 20)
  colnames(mat) <- c("A1", "A2", "B1", "B2", "C1", "C2")

  design <- data.frame(
    sample = colnames(mat),
    group = c(rep("A", 2), rep("B", 2), rep("C", 2)),
    batch = rep(c("b1", "b2"), 3),
    stringsAsFactors = FALSE
  )
  rownames(design) <- design$sample

  se <- SummarizedExperiment::SummarizedExperiment(
    assay = mat,
    colData = design,
    rowData = data.frame(annot = paste0("Pep", seq_len(nrow(mat))))
  )

  nst <- NormalyzerStatistics(se)
  out <- suppressWarnings(calculateContrasts(
    nst,
    condCol = "group",
    batchCol = "batch",
    type = "limma",
    leastRepCount = 1,
    oneVsRest = TRUE
  ))

  restLabel <- NormalyzerDE:::chooseOneVsRestLabel(design$group)
  expected <- paste0(c("A", "B", "C"), "-", restLabel)

  folds <- pairwiseCompsFold(out)
  expect_true(all(expected %in% names(folds)))
  expect_equal(length(folds[[expected[1]]]), nrow(mat))
})

test_that("calculateContrasts_limma_allows_group_names_with_spaces", {
  test_data <- matrix(
    c(
      1,
      1,
      2,
      2,
      0,
      1,
      0,
      1
    ),
    nrow = 2,
    byrow = TRUE
  )
  colnames(test_data) <- c("s1", "s2", "s3", "s4")

  test_df <- data.frame(
    sample = colnames(test_data),
    group = c(rep("Group A", 2), rep("Group B", 2))
  )
  rownames(test_df) <- test_df$sample

  se <- SummarizedExperiment::SummarizedExperiment(
    assay = test_data,
    colData = test_df,
    rowData = data.frame(annot = paste0("Pep", seq_len(2)))
  )

  nst <- NormalyzerStatistics(se)
  out <- calculateContrasts(
    nst,
    comparisons = c("Group A-Group B"),
    condCol = "group",
    type = "limma"
  )

  folds <- pairwiseCompsFold(out)
  expect_true("Group A-Group B" %in% names(folds))
  expect_true(all.equal(-1, round(folds[["Group A-Group B"]][1], 5)))
  expect_true(all.equal(0, round(folds[["Group A-Group B"]][2], 5)))
})

test_that("calculateContrasts_limma_allows_group_names_starting_with_special_characters", {
  test_data <- matrix(
    c(
      1,
      1,
      2,
      2
    ),
    nrow = 1,
    byrow = TRUE
  )
  colnames(test_data) <- c("s1", "s2", "s3", "s4")

  test_df <- data.frame(
    sample = colnames(test_data),
    group = c(rep("#A", 2), rep("B", 2))
  )
  rownames(test_df) <- test_df$sample

  se <- SummarizedExperiment::SummarizedExperiment(
    assay = test_data,
    colData = test_df,
    rowData = data.frame(annot = "Pep1")
  )

  nst <- NormalyzerStatistics(se)
  out <- calculateContrasts(
    nst,
    comparisons = c("#A-B"),
    condCol = "group",
    type = "limma"
  )

  folds <- pairwiseCompsFold(out)
  expect_true("#A-B" %in% names(folds))
  expect_true(all.equal(-1, round(folds[["#A-B"]][1], 5)))
})

test_that("calculateContrasts_oneVsRest_allows_group_named_rest", {
  test_data <- matrix(
    c(
      1,
      1,
      2,
      2
    ),
    nrow = 1,
    byrow = TRUE
  )
  colnames(test_data) <- c("r1", "r2", "a1", "a2")

  test_df <- data.frame(
    sample = colnames(test_data),
    group = c(rep("rest", 2), rep("A", 2))
  )
  rownames(test_df) <- test_df$sample

  se <- SummarizedExperiment::SummarizedExperiment(
    assay = test_data,
    colData = test_df,
    rowData = data.frame(annot = "Pep1")
  )

  nst <- NormalyzerStatistics(se)
  out <- calculateContrasts(
    nst,
    condCol = "group",
    type = "limma",
    oneVsRest = TRUE
  )

  folds <- pairwiseCompsFold(out)
  expect_true(all(c("rest-others", "A-others") %in% names(folds)))
  expect_true(all.equal(-1, round(folds[["rest-others"]][1], 5)))
  expect_true(all.equal(1, round(folds[["A-others"]][1], 5)))
})
