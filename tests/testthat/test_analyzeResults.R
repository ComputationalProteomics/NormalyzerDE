context("analyzeResults.R")

data("example_data_only_values")
data("example_design")
test_design <- example_design[example_design$group %in% c("1", "2", "3"), ]
test_data <- example_data_only_values[, as.character(test_design$sample)]
group_header <- test_design$group
unique_groups <- unique(group_header)

# Setup equivalent to regression_test_nr
data("example_wide_data")
data("example_wide_design")
sdf <- example_wide_data[, example_wide_design$sample]
adf <- example_wide_data[, !colnames(example_wide_data) %in% example_wide_design$sample]
se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(raw=as.matrix(sdf)),
    colData=example_wide_design,
    rowData=adf
)
SummarizedExperiment::metadata(se) <- list(
    sample="sample",
    group="group"
)

nds <- NormalyzerDE::getVerifiedNormalyzerObject(
    "unit_test_run_norm",
    summarizedExp = se
)
regression_test_nr <- NormalyzerDE::normMethods(nds)
regression_test_nr <- NormalyzerDE::analyzeNormalizations(regression_test_nr)
regression_test_ner <- regression_test_nr@ner

#data("regression_test_nr")
#regression_test_ner <- ner(regression_test_nr)
#nds <- nds(regression_test_nr)

# browser()

sampleReplicateGroups <- sampleReplicateGroups(nds)
normMatrices <- normalizations(regression_test_nr)
anova_pvalues <- calculateANOVAPValues(
    normMatrices,
    sampleReplicateGroups,
    categoricalANOVA=FALSE)

test_that("calculateCorrSum gives same Pearson output", {

    expected_group_pearson <- c(
        0.9890682, 0.9885302, 0.9972236, 
        0.9937719, 0.9976159, 0.9948526, 
        0.6108868, 0.9951509, 0.5691397
    )
    
    pear_out <- calculateCorrSum(test_data, group_header, unique_groups, "pearson")
    expect_true(
        all.equal(
            expected_group_pearson,
            round(pear_out, 7)
        )
    )
})

test_that("calculateCorrSum gives same Spearman output", {

    expected_group_spearman <- c(
        0.9234304, 0.9665739, 0.9497886,
        0.9686747, 0.9551886, 0.9420609, 
        0.9250426, 0.9549601, 0.9315987
    )
    
    spear_out <- calculateCorrSum(test_data, group_header, unique_groups, "spearman")
    expect_true(
        all.equal(
            expected_group_spearman,
            round(spear_out, 7)
        )
    )
})

test_that("calculateSummarizedCorrelationVector_Spearman", {
    
    # expected_out <- repCorSpear(regression_test_ner)
    expected_out <- read.csv(system.file(package="NormalyzerDE", "testdata", "calculateSummarizedCorrelationVector_spearman.csv"), check.names = FALSE)
    out <- calculateSummarizedCorrelationVector(
        normMatrices, 
        sampleReplicateGroups,
        unique(sampleReplicateGroups),
        "spearman"
    )
    
    expect_true(
        all.equal(
            round(as.matrix(expected_out), 5),
            round(out, 5)
        )
    )
})

test_that("calculateSummarizedCorrelationVector_Pearson", {
    
    #expected_out <- repCorPear(regression_test_ner)
    expected_out <- read.csv(system.file(package="NormalyzerDE", "testdata", "calculateSummarizedCorrelationVector_pearson.csv"), check.names = FALSE)
    out <- as.matrix(unlist(calculateSummarizedCorrelationVector(
        normMatrices, 
        sampleReplicateGroups,
        unique(sampleReplicateGroups),
        "pearson"
    )))

    expect_true(
        all.equal(
            round(as.matrix(expected_out), 5),
            round(out, 5)
        )
    )
})

test_that("calculateReplicateCV", {
    
    expected_out <- avgcvmem(regression_test_ner)
    out <- calculateReplicateCV(normMatrices, sampleReplicateGroups)

    expect_true(
        all.equal(
            expected_out,
            out
        )
    )
})

test_that("calculateFeatureCV", {
    
    expected_out <- featureCVPerMethod(regression_test_ner)
    out <- calculateFeatureCV(normMatrices)
    
    expect_true(
        all.equal(
            expected_out,
            out
        )
    )
})

test_that("calculateAvgMadMem", {
    
    expected_out <- avgmadmem(regression_test_ner)
    out <- calculateAvgMadMem(normMatrices, sampleReplicateGroups)

    expect_true(
        all.equal(
            expected_out,
            out
        )
    )
})

test_that("calculateAvgReplicateVariation", {
    
    expected_out <- avgvarmem(regression_test_ner)
    out <- calculateAvgReplicateVariation(normMatrices, sampleReplicateGroups)

    expect_true(
        all.equal(
            expected_out,
            out
        )
    )
})


test_that("calculateANOVAPValues", {
    
    # expected_out <- anovaP(regression_test_ner)
    expected_out <- read.csv(system.file(package="NormalyzerDE", "testdata", "anova.csv"), check.names = FALSE)
    anova_pvalues_copy <- anova_pvalues

    expect_true(
        all.equal(
            round(as.matrix(expected_out), 4),
            round(anova_pvalues_copy, 4)
        )
    )
})

test_that("findLowlyVariableFeaturesCVs", {
    
    # expected_out <- lowVarFeaturesCVs(regression_test_ner)
    expected_out <- c(
        "log2" = 1.153562, 
        "VSN" = 1.089077, 
        "GI" = 1.159375, 
        "median" = 1.320822, 
        "mean" = 1.163794, 
        "Quantile" = 1.120080, 
        "CycLoess" = 1.032822, 
        "RLR" = 1.048589,
        "RT-median" = 1.320822, 
        "RT-mean" = 1.163794, 
        "RT-Loess" = 1.032822, 
        "RT-VSN" = 1.089077
    )
    
    anova_pvalues_contrast_nonna <- !is.na(anova_pvalues[, 1])
    log2AnovaFDR <- rep(NA, length(anova_pvalues_contrast_nonna))
    log2AnovaFDR[anova_pvalues_contrast_nonna] <- stats::p.adjust(
        anova_pvalues[, 1][anova_pvalues_contrast_nonna], 
        method="BH")
    
    out <- findLowlyVariableFeaturesCVs(log2AnovaFDR, normMatrices)

    expect_true(
        all.equal(
            round(expected_out, 6),
            round(out, 6)
        )
    )
})

test_that("calculatePercentageAvgDiffInMat_small_test", {
    
    test_mat <- data.frame(c(1,3), c(1,5), c(1,7))
    expect_out <- c(100, 150, 200)
    out <- calculatePercentageAvgDiffInMat(test_mat)
    
    expect_true(
        all.equal(
            out,
            expect_out
        )
    )
})

test_that("calculatePercentageAvgDiffInMat_MAD", {
    
    expected_out <- avgmadmempdiff(regression_test_ner)
    avgMadMemMat <- calculateAvgMadMem(normMatrices, sampleReplicateGroups)
    out <- calculatePercentageAvgDiffInMat(avgMadMemMat)
    
    expect_true(
        all.equal(
            expected_out,
            out
        )
    )
})

test_that("calculatePercentageAvgDiffInMat_CV", {
    
    expected_out <- avgcvmempdiff(regression_test_ner)
    avgCVMat <- calculateReplicateCV(normMatrices, sampleReplicateGroups)
    out <- calculatePercentageAvgDiffInMat(avgCVMat)
    
    expect_true(
        all.equal(
            expected_out,
            out
        )
    )
})



