context("analyzeResults.R")

# browser()

data("example_data_only_values")
data("example_design")
test_design <- example_design[example_design$group %in% c("1", "2", "3"), ]
test_data <- example_data_only_values[, as.character(test_design$sample)]
group_header <- test_design$group
unique_groups <- unique(group_header)

data("regression_test_nr")
regression_test_ner <- ner(regression_test_nr)
nds <- nds(regression_test_nr)

sampleReplicateGroups <- sampleReplicateGroups(nds)
normMatrices <- getNormalizationMatrices(regression_test_nr)
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
    expect_that(
        all.equal(
            expected_group_pearson,
            round(pear_out, 7)
        ),
        is_true()
    )
})

test_that("calculateCorrSum gives same Spearman output", {

    expected_group_spearman <- c(
        0.9234304, 0.9665739, 0.9497886,
        0.9686747, 0.9551886, 0.9420609, 
        0.9250426, 0.9549601, 0.9315987
    )
    
    spear_out <- calculateCorrSum(test_data, group_header, unique_groups, "spearman")
    expect_that(
        all.equal(
            expected_group_spearman,
            round(spear_out, 7)
        ),
        is_true()
    )
})

test_that("calculateSummarizedCorrelationVector_Spearman", {
    
    expected_out <- matrix(unlist(repCorSpear(regression_test_ner)), nrow=length(sampleReplicateGroups), byrow=FALSE)
    out <- calculateSummarizedCorrelationVector(
        normMatrices, 
        sampleReplicateGroups,
        unique(sampleReplicateGroups),
        "spearman"
    )
    dimnames(out) <- NULL
    
    expect_that(
        all.equal(
            expected_out,
            out
        ),
        is_true()
    )
})

test_that("calculateSummarizedCorrelationVector_Pearson", {
    
    expected_out <- matrix(unlist(repCorPear(regression_test_ner)), nrow=length(sampleReplicateGroups), byrow=FALSE)
    out <- as.matrix(unlist(calculateSummarizedCorrelationVector(
        normMatrices, 
        sampleReplicateGroups,
        unique(sampleReplicateGroups),
        "pearson"
    )))
    dimnames(out) <- NULL
    
    expect_that(
        all.equal(
            expected_out,
            out
        ),
        is_true()
    )
})

test_that("calculateReplicateCV", {
    
    expected_out <- avgcvmem(regression_test_ner)
    out <- calculateReplicateCV(normMatrices, sampleReplicateGroups)
    dimnames(out) <- NULL
    
    expect_that(
        all.equal(
            expected_out,
            out
        ),
        is_true()
    )
})

test_that("calculateFeatureCV", {
    
    expected_out <- featureCVPerMethod(regression_test_ner)
    out <- calculateFeatureCV(normMatrices)
    
    print("Comparing without header for now")
    dimnames(out) <- NULL
    
    expect_that(
        all.equal(
            expected_out,
            out
        ),
        is_true()
    )
})

test_that("calculateAvgMadMem", {
    
    expected_out <- avgmadmem(regression_test_ner)
    out <- calculateAvgMadMem(normMatrices, sampleReplicateGroups)
    dimnames(out) <- NULL
    
    expect_that(
        all.equal(
            expected_out,
            out
        ),
        is_true()
    )
})

test_that("calculateAvgReplicateVariation", {
    
    expected_out <- avgvarmem(regression_test_ner)
    out <- calculateAvgReplicateVariation(normMatrices, sampleReplicateGroups)
    dimnames(out) <- NULL
    
    expect_that(
        all.equal(
            expected_out,
            out
        ),
        is_true()
    )
})


test_that("calculateANOVAPValues", {
    
    expected_out <- anovaP(regression_test_ner)
    anova_pvalues_copy <- anova_pvalues
    dimnames(anova_pvalues_copy) <- NULL

    expect_that(
        all.equal(
            expected_out,
            anova_pvalues_copy
        ),
        is_true()
    )
}) 

test_that("findLowlyVariableFeaturesCVs", {
    
    expected_out <- lowVarFeaturesCVs(regression_test_ner)
    
    anova_pvalues_contrast_nonna <- !is.na(anova_pvalues[, 1])
    log2AnovaFDR <- rep(NA, length(anova_pvalues_contrast_nonna))
    log2AnovaFDR[anova_pvalues_contrast_nonna] <- stats::p.adjust(
        anova_pvalues[, 1][anova_pvalues_contrast_nonna], 
        method="BH")
    
    out <- findLowlyVariableFeaturesCVs(log2AnovaFDR, normMatrices)
    dimnames(out) <- NULL
    
    expect_that(
        all.equal(
            expected_out,
            out
        ),
        is_true()
    )
})

test_that("calculatePercentageAvgDiffInMat_small_test", {
    
    test_mat <- data.frame(c(1,3), c(1,5), c(1,7))
    expect_out <- c(100, 150, 200)
    out <- calculatePercentageAvgDiffInMat(test_mat)
    
    expect_that(
        all.equal(
            out,
            expect_out
        ),
        is_true()
    )
})

test_that("calculatePercentageAvgDiffInMat_MAD", {
    
    expected_out <- avgmadmempdiff(regression_test_ner)
    avgMadMemMat <- calculateAvgMadMem(normMatrices, sampleReplicateGroups)
    out <- calculatePercentageAvgDiffInMat(avgMadMemMat)
    
    expect_that(
        all.equal(
            expected_out,
            out
        ),
        is_true()
    )
})

test_that("calculatePercentageAvgDiffInMat_CV", {
    
    expected_out <- avgcvmempdiff(regression_test_ner)
    avgCVMat <- calculateReplicateCV(normMatrices, sampleReplicateGroups)
    out <- calculatePercentageAvgDiffInMat(avgCVMat)
    
    expect_that(
        all.equal(
            expected_out,
            out
        ),
        is_true()
    )
})



