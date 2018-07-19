context("analyzeResults.R")

data("example_data_only_values")
data("example_design")
test_design <- example_design[which(example_design$group %in% c("1", "2", "3")), ]
test_data <- example_data_only_values[, as.character(test_design$sample)]
group_header <- test_design$group
unique_groups <- unique(group_header)

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


