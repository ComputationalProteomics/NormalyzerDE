context("calculateStatistics")

# Transfer to utils
test_that("filterLowRep", {

    test_data <- data.frame(
        c(NA, 1,  1,  NA), 
        c(1,  NA, 1,  NA), 
        c(3,  3,  3,  0), 
        c(5,  3,  NA, 0), 
        c(NA, 5,  NA, 0), 
        c(7,  9,  NA, 0))
    colnames(test_data) <- c("a1", "a2", "a3", "b1", "b2", "b3")
    groups <- c(rep("A", 3), rep("B", 3))
    
    expected_out_data <- data.frame(
        "a1"= c(NA, 1), 
        "a2"=c(1,  NA), 
        "a3"=c(3,  3), 
        "b1"=c(5,  3), 
        "b2"=c(NA, 5), 
        "b3"=c(7,  9))
    
    out <- filterLowRep(test_data, groups, leastRep=2)
    
    expect_that(
        all.equal(
            expected_out_data,
            out
        ),
        is_true()
    )
})

test_that("reduceTechnicalReplicates", {
    
    tech_rep <- c("a", "a", "b", "b", "c", "c", "d", "d")
    test_data <- data.frame(
        c(1,1,1), 
        c(1,2,1), 
        c(3,3,3), 
        c(5,3,3), 
        c(5,5,4), 
        c(5,5,5), 
        c(7,7,7), 
        c(7,9,7))
    colnames(test_data) <- c("a1", "a2", "b1", "b2", "c1", "c2", "d1", "d2")
    
    expect_out_data <- as.matrix(data.frame(
        "a"=c(1,1.5,1), 
        "b"=c(4,3,3), 
        "c"=c(5,5,4.5), 
        "d"=c(7,8,7)))
    
    out <- reduceTechnicalReplicates(test_data, tech_rep)
    
    expect_that(
        all.equal(
            expect_out_data,
            out
        ),
        is_true()
    )
})

test_that("reduceDesignTechRep", {
    
    test_df <- data.frame(
        sample=c("a1", "a2", "a3", "b1", "b2", "c1", "c2", "d1"),
        group=c(rep("A", 5), rep("B", 3)),
        techrep=c("a", "a", "a", "b", "b", "c", "c", "d")
    )
    
    expected_out_df <- data.frame(
        sample=c("a1", "b1", "c1", "d1"),
        group=c(rep("A", 2), rep("B", 2)),
        techrep=c("a", "b", "c", "d")
    )
    
    out <- reduceDesignTechRep(test_df, test_df$techrep)
    
    expect_that(
        all.equal(
            expected_out_df,
            out
        ),
        is_true()
    )
})

# Statistics
test_that("setupStatisticsObject", {
    
    # expected_out <- regression_test_ner@avgvarmem
    # out <- calculateAvgReplicateVariation(normMatrices, sampleReplicateGroups)
    
    # expect_that(
    #     all.equal(
    #         expected_out,
    #         out
    #     ),
    #     is_true()
    # )
})

test_that("generateAnnotatedMatrix", {
    
    # expected_out <- regression_test_ner@avgvarmem
    # out <- calculateAvgReplicateVariation(normMatrices, sampleReplicateGroups)
    
    # expect_that(
    #     all.equal(
    #         expected_out,
    #         out
    #     ),
    #     is_true()
    # )
})

test_that("generateStatsReport", {
    
    # expected_out <- regression_test_ner@avgvarmem
    # out <- calculateAvgReplicateVariation(normMatrices, sampleReplicateGroups)
    
    # expect_that(
    #     all.equal(
    #         expected_out,
    #         out
    #     ),
    #     is_true()
    # )
})

test_that("calculateWelch", {
    
    # expected_out <- regression_test_ner@avgvarmem
    # out <- calculateAvgReplicateVariation(normMatrices, sampleReplicateGroups)
    
    # expect_that(
    #     all.equal(
    #         expected_out,
    #         out
    #     ),
    #     is_true()
    # )
})

test_that("calculateLimmaContrast", {
    
    # expected_out <- regression_test_ner@avgvarmem
    # out <- calculateAvgReplicateVariation(normMatrices, sampleReplicateGroups)
    
    # expect_that(
    #     all.equal(
    #         expected_out,
    #         out
    #     ),
    #     is_true()
    # )
})


