context("Global test")

designPath <- "../../vignettes/design.tsv"
dataPath <- "../../vignettes/data.tsv"
tempOut <- tempdir()

test_that("Normalization run succeeds without errors", {
    
    expect_silent(
        normalyzer(
            jobName="unit_test_run_norm",
            dataPath=dataPath,
            designPath=designPath,
            outputDir=tempOut, 
            quiet=TRUE
        )
    )
})

test_that("Statistics run succeeds without errors", {
    
    expect_silent(
        normalyzerDE(
            jobName="unit_test_run_stat",
            dataPath=dataPath,
            designPath=designPath,
            outputDir=tempOut,
            comparisons=c("1-2", "2-3"),
            logTrans=TRUE,
            quiet=TRUE
        )   
    )
})
