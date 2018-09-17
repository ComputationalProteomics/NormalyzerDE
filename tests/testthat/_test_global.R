context("Global test")

designPath <- system.file(package="NormalyzerDE", "extdata", "tiny_design.tsv")
dataPath <- system.file(package="NormalyzerDE", "extdata", "tiny_data.tsv")
# designPath <- "../../vignettes/design.tsv"
# dataPath <- "../../vignettes/data.tsv"
tempOut <- tempdir()

# Note: During testing the expect_silent uses sink() to capture all output
# For debugging using the browser() statement while in sink - run
# closeAllConnections()

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
            comparisons=c("4-5"),
            logTrans=TRUE,
            quiet=TRUE
        )   
    )
})

test_that("Normalization run from SummarizedExperiment succeeds without errors", {
    
    rdf <- read.csv(dataPath, sep="\t")
    ddf <- read.csv(designPath, sep="\t")
    ddf$sample <- as.character(ddf$sample)
    sdf <- rdf[, ddf$sample]
    adf <- rdf[, !(colnames(rdf) %in% ddf$sample)]
    se <- SummarizedExperiment::SummarizedExperiment(
        assays=list(raw=as.matrix(sdf)),
        colData=ddf,
        rowData=adf
    )
    
    expect_silent(
        normalyzer(
            jobName="unit_test_run_norm",
            experimentObj=se,
            outputDir=tempOut, 
            quiet=TRUE
        )
    )
})

test_that("Statistics run from SummarizedExperiment succeeds without errors", {
    
    rdf <- read.csv(dataPath, sep="\t")
    ddf <- read.csv(designPath, sep="\t")
    ddf$sample <- as.character(ddf$sample)
    sdf <- rdf[, ddf$sample]
    adf <- rdf[, !(colnames(rdf) %in% ddf$sample)]
    se <- SummarizedExperiment::SummarizedExperiment(
        assays=list(raw=as.matrix(sdf)),
        colData=ddf,
        rowData=adf
    )
    
    expect_silent(
        normalyzerDE(
            jobName="unit_test_run_stat",
            experimentObj=se,
            outputDir=tempOut,
            comparisons=c("4-5"),
            logTrans=TRUE,
            quiet=TRUE
        )   
    )
})
