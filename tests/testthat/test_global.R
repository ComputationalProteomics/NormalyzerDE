context("Global test")

designPath <- system.file(package="NormalyzerDE", "extdata", "tiny_design.tsv")
designSingleRepPath <- system.file(package="NormalyzerDE", "extdata", "tiny_design_singlerep.tsv")
dataPath <- system.file(package="NormalyzerDE", "extdata", "tiny_data.tsv")
# designPath <- "../../vignettes/design.tsv"
# dataPath <- "../../vignettes/data.tsv"
tempOut <- tempdir()

referenceNormResultsDir <- system.file(package="NormalyzerDE", "extdata", "unit_test_run_norm_reference")
referenceStatResultsDir <- system.file(package="NormalyzerDE", "extdata", "unit_test_run_stat_reference")

# Note: During testing the expect_silent uses sink() to capture all output
# For debugging using the browser() statement while in sink - run
# closeAllConnections()

are_matrices_identical <- function(label, samples, path1, path2) {
    
    # browser()
    
    df1 <- read.csv(path1, sep="\t")
    df2 <- read.csv(path2, sep="\t")
    
    expect_true(
        all(dim(df1) == dim(df2)), 
        paste("Expected identical dimensions for", label, 
              "found:", paste(dim(df1), collapse=", "), "and", paste(dim(df2), collapse=", "))
    )
    
    expect_true(
        all(colnames(df1) == colnames(df2)),
        paste("Expected identical column names for", label,
              "found:", paste(colnames(df1), collapse=", "), "and", paste(colnames(df2), collapse=", "))
    )

    df1colSums <- round(colSums(df1[, samples], na.rm=TRUE), 4)
    df2colSums <- round(colSums(df2[, samples], na.rm=TRUE), 4)

    expect_true(
        all(df1colSums == df2colSums), 
        paste("Expected identical column sums for", label, 
              "found:", paste(df1colSums, collapse=", "), "and", paste(df2colSums, collapse=", "))
    )
    
    df1annots <- df1[, !(colnames(df1) %in% samples)]
    df2annots <- df2[, !(colnames(df2) %in% samples)]
    
    expect_true(
        all(df1annots == df2annots), 
        paste("Expected identical column sums for", label, 
              "found:", paste(df1colSums, collapse=", "), "and", paste(df2colSums, collapse=", "))
    )
}

compare_output_directories <- function(label, samples, dir1, dir2) {
    
    currFiles <- list.files(dir1 , pattern="*.txt", full.names=TRUE)
    refFiles <- list.files(dir2, pattern="*.txt", full.names=TRUE)
    expect_true(length(currFiles) == length(refFiles), "Number of normalized matrices should be same as reference")
    
    designDf <- read.csv(designPath, sep="\t")
    samples <- as.character(designDf$sample)
    
    for (i in seq_len(length(currFiles))) {
        label <- basename(refFiles[[1]])
        are_matrices_identical(label, samples, currFiles[[i]], refFiles[[i]])
    }
}

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

test_that("Normalization results are identical to previous", {

    designDf <- read.csv(designPath, sep="\t")
    samples <- as.character(designDf$sample)
    currOutDir <- paste0(tempOut, "/unit_test_run_norm")
    compare_output_directories("NormalRun", samples, currOutDir, referenceNormResultsDir)
})

test_that("Normalization run succeeds without errors (single replicates)", {

    expect_silent(
        normalyzer(
            jobName="unit_test_run_norm_singlerep",
            dataPath=dataPath,
            designPath=designSingleRepPath,
            outputDir=tempOut,
            quiet=TRUE,
            requireReplicates=FALSE
        )
    )
})

test_that("Single-replicate run output has identical normalizations", {
    
    designDf <- read.csv(designPath, sep="\t")
    samples <- as.character(designDf$sample)
    currOutDir <- paste0(tempOut, "/unit_test_run_norm_singlerep")
    compare_output_directories("SingleRepRun", samples, currOutDir, referenceNormResultsDir)
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

test_that("Statistics with custom thresholding", {
    
    expect_silent(
        normalyzerDE(
            jobName="unit_test_run_stat_custom_thres",
            dataPath=dataPath,
            designPath=designPath,
            outputDir=tempOut,
            comparisons=c("4-5"),
            logTrans=TRUE,
            quiet=TRUE,
            sigThres=0.05,
            sigThresType="p",
            log2FoldThres="1"
        )   
    )
})

test_that("Statistics results are identical to previous", {
    
    currOutDir <- paste0(tempOut, "/unit_test_run_stat")
    currStatFiles <- list.files(currOutDir , pattern="*.tsv", full.names=TRUE)
    refStatFiles <- list.files(referenceStatResultsDir, pattern="*.tsv", full.names=TRUE)
    expect_true(length(currStatFiles) == length(refStatFiles), "Number of normalized matrices should be same as reference")
    
    currStatMd5 <- tools::md5sum(currStatFiles)
    refStatMd5 <- tools::md5sum(refStatFiles)
    
    expect_true(all(currStatMd5 == refStatMd5), "MD5-sums for stat matrices should be equal")
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
