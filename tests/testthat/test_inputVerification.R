context("Global test")

designPath <- system.file(package="NormalyzerDE", "extdata", "tiny_design.tsv")
dataPath <- system.file(package="NormalyzerDE", "extdata", "tiny_data.tsv")
tempOut <- tempdir()

design <- read.csv(designPath, sep="\t")
sample_names <- as.character(design$sample)

rdf <- read.csv(
    dataPath, 
    stringsAsFactors=FALSE, 
    quote="", 
    comment.char="", 
    sep="\t",
    check.names=FALSE
)

sdf <- rdf[, sample_names]
adf <- rdf[, !(colnames(rdf) %in% sample_names)]

test_that("Parsed input data is same as raw input", {
    
    se <- setupRawDataObject(
        dataPath, 
        designPath,
        zeroToNA=TRUE,
        sampleColName="sample",
        groupColName="group",
        inputFormat="default"
    )
    
    expect_true(all(design == data.frame(SummarizedExperiment::colData(se))), "colData should be identical to data frame")

    se_data <- SummarizedExperiment::assay(se)
    class(se_data) <- "numeric"
    
    expect_true(all(colSums(se_data, na.rm=TRUE) == colSums(sdf, na.rm=TRUE)), "Sum of data should be same")
    expect_true(all(colnames(se_data) == colnames(sdf)), "Column names should be the same")
    expect_true(all(adf == data.frame(SummarizedExperiment::rowData(se))), "Annotation information should be the same")
})

