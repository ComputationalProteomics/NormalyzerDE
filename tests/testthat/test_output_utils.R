context("outputUtils")

test_that("writeNormalizedDatasets can include anova and CV columns", {

    data(example_summarized_experiment)

    se <- example_summarized_experiment
    SummarizedExperiment::rowData(se) <- S4Vectors::DataFrame(row.names=seq_len(nrow(se)))

    normObj <- getVerifiedNormalyzerObject(
        jobName=paste0("oututils_", sample.int(1e9, 1)),
        summarizedExp=se,
        threshold=0,
        omitSamples=FALSE,
        requireReplicates=TRUE,
        quiet=TRUE
    )
    nr <- normMethods(normObj, normalizeRetentionTime=FALSE, quiet=TRUE)
    nr <- analyzeNormalizations(nr)

    outDir <- file.path(tempdir(), paste0("oututils_", sample.int(1e9, 1)))
    dir.create(outDir, recursive=TRUE)
    on.exit(unlink(outDir, recursive=TRUE, force=TRUE), add=TRUE)

    expect_silent(suppressWarnings(writeNormalizedDatasets(
        nr,
        jobdir=outDir,
        includeCvCol=TRUE,
        includeAnovaP=TRUE,
        includePairwiseComparisons=TRUE,
        normSuffix="-test-normalized.txt",
        rawdataName="raw_test.txt"
    )))

    methodNames <- names(normalizations(nr))
    expect_true(length(methodNames) > 0)

    filePath <- file.path(outDir, paste0(methodNames[1], "-test-normalized.txt"))
    expect_true(file.exists(filePath))

    outDf <- utils::read.table(filePath, sep="\t", header=TRUE, check.names=FALSE)
    expect_true("anovaP" %in% colnames(outDf))
    expect_true("CV" %in% colnames(outDf))
    expect_true("comp_4-5_p" %in% colnames(outDf))
    expect_true("comp_4-5_fdr" %in% colnames(outDf))
    expect_true(file.exists(file.path(outDir, "raw_test.txt")))
})
