context("entrypoints")

test_that("normalyzer writes outputs when given SummarizedExperiment", {

    data(example_summarized_experiment)

    jobName <- paste0("entry_norm_", sample.int(1e9, 1))
    outDir <- tempdir()
    expectedDir <- file.path(outDir, NormalyzerDE:::sanitizeJobName(jobName))
    if (dir.exists(expectedDir)) unlink(expectedDir, recursive=TRUE, force=TRUE)
    on.exit(unlink(expectedDir, recursive=TRUE, force=TRUE), add=TRUE)

    out <- suppressWarnings(normalyzer(
        jobName=jobName,
        experimentObj=example_summarized_experiment,
        outputDir=outDir,
        normalizeRetentionTime=FALSE,
        skipAnalysis=TRUE,
        quiet=TRUE
    ))

    expect_null(out)
    expect_true(dir.exists(expectedDir))
    expect_true(file.exists(file.path(expectedDir, "submitted_rawdata.txt")))
    expect_gt(length(list.files(expectedDir, pattern="-normalized\\.txt$", full.names=TRUE)), 0)
})

test_that("normalyzerDE can compute one-vs-rest without explicit comparisons", {

    mat <- matrix(stats::rnorm(20 * 4, mean=10, sd=1), nrow=20)
    colnames(mat) <- paste0("s", seq_len(ncol(mat)))

    design <- data.frame(
        sample=colnames(mat),
        group=c("A", "A", "B", "B"),
        stringsAsFactors=FALSE
    )
    rownames(design) <- design$sample

    se <- SummarizedExperiment::SummarizedExperiment(
        assay=mat,
        colData=design,
        rowData=data.frame(feature=paste0("f", seq_len(nrow(mat))))
    )

    jobName <- paste0("onevsrest_", sample.int(1e9, 1))
    outDir <- tempdir()
    expectedDir <- file.path(outDir, NormalyzerDE:::sanitizeJobName(jobName))
    if (dir.exists(expectedDir)) unlink(expectedDir, recursive=TRUE, force=TRUE)
    on.exit(unlink(expectedDir, recursive=TRUE, force=TRUE), add=TRUE)

    out <- suppressWarnings(normalyzerDE(
        jobName=jobName,
        comparisons=NULL,
        experimentObj=se,
        outputDir=outDir,
        type="limma",
        oneVsRest=TRUE,
        quiet=TRUE
    ))

    expect_null(out)

    outStatsPath <- file.path(expectedDir, paste0(basename(expectedDir), "_stats.tsv"))
    expect_true(file.exists(outStatsPath))

    restLabel <- NormalyzerDE:::chooseOneVsRestLabel(design$group)
    expectedPCols <- paste0(c("A", "B"), "-", restLabel, "_PValue")

    outDf <- utils::read.table(outStatsPath, sep="\t", header=TRUE, check.names=FALSE)
    expect_true(all(expectedPCols %in% colnames(outDf)))
})

