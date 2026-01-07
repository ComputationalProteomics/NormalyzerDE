context("technical replicates in DE")

test_that("normalyzerDE can reduce technical replicates before contrasts", {

    mat <- matrix(stats::rnorm(20 * 8, mean=10, sd=1), nrow=20)
    colnames(mat) <- paste0("s", seq_len(ncol(mat)))

    design <- data.frame(
        sample=colnames(mat),
        group=c("A", "A", "A", "A", "B", "B", "B", "B"),
        techrep=c("t1", "t1", "t2", "t2", "t3", "t3", "t4", "t4"),
        stringsAsFactors=FALSE
    )
    rownames(design) <- design$sample

    se <- SummarizedExperiment::SummarizedExperiment(
        assay=mat,
        colData=design,
        rowData=data.frame(feature=paste0("f", seq_len(nrow(mat))))
    )

    jobName <- paste0("techrep_de_", sample.int(1e9, 1))
    outDir <- tempdir()
    expectedDir <- file.path(outDir, NormalyzerDE:::sanitizeJobName(jobName))
    if (dir.exists(expectedDir)) unlink(expectedDir, recursive=TRUE, force=TRUE)
    on.exit(unlink(expectedDir, recursive=TRUE, force=TRUE), add=TRUE)

    out <- suppressWarnings(normalyzerDE(
        jobName=jobName,
        comparisons="A-B",
        experimentObj=se,
        sampleCol="sample",
        condCol="group",
        techRepCol="techrep",
        type="limma",
        outputDir=outDir,
        quiet=TRUE
    ))

    expect_null(out)

    outStatsPath <- file.path(expectedDir, paste0(basename(expectedDir), "_stats.tsv"))
    expect_true(file.exists(outStatsPath))

    outDf <- utils::read.table(outStatsPath, sep="\t", header=TRUE, check.names=FALSE)
    expectedCollapsedNames <- c("s1.s2", "s3.s4", "s5.s6", "s7.s8")
    expect_true(all(expectedCollapsedNames %in% colnames(outDf)))
    expect_false(any(paste0("s", seq_len(8)) %in% colnames(outDf)))
})
