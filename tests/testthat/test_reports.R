context("reports")

test_that("generateStatsReport writes a report for multiple contrasts", {

    data(example_stat_summarized_experiment)

    jobName <- paste0("de_report_", sample.int(1e9, 1))
    outDir <- file.path(tempdir(), paste0("report_", sample.int(1e9, 1)))
    dir.create(outDir, recursive=TRUE)
    on.exit(unlink(outDir, recursive=TRUE, force=TRUE), add=TRUE)

    nst <- NormalyzerStatistics(example_stat_summarized_experiment, logTrans=FALSE)
    nst <- calculateContrasts(
        nst,
        comparisons=c("1-2", "2-3"),
        condCol="group",
        type="limma"
    )

    safeJobName <- NormalyzerDE:::sanitizeJobName(jobName)
    pdfPath <- file.path(outDir, paste0("Norm-stats-report-", safeJobName, ".pdf"))
    pngDir <- file.path(outDir, "de_pngs")

    expect_silent(suppressWarnings(generateStatsReport(
        nst,
        jobName=jobName,
        jobDir=outDir,
        sigThres=0.5,
        sigThresType="fdr",
        log2FoldThres=0,
        writeAsPngs=TRUE
    )))

    expect_true(file.exists(pdfPath))
    expect_true(dir.exists(pngDir))
    expect_true(file.exists(file.path(pngDir, "1_front.png")))
    expect_gt(length(list.files(pngDir, pattern="\\.png$", full.names=TRUE)), 0)

    # Exercise the non-PNG branch too
    outDir2 <- file.path(tempdir(), paste0("report2_", sample.int(1e9, 1)))
    dir.create(outDir2, recursive=TRUE)
    on.exit(unlink(outDir2, recursive=TRUE, force=TRUE), add=TRUE)

    pdfPath2 <- file.path(outDir2, paste0("Norm-stats-report-", safeJobName, ".pdf"))
    pngDir2 <- file.path(outDir2, "de_pngs")

    expect_silent(suppressWarnings(generateStatsReport(
        nst,
        jobName=jobName,
        jobDir=outDir2,
        sigThres=0.5,
        sigThresType="fdr",
        log2FoldThres=0,
        writeAsPngs=FALSE
    )))
    expect_true(file.exists(pdfPath2))
    expect_false(dir.exists(pngDir2))
})

test_that("getSigs supports p-value and fold thresholds", {

    data(example_stat_summarized_experiment)

    nst <- NormalyzerStatistics(example_stat_summarized_experiment, logTrans=FALSE)
    nst <- calculateContrasts(
        nst,
        comparisons=c("1-2", "2-3"),
        condCol="group",
        type="limma"
    )

    sigs <- NormalyzerDE:::getSigs(
        nst,
        sigThresType="p",
        sigThres=0.5,
        log2FoldThres=0.1
    )
    expect_equal(length(sigs), length(comparisons(nst)))
    expect_true(all(vapply(sigs, length, integer(1)) == nrow(dataMat(nst))))
})

test_that("getSigs errors for unknown threshold type", {

    data(example_stat_summarized_experiment)

    nst <- NormalyzerStatistics(example_stat_summarized_experiment, logTrans=FALSE)
    nst <- calculateContrasts(
        nst,
        comparisons=c("1-2", "2-3"),
        condCol="group",
        type="limma"
    )

    expect_error(
        NormalyzerDE:::getSigs(nst, sigThresType="bad", sigThres=0.1, log2FoldThres=0),
        "Unknown significance threshold type"
    )
})

test_that("generatePlots writes a Normalyzer report", {

    data(example_summarized_experiment)

    jobName <- paste0("norm_report_", sample.int(1e9, 1))
    normObj <- getVerifiedNormalyzerObject(
        jobName=jobName,
        summarizedExp=example_summarized_experiment,
        threshold=0,
        omitSamples=FALSE,
        requireReplicates=TRUE,
        quiet=TRUE
    )
    nr <- suppressWarnings(normMethods(normObj, normalizeRetentionTime=FALSE, quiet=TRUE))
    nr <- suppressWarnings(analyzeNormalizations(nr))

    outDir <- file.path(tempdir(), paste0("normplots_", sample.int(1e9, 1)))
    dir.create(outDir, recursive=TRUE)
    on.exit(unlink(outDir, recursive=TRUE, force=TRUE), add=TRUE)

    safeJobName <- NormalyzerDE:::sanitizeJobName(jobName)
    pdfPath <- file.path(outDir, paste0("Norm-report-", safeJobName, ".pdf"))
    pngDir <- file.path(outDir, "pngs")

    expect_silent(suppressWarnings(generatePlots(nr, outDir, writeAsPngs=TRUE)))

    expect_true(file.exists(pdfPath))
    expect_true(dir.exists(pngDir))
    expect_true(file.exists(file.path(pngDir, "1_front_page.png")))
    expect_gt(length(list.files(pngDir, pattern="\\.png$", full.names=TRUE)), 0)

    # Exercise the non-PNG branch and tiny run mode too
    tiny_mat <- matrix(stats::rnorm(10 * 4, mean=10, sd=1), nrow=10)
    colnames(tiny_mat) <- paste0("s", seq_len(ncol(tiny_mat)))
    tiny_design <- data.frame(
        sample=colnames(tiny_mat),
        group=c("A", "A", "B", "B"),
        stringsAsFactors=FALSE
    )
    rownames(tiny_design) <- tiny_design$sample

    tiny_se <- SummarizedExperiment::SummarizedExperiment(
        assay=tiny_mat,
        colData=tiny_design,
        rowData=data.frame(feature=paste0("f", seq_len(nrow(tiny_mat))))
    )
    S4Vectors::metadata(tiny_se) <- list(sample="sample", group="group")

    tiny_normObj <- getVerifiedNormalyzerObject(
        jobName=jobName,
        summarizedExp=tiny_se,
        threshold=0,
        omitSamples=FALSE,
        requireReplicates=TRUE,
        quiet=TRUE,
        tinyRunThres=50
    )
    tiny_nr <- suppressWarnings(normMethods(tiny_normObj, normalizeRetentionTime=FALSE, quiet=TRUE))
    tiny_nr <- suppressWarnings(analyzeNormalizations(tiny_nr))

    outDir2 <- file.path(tempdir(), paste0("normplots2_", sample.int(1e9, 1)))
    dir.create(outDir2, recursive=TRUE)
    on.exit(unlink(outDir2, recursive=TRUE, force=TRUE), add=TRUE)

    pdfPath2 <- file.path(outDir2, paste0("Norm-report-", safeJobName, ".pdf"))
    pngDir2 <- file.path(outDir2, "pngs")

    expect_silent(suppressWarnings(generatePlots(tiny_nr, outDir2, writeAsPngs=FALSE)))
    expect_true(file.exists(pdfPath2))
    expect_false(dir.exists(pngDir2))
})

test_that("generateStatsReport handles missing PC3/PC4", {

    mat <- matrix(
        c(
            1, 2, 3, 4,
            4, 3, 2, 1
        ),
        nrow=2,
        byrow=TRUE,
        dimnames=list(c("f1", "f2"), c("s1", "s2", "s3", "s4"))
    )
    design <- data.frame(
        sample=colnames(mat),
        group=c("A", "A", "B", "B"),
        stringsAsFactors=FALSE
    )
    rownames(design) <- design$sample

    se <- SummarizedExperiment::SummarizedExperiment(
        assay=mat,
        colData=design,
        rowData=data.frame(feature=rownames(mat))
    )

    nst <- NormalyzerStatistics(se, logTrans=FALSE)
    nst <- calculateContrasts(
        nst,
        condCol="group",
        type="limma",
        leastRepCount=0,
        oneVsRest=TRUE,
        oneVsRestGroups="A"
    )

    jobName <- paste0("de_report_pca_", sample.int(1e9, 1))
    outDir <- file.path(tempdir(), paste0("report_pca_", sample.int(1e9, 1)))
    dir.create(outDir, recursive=TRUE)
    on.exit(unlink(outDir, recursive=TRUE, force=TRUE), add=TRUE)

    expect_silent(suppressWarnings(generateStatsReport(
        nst,
        jobName=jobName,
        jobDir=outDir,
        sigThres=0.5,
        sigThresType="fdr",
        log2FoldThres=0,
        writeAsPngs=TRUE
    )))
})
