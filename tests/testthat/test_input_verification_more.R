context("inputVerification (more)")

test_that("loadData supports MaxQuant input formats", {

    pepPath <- system.file("extdata", "mq_peptides_100.txt", package="NormalyzerDE")
    protPath <- system.file("extdata", "mq_proteinGroups_100.txt", package="NormalyzerDE")
    expect_true(nzchar(pepPath))
    expect_true(nzchar(protPath))

    pepMat <- loadData(pepPath, inputFormat="maxquantpep")
    expect_true(is.matrix(pepMat))
    expect_true(nrow(pepMat) > 1)
    expect_true(ncol(pepMat) > 1)

    protMat <- loadData(protPath, inputFormat="maxquantprot")
    expect_true(is.matrix(protMat))
    expect_true(nrow(protMat) > 1)
    expect_true(ncol(protMat) > 1)
})

test_that("loadData supports Proteios input format", {

    proteiosPath <- system.file("extdata", "tiny_data_proteios.tsv", package="NormalyzerDE")
    expect_true(nzchar(proteiosPath))

    protMat <- loadData(proteiosPath, inputFormat="proteios")
    expect_true(is.matrix(protMat))
    expect_true(nrow(protMat) > 0)
    expect_true(ncol(protMat) > 0)
})

test_that("loadData errors for unknown input formats", {
    expect_error(loadData("dummy", inputFormat="bad"), "Unknown inputFormat")
})

test_that("loadDesign errors when sample/group columns are missing", {

    fp <- tempfile(fileext=".tsv")
    on.exit(unlink(fp), add=TRUE)

    utils::write.table(
        data.frame(a=1, b=2),
        file=fp,
        sep="\t",
        row.names=FALSE,
        quote=FALSE
    )

    expect_error(
        loadDesign(fp, sampleCol="sample", groupCol="group"),
        "Both sampleCol value and groupCol value must be present"
    )
})

test_that("verifyValidNumbers errors for below-one values when log transforming", {

    mat <- matrix(c("2", "0.5"), nrow=1)
    expect_error(
        NormalyzerDE:::verifyValidNumbers(mat, groups=c("A", "B"), noLogTransform=FALSE, quiet=TRUE),
        "below-one values"
    )
})

test_that("getLowCountSampleFiltered can omit low-count samples", {

    mat <- matrix(c(1, 3, NA, 2, 4, NA), nrow=2, byrow=TRUE)
    colnames(mat) <- c("s1", "s2", "s3")
    out <- expect_warning(
        NormalyzerDE:::getLowCountSampleFiltered(mat, groups=c("A", "B", "C"), threshold=2, stopIfTooFew=FALSE),
        "does not contain enough"
    )

    expect_equal(colnames(out), c("s1", "s2"))
})

test_that("getLowCountSampleFiltered errors when all samples fail threshold", {

    mat <- matrix(NA_real_, nrow=2, ncol=2)
    colnames(mat) <- c("s1", "s2")
    expect_error(
        NormalyzerDE:::getLowCountSampleFiltered(mat, groups=c("A", "B"), threshold=1, stopIfTooFew=TRUE),
        "None of the samples had enough"
    )
})

