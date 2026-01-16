context("DIANN input")

test_that("setupRawDataObject reads DIANN pg_matrix and aligns sample names", {

    tmpDir <- tempdir()
    dataPath <- file.path(tmpDir, "diann_pg_matrix.tsv")
    designPath <- file.path(tmpDir, "diann_design.tsv")

    diannMatrix <- data.frame(
        `Protein.Group`=c("P1", "P2"),
        `Protein.Names`=c("Prot1", "Prot2"),
        Genes=c("G1", "G2"),
        `First.Protein.Description`=c("Desc1", "Desc2"),
        `N.Sequences`=c(1, 2),
        `N.Proteotypic.Sequences`=c(1, 2),
        `D:\\path\\S1.raw`=c(100, 0),
        `D:\\path\\S2.raw`=c(200, 300),
        check.names=FALSE
    )

    design <- data.frame(
        sample=c("S1", "S2"),
        group=c("A", "B"),
        stringsAsFactors=FALSE,
        check.names=FALSE
    )

    utils::write.table(diannMatrix, file=dataPath, sep="\t", row.names=FALSE, quote=FALSE)
    utils::write.table(design, file=designPath, sep="\t", row.names=FALSE, quote=FALSE)

    se <- setupRawDataObject(
        dataPath=dataPath,
        designPath=designPath,
        inputFormat="diann",
        sampleColName="sample",
        groupColName="group"
    )

    expect_true(identical(colnames(SummarizedExperiment::assay(se)), c("S1", "S2")))

    rowDf <- as.data.frame(SummarizedExperiment::rowData(se), optional=TRUE)
    expect_true(all(c("Protein.Group", "Protein.Names") %in% colnames(rowDf)))
})

test_that("setupRawDataObject reads DIANN report.tsv and aggregates duplicates", {

    tmpDir <- tempdir()
    dataPath <- file.path(tmpDir, "diann_report.tsv")
    designPath <- file.path(tmpDir, "diann_design_report.tsv")

    diannReport <- data.frame(
        `File.Name`=c("D:\\path\\S1.raw", "D:\\path\\S1.raw", "D:\\path\\S2.raw", "D:\\path\\S2.raw"),
        Run=c("S1", "S1", "S2", "S2"),
        `Protein.Group`=c("P1", "P1", "P1", "P2"),
        `Protein.Names`=c("Prot1", "Prot1", "Prot1", "Prot2"),
        Genes=c("G1", "G1", "G1", "G2"),
        `PG.Quantity`=c(100, 150, 200, 50),
        stringsAsFactors=FALSE,
        check.names=FALSE
    )

    design <- data.frame(
        sample=c("S1", "S2"),
        group=c("A", "B"),
        stringsAsFactors=FALSE,
        check.names=FALSE
    )

    utils::write.table(diannReport, file=dataPath, sep="\t", row.names=FALSE, quote=FALSE)
    utils::write.table(design, file=designPath, sep="\t", row.names=FALSE, quote=FALSE)

    se <- setupRawDataObject(
        dataPath=dataPath,
        designPath=designPath,
        inputFormat="diann",
        sampleColName="sample",
        groupColName="group"
    )

    dataMatrix <- SummarizedExperiment::assay(se)
    rowDf <- as.data.frame(SummarizedExperiment::rowData(se), optional=TRUE)

    p1Index <- which(rowDf$Protein.Group == "P1")
    p2Index <- which(rowDf$Protein.Group == "P2")

    expect_true(length(p1Index) == 1)
    expect_true(length(p2Index) == 1)

    expect_equal(as.numeric(dataMatrix[p1Index, "S1"]), 150)
    expect_equal(as.numeric(dataMatrix[p1Index, "S2"]), 200)
    expect_equal(as.numeric(dataMatrix[p2Index, "S1"]), NA_real_)
    expect_equal(as.numeric(dataMatrix[p2Index, "S2"]), 50)
})

test_that("DIANN report q-value filtering removes failing rows/features", {

    tmpDir <- tempdir()
    dataPath <- file.path(tmpDir, "diann_report_qfilter.tsv")
    designPath <- file.path(tmpDir, "diann_design_qfilter.tsv")

    diannReport <- data.frame(
        Run=c("S1", "S2", "S1", "S2"),
        `Precursor.Id`=c("pep1", "pep1", "pep2", "pep2"),
        `Protein.Group`=c("P1", "P1", "P2", "P2"),
        `Precursor.Quantity`=c(100, 200, 300, 400),
        `Q.Value`=c(0.005, 0.005, 0.02, 0.02),
        stringsAsFactors=FALSE,
        check.names=FALSE
    )

    design <- data.frame(
        sample=c("S1", "S2"),
        group=c("A", "B"),
        stringsAsFactors=FALSE,
        check.names=FALSE
    )

    utils::write.table(diannReport, file=dataPath, sep="\t", row.names=FALSE, quote=FALSE)
    utils::write.table(design, file=designPath, sep="\t", row.names=FALSE, quote=FALSE)

    se <- setupRawContrastObject(
        dataPath=dataPath,
        designPath=designPath,
        sampleColName="sample",
        inputFormat="diann",
        diannLevel="precursor",
        diannFilterQValue=TRUE,
        diannQValueCols=c("Q.Value"),
        diannQValueCutoffs=0.01
    )

    mat <- SummarizedExperiment::assay(se)
    ann <- as.data.frame(SummarizedExperiment::rowData(se), stringsAsFactors=FALSE)

    expect_true(nrow(mat) == 1)
    expect_true("Precursor.Id" %in% colnames(ann))
    expect_true(ann$Precursor.Id[1] == "pep1")
})

test_that("DIANN report precursor-level reading adds median RT annotation", {

    tmpDir <- tempdir()
    dataPath <- file.path(tmpDir, "diann_report_rt.tsv")
    designPath <- file.path(tmpDir, "diann_design_rt.tsv")

    diannReport <- data.frame(
        Run=c("S1", "S1", "S2", "S2"),
        `Precursor.Id`=c("pep1", "pep2", "pep1", "pep2"),
        `Protein.Group`=c("P1", "P2", "P1", "P2"),
        `Precursor.Quantity`=c(100, 200, 300, 400),
        RT=c(10, 20, 30, 40),
        `Q.Value`=c(0.001, 0.001, 0.001, 0.001),
        stringsAsFactors=FALSE,
        check.names=FALSE
    )

    design <- data.frame(
        sample=c("S1", "S2"),
        group=c("A", "B"),
        stringsAsFactors=FALSE,
        check.names=FALSE
    )

    utils::write.table(diannReport, file=dataPath, sep="\t", row.names=FALSE, quote=FALSE)
    utils::write.table(design, file=designPath, sep="\t", row.names=FALSE, quote=FALSE)

    se <- setupRawDataObject(
        dataPath=dataPath,
        designPath=designPath,
        inputFormat="diann",
        sampleColName="sample",
        groupColName="group",
        diannLevel="precursor",
        diannQuantityCol="Precursor.Quantity",
        diannFilterQValue=TRUE,
        diannQValueCols=c("Q.Value"),
        diannQValueCutoffs=0.01,
        diannRTCol="RT"
    )

    rowDf <- as.data.frame(SummarizedExperiment::rowData(se), optional=TRUE)
    expect_true("RT" %in% colnames(rowDf))

    pep1 <- rowDf[rowDf$Precursor.Id == "pep1", , drop=FALSE]
    pep2 <- rowDf[rowDf$Precursor.Id == "pep2", , drop=FALSE]
    expect_true(nrow(pep1) == 1)
    expect_true(nrow(pep2) == 1)

    expect_equal(as.numeric(pep1$RT), 20)
    expect_equal(as.numeric(pep2$RT), 30)
})

test_that("DIANN min-positive threshold converts tiny values to NA", {

    tmpDir <- tempdir()
    dataPath <- file.path(tmpDir, "diann_report_minpos.tsv")
    designPath <- file.path(tmpDir, "diann_design_minpos.tsv")

    diannReport <- data.frame(
        Run=c("S1", "S2"),
        `Precursor.Id`=c("pep1", "pep1"),
        `Protein.Group`=c("P1", "P1"),
        `Precursor.Quantity`=c(0.005, 0.02),
        `Q.Value`=c(0.001, 0.001),
        stringsAsFactors=FALSE,
        check.names=FALSE
    )

    design <- data.frame(
        sample=c("S1", "S2"),
        group=c("A", "B"),
        stringsAsFactors=FALSE,
        check.names=FALSE
    )

    utils::write.table(diannReport, file=dataPath, sep="\t", row.names=FALSE, quote=FALSE)
    utils::write.table(design, file=designPath, sep="\t", row.names=FALSE, quote=FALSE)

    se <- setupRawContrastObject(
        dataPath=dataPath,
        designPath=designPath,
        sampleColName="sample",
        inputFormat="diann",
        diannLevel="precursor",
        diannQuantityCol="Precursor.Quantity",
        diannMinPositive=0.01,
        diannFilterQValue=TRUE,
        diannQValueCols=c("Q.Value"),
        diannQValueCutoffs=0.01
    )

    mat <- SummarizedExperiment::assay(se)
    ann <- as.data.frame(SummarizedExperiment::rowData(se), stringsAsFactors=FALSE)

    pep1 <- which(ann$Precursor.Id == "pep1")
    expect_true(length(pep1) == 1)
    expect_true(is.na(as.numeric(mat[pep1, "S1"])))
    expect_equal(as.numeric(mat[pep1, "S2"]), 0.02)
})

test_that("DIANN level option selects precursor vs protein", {

    tmpDir <- tempdir()
    dataPath <- file.path(tmpDir, "diann_report_level.tsv")
    designPath <- file.path(tmpDir, "diann_design_level.tsv")

    diannReport <- data.frame(
        Run=c("S1", "S1", "S2", "S2"),
        `Protein.Group`=c("P1", "P1", "P1", "P1"),
        `Protein.Names`=c("Prot1", "Prot1", "Prot1", "Prot1"),
        Genes=c("G1", "G1", "G1", "G1"),
        `Precursor.Id`=c("pep1", "pep2", "pep1", "pep2"),
        `Precursor.Normalised`=c(10, 0, 20, 30),
        `PG.Quantity`=c(100, 150, 200, 50),
        stringsAsFactors=FALSE,
        check.names=FALSE
    )

    design <- data.frame(
        sample=c("S1", "S2"),
        group=c("A", "B"),
        stringsAsFactors=FALSE,
        check.names=FALSE
    )

    utils::write.table(diannReport, file=dataPath, sep="\t", row.names=FALSE, quote=FALSE)
    utils::write.table(design, file=designPath, sep="\t", row.names=FALSE, quote=FALSE)

    seProtein <- setupRawContrastObject(
        dataPath=dataPath,
        designPath=designPath,
        sampleColName="sample",
        inputFormat="diann",
        diannLevel="protein"
    )

    protMat <- SummarizedExperiment::assay(seProtein)
    protAnn <- as.data.frame(SummarizedExperiment::rowData(seProtein), stringsAsFactors=FALSE)
    expect_true(nrow(protMat) == 1)
    expect_true("Protein.Group" %in% colnames(protAnn))
    expect_equal(as.numeric(protMat[1, "S1"]), 150)
    expect_equal(as.numeric(protMat[1, "S2"]), 200)

    sePrec <- setupRawContrastObject(
        dataPath=dataPath,
        designPath=designPath,
        sampleColName="sample",
        inputFormat="diann",
        diannLevel="precursor"
    )

    precMat <- SummarizedExperiment::assay(sePrec)
    precAnn <- as.data.frame(SummarizedExperiment::rowData(sePrec), stringsAsFactors=FALSE)
    expect_true(nrow(precMat) == 2)
    expect_true("Precursor.Id" %in% colnames(precAnn))

    pep1 <- which(precAnn$Precursor.Id == "pep1")
    pep2 <- which(precAnn$Precursor.Id == "pep2")
    expect_true(length(pep1) == 1)
    expect_true(length(pep2) == 1)

    expect_equal(as.numeric(precMat[pep1, "S1"]), 10)
    expect_equal(as.numeric(precMat[pep1, "S2"]), 20)
    expect_true(is.na(as.numeric(precMat[pep2, "S1"])))
    expect_equal(as.numeric(precMat[pep2, "S2"]), 30)
})

test_that("setupRawContrastObject reads DIANN report.parquet", {
    testthat::skip_if_not_installed("arrow")

    tmpDir <- tempdir()
    dataPath <- file.path(tmpDir, "diann_report.parquet")
    designPath <- file.path(tmpDir, "diann_design_parquet.tsv")

    diannParquet <- data.frame(
        Run=c("S1", "S2", "S1"),
        `Protein.Group`=c("P1", "P1", "P2"),
        `Protein.Names`=c("Prot1", "Prot1", "Prot2"),
        Genes=c("G1", "G1", "G2"),
        `PG.MaxLFQ`=c(1000, 2000, 3000),
        stringsAsFactors=FALSE,
        check.names=FALSE
    )

    arrow::write_parquet(diannParquet, dataPath)

    design <- data.frame(
        sample=c("S1", "S2"),
        group=c("A", "B"),
        stringsAsFactors=FALSE,
        check.names=FALSE
    )
    utils::write.table(design, file=designPath, sep="\t", row.names=FALSE, quote=FALSE)

    se <- setupRawContrastObject(
        dataPath=dataPath,
        designPath=designPath,
        sampleColName="sample",
        inputFormat="diann"
    )

    dataMatrix <- SummarizedExperiment::assay(se)
    expect_true(is.matrix(dataMatrix))
    expect_true(is.numeric(dataMatrix))
    expect_true(identical(colnames(dataMatrix), c("S1", "S2")))
})
