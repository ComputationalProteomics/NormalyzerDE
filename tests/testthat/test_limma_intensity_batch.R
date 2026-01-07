context("limma intensity + batch")

test_that("calculateContrasts supports type='limma_intensity' with batch column", {

    set.seed(1)
    mat <- matrix(stats::rnorm(50 * 6, mean=10, sd=1), nrow=50)
    colnames(mat) <- paste0("s", seq_len(ncol(mat)))

    design <- data.frame(
        sample=colnames(mat),
        group=c("A", "A", "A", "B", "B", "B"),
        batch=c("b1", "b2", "b1", "b2", "b1", "b2"),
        stringsAsFactors=FALSE
    )
    rownames(design) <- design$sample

    se <- SummarizedExperiment::SummarizedExperiment(
        assay=mat,
        colData=design,
        rowData=data.frame(feature=paste0("f", seq_len(nrow(mat))))
    )

    nst <- NormalyzerStatistics(se, logTrans=FALSE)
    out <- calculateContrasts(
        nst,
        comparisons="A-B",
        condCol="group",
        batchCol="batch",
        type="limma_intensity",
        leastRepCount=1
    )

    expect_true("A-B" %in% names(pairwiseCompsP(out)))
    expect_equal(length(pairwiseCompsP(out)[["A-B"]]), nrow(mat))
})

