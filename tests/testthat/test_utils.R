context("Evaluating util functions")

test_that("getIndexList", {
    
    out <- getIndexList(c(1,2,3))
    expect_that(all.equal(out, list("1"=1, "2"=2, "3"=3)), is_true())
    
    out <- getIndexList(c(1,2,1,2,3))
    expect_that(all.equal(out, list("1"=c(1, 3), "2"=c(2, 4), "3"=5)), is_true())
})

test_that("getRowNAFilterContrast", {

    data("example_data_only_values")
    data("example_design")
    levels <- example_design$group
    out <- getRowNAFilterContrast(head(example_stat_data), levels, minCount=3)
    expect_that(
        all.equal(
            out,
            c("1"=TRUE, "2"=TRUE, "3"=FALSE, "4"=FALSE, "5"=TRUE, "6"=TRUE)),
        is_true()
    )
})
