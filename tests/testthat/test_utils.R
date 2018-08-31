context("utils.R")

test_that("getIndexList", {
    
    out <- getIndexList(c(1,2,3))
    expect_that(all.equal(out, list("1"=1, "2"=2, "3"=3)), is_true())
    
    out <- getIndexList(c(1,2,1,2,3))
    expect_that(all.equal(out, list("1"=c(1, 3), "2"=c(2, 4), "3"=5)), is_true())
})

test_that("getRowNAFilterContrast", {

    data("example_stat_data")
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

test_that("getReplicateSortedData_constant", {

    rawMat <- as.matrix(data.frame(
        "A"=c(1, 1, 1), 
        "A"=c(2, 2, 2),
        "B"=c(3, 3, 3),
        "B"=c(4, 4, 4)
    ))
    
    groups <- c("A", "A", "B", "B")
    sortedMat <- getReplicateSortedData(rawMat, groups)
            
    expect_that(
        all.equal(
            rawMat,
            sortedMat),
        is_true()
    )
})

test_that("getReplicateSortedData_reordering", {
    
    rawMat <- as.matrix(data.frame(
        "B"=c(3, 3, 3),
        "A"=c(1, 1, 1), 
        "B"=c(4, 4, 4),
        "A"=c(2, 2, 2)
    ))
    
    expectedMat <- as.matrix(data.frame(
        "A"=c(1, 1, 1), 
        "A"=c(2, 2, 2),
        "B"=c(3, 3, 3),
        "B"=c(4, 4, 4)
    ))
    
    groups <- c("B", "A", "B", "A")
    sortedMat <- getReplicateSortedData(rawMat, groups)
    
    expect_that(
        all.equal(
            sortedMat,
            expectedMat),
        is_true()
    )
})
