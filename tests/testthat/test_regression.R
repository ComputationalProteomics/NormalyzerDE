library(NormalyzerDE)


normalyzer("data/regression_cases/web_dataset/test_data.tsv")

test_that("true is true", {
    expect_equal(TRUE, TRUE)
})

test_that("two isn't three", {
    expect_equal(2 != 3, TRUE)
    expect_equal(a != b, TRUE)
})

all_NA_one_replicate_set_m <- as.matrix(read.table("data/corner_cases/all_NA_one_replicate_set.tsv", header=FALSE, sep="\t", stringsAsFactors=FALSE, quote="", comment.char=""))
test_that("Raw loading for standard file works", {
    expect_equal(all_NA_one_replicate_set_m, loadRawDataFromFile("data/corner_cases/all_NA_one_replicate_set.tsv"))
})

