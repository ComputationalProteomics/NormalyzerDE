context("test testing")

test_that("Dummy function returns one", {
    value <- returnOne()
    expect_that(value, equals(1))
})