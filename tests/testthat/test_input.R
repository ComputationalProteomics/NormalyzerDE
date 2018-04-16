a <- 2
b <- 3

#print(getwd())
#path <- base::file.path(devtools::inst(name="normalyzer"), "testdata/data.csv")
# path <- system.file("testdata", package="normalyzer")
# print(path)
df <- read.csv("data.csv")
#df <- read.csv("../../inst/testdata/data.csv")
print(df)
#df <- read.csv(system.file("testdata/data.tsv", ..., package="normalyzer"))
#print(df)

test_that("true is true", {
    expect_equal(TRUE, TRUE)
})

test_that("two isn't three", {
    expect_equal(2 != 3, TRUE)
    expect_equal(a != b, TRUE)
})