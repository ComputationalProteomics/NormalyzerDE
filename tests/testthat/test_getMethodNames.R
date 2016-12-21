context("Tests")
test_that("getMethodNames function returns given set of method names without 
          provided flag", {
    expectedMethodNames <- c("Log2", "TI-G", "MedI-G", "AI-G")
    houseKeepingFlag <- TRUE
    
    receivedMethodNames <- getMethodNames(houseKeepingFlag)
    
    for(i in 1:length(expectedMethodNames)) {
        expect_equal(expectedMethodNames[i], receivedMethodNames[i])
    }
})

test_that("Dummy succeed test", {
    expect_equal(1, 1)
})