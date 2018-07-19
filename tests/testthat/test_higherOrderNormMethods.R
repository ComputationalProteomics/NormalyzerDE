context("higherOrderNormMethods.R")

data("example_design")
data("example_data_only_values")
data("example_data")

# Subset the data to only look at first three conditions
test_design <- example_design[which(example_design$group %in% c("1", "2", "3")), ]
test_data <- example_data_only_values[, as.character(test_design$sample)]

# Remove rows with only NA-values
test_data <- test_data[rowSums(is.na(test_data)) != ncol(test_data), ]

test_that("getRTNormalizedMatrix", {})

test_that("getWidenedRTRange_natural_numbers", {
    
    rt_times <- seq_len(10)
    
    expect_that(
        all.equal(
            getWidenedRTRange(4, 5, 2, rt_times),
            c(4, 5)
        ),
        is_true()
    )
    
    expect_that(
        all.equal(
            getWidenedRTRange(4, 5, 3, rt_times),
            c(3, 5)
        ),
        is_true()
    )
    
    expect_that(
        all.equal(
            getWidenedRTRange(4, 5, 4, rt_times),
            c(3, 6)
        ),
        is_true()
    )
})

test_that("getWidenedRTRange_float_numbers", {
    
    rt_times <- c(1.1, 2.2, 3.3, 4.2, 5.4, 6.7, 7.1, 8.2, 9.2, 10.3)
    
    expect_that(
        all.equal(
            getWidenedRTRange(2.2, 3.3, 2, rt_times),
            c(2.2, 3.3)
        ),
        is_true()
    )
    
    expect_that(
        all.equal(
            getWidenedRTRange(2.2, 3.3, 6, rt_times),
            c(1.1, 6.7)
        ),
        is_true()
    )
    
    expect_that(
        all.equal(
            getWidenedRTRange(2.2, 3.3, 10, rt_times),
            c(1.1, 10.3)
        ),
        is_true()
    )
    
    expect_error(
        getWidenedRTRange(2.2, 3.3, 11, rt_times)
    )
})

test_that("getWidenedRTRange_real_numbers", {
    
    real_rt_vals <- round(example_stat_data$Average.RT, 5)
    
    expect_that(
        all.equal(
            getWidenedRTRange(34.49672, 45.08198, 20, real_rt_vals),
            c(32.61380, 45.33467)
        ),
        is_true()
    )
    
    expect_that(
        all.equal(
            getWidenedRTRange(34.49672, 45.08198, 30, real_rt_vals),
            c(28.60213, 48.11109)
        ),
        is_true()
    )
    
    expect_that(
        all.equal(
            getWidenedRTRange(34.49672, 45.08198, 100, real_rt_vals),
            c(9.02435, 128.70183)
        ),
        is_true()
    )
    
    expect_error(
        getWidenedRTRange(34.49672, 45.08198, 101, real_rt_vals)
    )
    
    expect_error(
        getWidenedRTRange(34.49672, 45.08198, 10, real_rt_vals)
    )
    
    expect_that(
        all.equal(
            getWidenedRTRange(
                34.49672, 45.08198, 10, real_rt_vals, allowTooWideData=TRUE),
            c(34.49672, 45.08198)
        ),
        is_true()
    )
})

test_that("getSmoothedRTNormalizedMatrix", {})

test_that("getCombinedMatrix_minimal_symmetric", {
    
    expect_out <- as.matrix(data.frame(
        a=c(1.5, 1.5), 
        b=c(1.5, 1.5)))
    
    df1 <- as.matrix(data.frame(a=c(1, 1), b=c(1, 1)))
    df2 <- as.matrix(data.frame(a=c(2, 2), b=c(2, 2)))
    l <- list(df1, df2)
    out <- getCombinedMatrix(l, mean)
    
    expect_that(
        all.equal(
            expect_out,
            out
        ),
        is_true()
    )
})

test_that("getCombinedMatrix_minimal_asymmetric", {
    
    expect_out <- as.matrix(data.frame(
        a=c(0.5, 1, 1.5), 
        b=c(2, 2.5, 3)))
    
    df1 <- as.matrix(data.frame(a=c(0, 0, 0), b=c(0, 0, 0)))
    df2 <- as.matrix(data.frame(a=c(1, 2, 3), b=c(4, 5, 6)))
    l <- list(df1, df2)
    out <- getCombinedMatrix(l, mean)
    
    expect_that(
        all.equal(
            expect_out,
            out
        ),
        is_true()
    )
})

test_that("getCombinedMatrix_three_layer_mean", {
    
    expect_out <- as.matrix(data.frame(
        a=c(1.00000, 3.33333, 5.66667), 
        b=c(8.00000, 10.33333, 12.66667)))
    
    df1 <- as.matrix(data.frame(a=c(1, 2, 3), b=c(4, 5, 6)))
    df2 <- as.matrix(data.frame(a=c(1, 3, 5), b=c(7, 9, 11)))
    df3 <- as.matrix(data.frame(a=c(1, 5, 9), b=c(13, 17, 21)))
    
    l <- list()
    l[[1]] <- df1
    l[[2]] <- df2
    l[[3]] <- df3
    
    out <- round(getCombinedMatrix(l, mean), 5)
    
    expect_that(
        all.equal(
            expect_out,
            out
        ),
        is_true()
    )
})

test_that("getCombinedMatrix_three_layer_median", {
    
    expect_out <- as.matrix(data.frame(a=c(1, 3, 5), b=c(7, 9, 11)))
    
    df1 <- as.matrix(data.frame(a=c(1, 2, 3), b=c(4, 5, 6)))
    df2 <- as.matrix(data.frame(a=c(1, 3, 5), b=c(7, 9, 11)))
    df3 <- as.matrix(data.frame(a=c(1, 5, 9), b=c(13, 17, 21)))
    
    l <- list()
    l[[1]] <- df1
    l[[2]] <- df2
    l[[3]] <- df3
    
    out <- getCombinedMatrix(l, median)
    
    expect_that(
        all.equal(
            expect_out,
            out
        ),
        is_true()
    )
})

