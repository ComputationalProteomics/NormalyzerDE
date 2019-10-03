context("higherOrderNormMethods.R")

data("example_design")
data("example_data_only_values")
data("example_data")
data("example_stat_data")

# Subset the data to only look at first three conditions
test_design <- example_design[example_design$group %in% c("1", "2", "3"), ]
test_data <- example_data_only_values[, as.character(test_design$sample)]

# Remove rows with only NA-values
non_na_rows <- which(rowSums(is.na(test_data)) != ncol(test_data))
test_data <- test_data[non_na_rows, ]
real_rt_vals <- round(example_stat_data$Average.RT, 5)[non_na_rows]


test_that("getRTNormalizedMatrix_global_intensity", {
    
    expect_dim <- c(98, 9)
    expect_colsum <- c(
        2482.77553, 2420.75378, 2497.33571, 
        2464.09042, 2424.41947, 2472.58042, 
        2515.43627, 2419.77311, 2393.76557
    )
    
    out <- getRTNormalizedMatrix(
        rawMatrix=as.matrix(test_data), 
        retentionTimes=real_rt_vals, 
        normMethod=globalIntensityNormalization, 
        stepSizeMinutes=1, 
        windowMinCount=10)
    
    out_dim <- dim(out)
    out_colsum <- round(colSums(out, na.rm=TRUE), 5)
    
    expect_true(
        all.equal(
            expect_dim,
            out_dim
        )
    )
    
    expect_true(
        all(
            expect_colsum == out_colsum
        )
    )
})

test_that("getRTNormalizedMatrix_loess", {
    
    expect_dim <- c(98, 9)
    expect_colsum <- c(
        2472.8547, 2412.4723, 2496.2553, 
        2458.90508, 2419.86546, 2471.89521, 
        2515.73799, 2413.87385, 2392.8269
    )
    
    out <- getRTNormalizedMatrix(
        rawMatrix=as.matrix(test_data), 
        retentionTimes=real_rt_vals, 
        normMethod=performCyclicLoessNormalization, 
        stepSizeMinutes=1, 
        windowMinCount=10)
    
    out_dim <- dim(out)
    out_colsum <- round(colSums(out, na.rm=TRUE), 5)
    
    expect_true(
        all.equal(
            expect_dim,
            out_dim
        )
    )
    
    expect_true(
        all(
            expect_colsum == out_colsum
        )
    )
})

test_that("getSmoothedRTNormalizedMatrix_global_intensity", {
    
    expect_dim <- c(98, 9)
    expect_colsum <- c(
        2483.01579, 2420.61195, 2497.10932, 
        2463.82244, 2424.36461, 2472.49291, 
        2515.79797, 2421.82197, 2393.95915
    )
    
    out <- getSmoothedRTNormalizedMatrix(
        rawMatrix=as.matrix(test_data), 
        retentionTimes=real_rt_vals, 
        normMethod=globalIntensityNormalization, 
        stepSizeMinutes=1, 
        windowMinCount=10,
        mergeMethod="median",
        windowShifts=4
        )
    
    out_dim <- dim(out)
    out_colsum <- round(colSums(out, na.rm=TRUE), 5)
    
    expect_true(
        all.equal(
            expect_dim,
            out_dim
        )
    )
    
    expect_true(
        all(
            expect_colsum == out_colsum
        )
    )
})

test_that("getSmoothedRTNormalizedMatrix_loess", {
    
    expect_dim <- c(98, 9)
    expect_colsum <- c(
        2473.03148, 2412.86092, 2496.37835, 
        2458.6256, 2419.13168, 2471.7457, 
        2514.70862, 2413.37041, 2392.00512
    )
    
    out <- getSmoothedRTNormalizedMatrix(
        rawMatrix=as.matrix(test_data), 
        retentionTimes=real_rt_vals, 
        normMethod=performCyclicLoessNormalization, 
        stepSizeMinutes=1, 
        windowMinCount=10,
        mergeMethod="mean",
        windowShifts=3
        )
    
    out_dim <- dim(out)
    out_colsum <- round(colSums(out, na.rm=TRUE), 5)
    
    expect_true(
        all.equal(
            expect_dim,
            out_dim
        )
    )
    
    expect_true(
        all(
            expect_colsum == out_colsum
        )
    )
})

test_that("getWidenedRTRange_natural_numbers", {
    
    rt_times <- seq_len(10)
    
    expect_true(
        all.equal(
            getWidenedRTRange(4, 5, 2, rt_times),
            c(4, 5)
        )
    )
    
    expect_true(
        all.equal(
            getWidenedRTRange(4, 5, 3, rt_times),
            c(3, 5)
        )
    )
    
    expect_true(
        all.equal(
            getWidenedRTRange(4, 5, 4, rt_times),
            c(3, 6)
        )
    )
})

test_that("getWidenedRTRange_float_numbers", {
    
    rt_times <- c(1.1, 2.2, 3.3, 4.2, 5.4, 6.7, 7.1, 8.2, 9.2, 10.3)
    
    expect_true(
        all.equal(
            getWidenedRTRange(2.2, 3.3, 2, rt_times),
            c(2.2, 3.3)
        )
    )
    
    expect_true(
        all.equal(
            getWidenedRTRange(2.2, 3.3, 6, rt_times),
            c(1.1, 6.7)
        )
    )
    
    expect_true(
        all.equal(
            getWidenedRTRange(2.2, 3.3, 10, rt_times),
            c(1.1, 10.3)
        )
    )
    
    expect_error(
        getWidenedRTRange(2.2, 3.3, 11, rt_times)
    )
})

test_that("getWidenedRTRange_real_numbers", {
    
    expect_true(
        all.equal(
            getWidenedRTRange(34.49672, 45.08198, 20, real_rt_vals),
            c(30.93696, 45.33467)
        )
    )
    
    expect_true(
        all.equal(
            getWidenedRTRange(34.49672, 45.08198, 30, real_rt_vals),
            c(28.21106, 48.11109)
        )
    )
    
    expect_true(
        all.equal(
            getWidenedRTRange(34.49672, 45.08198, 98, real_rt_vals),
            c(9.02435, 128.70183)
        )
    )
    
    expect_error(
        getWidenedRTRange(34.49672, 45.08198, 100, real_rt_vals)
    )
    
    expect_error(
        getWidenedRTRange(34.49672, 45.08198, 10, real_rt_vals)
    )
    
    expect_true(
        all.equal(
            getWidenedRTRange(
                34.49672, 45.08198, 10, real_rt_vals, allowTooWideData=TRUE),
            c(34.49672, 45.08198)
        )
    )
})

test_that("getCombinedMatrix_minimal_symmetric", {
    
    expect_out <- as.matrix(data.frame(
        a=c(1.5, 1.5), 
        b=c(1.5, 1.5)))
    
    df1 <- as.matrix(data.frame(a=c(1, 1), b=c(1, 1)))
    df2 <- as.matrix(data.frame(a=c(2, 2), b=c(2, 2)))
    l <- list(df1, df2)
    out <- getCombinedMatrix(l, mean)
    
    expect_true(
        all.equal(
            expect_out,
            out
        )
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
    
    expect_true(
        all.equal(
            expect_out,
            out
        )
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
    
    expect_true(
        all.equal(
            expect_out,
            out
        )
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
    
    expect_true(
        all.equal(
            expect_out,
            out
        )
    )
})

