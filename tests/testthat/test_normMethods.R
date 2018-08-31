context("normMethods.R")

data("example_design")
data("example_data_only_values")

# Subset the data to only look at first three conditions
test_design <- example_design[example_design$group %in% c("1", "2", "3"), ]
test_data <- example_data_only_values[, as.character(test_design$sample)]

# Remove rows with only NA-values
test_data <- test_data[rowSums(is.na(test_data)) != ncol(test_data), ]

test_that("globalIntensityNormalization", {
    
    expect_dim <- c(98, 9)
    expect_colsum <- c(
        2481.86278, 2420.0029, 2493.85063, 
        2464.60255, 2422.10488, 2473.12488, 
        2516.12334, 2432.26214, 2390.14505
    )
    
    out <- globalIntensityNormalization(test_data)
    out_dim <- dim(out)
    out_colsum <- round(colSums(out, na.rm=TRUE), 5)
    
    expect_that(
        all.equal(
            expect_dim,
            out_dim
        ),
        is_true()
    )
    
    expect_that(
        all(
            expect_colsum == out_colsum
        ),
        is_true()
    )
})

test_that("medianNormalization", {
    
    expect_dim <- c(98, 9)
    expect_colsum <- c(
        2476.44406, 2416.74862, 2514.76548, 
        2480.75296, 2421.68804, 2460.60737, 
        2515.09438, 2403.50676, 2382.08166
    )
    
    out <- medianNormalization(test_data)
    out_dim <- dim(out)
    out_colsum <- round(colSums(out, na.rm=TRUE), 5)
    
    expect_that(
        all.equal(
            expect_dim,
            out_dim
        ),
        is_true()
    )
    
    expect_that(
        all(
            expect_colsum == out_colsum
        ),
        is_true()
    )
})

test_that("meanNormalization", {
    
    expect_dim <- c(98, 9)
    expect_colsum <- c(
        2480.02877, 2415.35575, 2493.44705, 
        2462.76854, 2417.45773, 2471.29087, 
        2517.16589, 2427.61499, 2384.11555
    )
    
    out <- meanNormalization(test_data)
    out_dim <- dim(out)
    out_colsum <- round(colSums(out, na.rm=TRUE), 5)
    
    expect_that(
        all.equal(
            expect_dim,
            out_dim
        ),
        is_true()
    )
    
    expect_that(
        all(
            expect_colsum == out_colsum
        ),
        is_true()
    )
})

test_that("performVSNNormalization", {
    
    expect_dim <- c(98, 9)
    expect_colsum <- c(
        2470.56447, 2407.87325, 2496.00488, 
        2461.14702, 2421.64088, 2477.71483, 
        2520.74195, 2413.61712, 2397.58554
    )
    
    out <- performVSNNormalization(test_data)
    out_dim <- dim(out)
    out_colsum <- round(colSums(out, na.rm=TRUE), 5)
    
    expect_that(
        all.equal(
            expect_dim,
            out_dim
        ),
        is_true()
    )
    
    expect_that(
        all(
            expect_colsum == out_colsum
        ),
        is_true()
    )
})

test_that("performQuantileNormalization", {
    
    expect_dim <- c(98, 9)
    expect_colsum <- c(
        2469.98486, 2415.69332, 2497.12467, 
        2469.98486, 2415.69332, 2469.98486, 
        2524.27466, 2415.69332, 2388.56021
    )
    
    out <- performQuantileNormalization(test_data)
    out_dim <- dim(out)
    out_colsum <- round(colSums(out, na.rm=TRUE), 5)
    
    expect_that(
        all.equal(
            expect_dim,
            out_dim
        ),
        is_true()
    )
    
    expect_that(
        all(
            expect_colsum == out_colsum
        ),
        is_true()
    )
})

test_that("performSMADNormalization", {
    
    expect_dim <- c(98, 9)
    expect_colsum <- c(
        2474.668, 2417.49851, 2507.3658, 
        2476.69406, 2419.97969, 2466.17639, 
        2520.90251, 2411.10946, 2386.40614
    )
    
    out <- performSMADNormalization(test_data)
    out_dim <- dim(out)
    out_colsum <- round(colSums(out, na.rm=TRUE), 5)
    
    expect_that(
        all.equal(
            expect_dim,
            out_dim
        ),
        is_true()
    )
    
    expect_that(
        all(
            expect_colsum == out_colsum
        ),
        is_true()
    )
})

test_that("performCyclicLoessNormalization", {
    
    expect_dim <- c(98, 9)
    expect_colsum <- c(
        2478.31304, 2410.73689, 2494.80824, 
        2454.38212, 2418.7555, 2471.09178, 
        2510.40413, 2402.29372, 2389.24695
    )
    
    out <- performCyclicLoessNormalization(test_data)
    out_dim <- dim(out)
    out_colsum <- round(colSums(out, na.rm=TRUE), 5)
    
    expect_that(
        all.equal(
            expect_dim,
            out_dim
        ),
        is_true()
    )
    
    expect_that(
        all(
            expect_colsum == out_colsum
        ),
        is_true()
    )
})

test_that("performGlobalRLRNormalization", {
    
    expect_dim <- c(98, 9)
    expect_colsum <- c(
        2481.46193, 2416.7949, 2499.26657, 
        2461.07399, 2424.5646, 2474.91485, 
        2517.3423, 2408.12205, 2393.68206
    )
    
    out <- performGlobalRLRNormalization(test_data)
    out_dim <- dim(out)
    out_colsum <- round(colSums(out, na.rm=TRUE), 5)
    
    expect_that(
        all.equal(
            expect_dim,
            out_dim
        ),
        is_true()
    )
    
    expect_that(
        all(
            expect_colsum == out_colsum
        ),
        is_true()
   )
})

test_that("performNoNormalization", {
    
    expect_dim <- c(98, 9)
    expect_colsum <- c(
        35539172278, 35144510958, 36832468728, 
        36954497268, 38286168846, 38114456004, 
        37366276438, 30257069174, 37091357660
    )
    
    out <- performNoNormalization(test_data)
    out_dim <- dim(out)
    out_colsum <- round(colSums(out, na.rm=TRUE), 5)
    
    expect_that(
        all.equal(
            expect_dim,
            out_dim
        ),
        is_true()
    )
    
    expect_that(
        all(
            expect_colsum == out_colsum
        ),
        is_true()
    )
})
