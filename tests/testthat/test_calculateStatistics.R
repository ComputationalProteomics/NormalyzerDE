context("calculateStatistics.R")

data("example_data_only_values")

# Transfer to utils
test_that("filterLowRep", {

    test_data <- data.frame(
        c(NA, 1,  1,  NA), 
        c(1,  NA, 1,  NA), 
        c(3,  3,  3,  0), 
        c(5,  3,  NA, 0), 
        c(NA, 5,  NA, 0), 
        c(7,  9,  NA, 0))
    colnames(test_data) <- c("a1", "a2", "a3", "b1", "b2", "b3")
    groups <- c(rep("A", 3), rep("B", 3))
    
    expected_out_data <- data.frame(
        "a1"= c(NA, 1), 
        "a2"=c(1,  NA), 
        "a3"=c(3,  3), 
        "b1"=c(5,  3), 
        "b2"=c(NA, 5), 
        "b3"=c(7,  9))
    
    out <- filterLowRep(test_data, groups, leastRep=2)
    
    expect_that(
        all.equal(
            expected_out_data,
            out
        ),
        is_true()
    )
})

# test_that("reduceTechnicalReplicates", {
#     
#     tech_rep <- c("a", "a", "b", "b", "c", "c", "d", "d")
#     test_data <- data.frame(
#         c(1,1,1), 
#         c(1,2,1), 
#         c(3,3,3), 
#         c(5,3,3), 
#         c(5,5,4), 
#         c(5,5,5), 
#         c(7,7,7), 
#         c(7,9,7))
#     colnames(test_data) <- c("a1", "a2", "b1", "b2", "c1", "c2", "d1", "d2")
#     
#     expect_out_data <- as.matrix(data.frame(
#         "a"=c(1,1.5,1), 
#         "b"=c(4,3,3), 
#         "c"=c(5,5,4.5), 
#         "d"=c(7,8,7)))
#     
#     out <- reduceTechnicalReplicates(test_data, tech_rep)
#     
#     expect_that(
#         all.equal(
#             expect_out_data,
#             out
#         ),
#         is_true()
#     )
# })
# 
# test_that("reduceDesignTechRep", {
#     
#     test_df <- data.frame(
#         sample=c("a1", "a2", "a3", "b1", "b2", "c1", "c2", "d1"),
#         group=c(rep("A", 5), rep("B", 3)),
#         techrep=c("a", "a", "a", "b", "b", "c", "c", "d")
#     )
#     
#     expected_out_df <- data.frame(
#         sample=c("a1", "b1", "c1", "d1"),
#         group=c(rep("A", 2), rep("B", 2)),
#         techrep=c("a", "b", "c", "d")
#     )
#     
#     out <- reduceDesignTechRep(test_df, test_df$techrep)
#     
#     expect_that(
#         all.equal(
#             expected_out_df,
#             out
#         ),
#         is_true()
#     )
# })

test_that("reduceTechnicalReplicates", {
    
    test_data <- data.frame(
        c(1,1,1), 
        c(1,1.5,NA), 
        c(4,2,3), 
        c(3,2.5,3), 
        c(5,3.5,3), 
        c(4,3,4), 
        c(6,7,5), 
        c(7,9,NA))
    colnames(test_data) <- c("a1", "a2", "a3", "b1", "b2", "c1", "c2", "d1")
    
    test_df <- data.frame(
        sample=c("a1", "a2", "a3", "b1", "b2", "c1", "c2", "d1"),
        group=c(rep("A", 5), rep("B", 3)),
        techrep=c("a", "a", "a", "b", "b", "c", "c", "d")
    )

    se <- SummarizedExperiment::SummarizedExperiment(
        assay=as.matrix(test_data),
        colData=test_df,
        rowData=data.frame(annot=paste0("Pep", seq_len(3)))
    )
        
    expect_out_data <- as.matrix(data.frame(
        "a1.a2.a3"=c(2,1.5,2), 
        "b1.b2"=c(4,3,3), 
        "c1.c2"=c(5,5,4.5), 
        "d1"=c(7,9,NA)))
    
    expect_out_design <- data.frame(
        sample=c("a1.a2.a3", "b1.b2", "c1.c2", "d1"),
        group=c(rep("A", 2), rep("B", 2)),
        techrep=c("a", "b", "c", "d")
    )
    expect_out_design$sample <- as.character(expect_out_design$sample)
    rownames(expect_out_design) <- expect_out_design$sample
    
    expect_out_annot <- data.frame(
        annot=paste0("Pep", seq_len(3))
    )
    
    out_se <- reduceTechnicalReplicates(se, "techrep", "sample")
    
    expect_that(
        all.equal(
            expect_out_data,
            SummarizedExperiment::assay(out_se)
        ),
        is_true(),
        "Reduced data check"
    )
    
    expect_that(
        all.equal(
            expect_out_design,
            data.frame(SummarizedExperiment::colData(out_se))
        ),
        is_true(),
        "Reduced design check"
    )
    
    expect_that(
        all.equal(
            expect_out_annot,
            data.frame(SummarizedExperiment::rowData(out_se))
        ),
        is_true(),
        "Annotation remains identical"
    )
})



# Statistics
test_that("calculateWelch_data_test", {
    
    small_df <- log2(head(example_data_only_values))
    
    # small_df content
    # 
    # s_500amol_1 s_500amol_2 s_500amol_3 s_2500amol_1 s_2500amol_2 s_2500amol_3
    # [1,]    26.54878    26.52997    26.63352     26.89629     26.79226     26.74493
    # [2,]    26.62507    26.39669    26.75476     26.56381     26.36598     26.17868
    # [3,]    28.42273    28.45091    28.58193           NA     28.56563     28.66707
    # [4,]    25.53961    25.90364    25.38250     23.43634     24.55913     24.19715
    # [5,]    27.36890    27.29489    27.37707     26.60043     26.34719     26.22433
    # [6,]    28.90688    28.83352    28.85170     28.12027     28.31122     28.08823

    expected_p <- c(0.01484, 0.21877, 0.17375, 0.0271, 0.01001, 0.00592)
    expected_fdr <- c(0.02968, 0.21877, 0.20851, 0.04065, 0.02968, 0.02968)
    expected_fold <- c(-0.2404, 0.22268, -0.13116, 1.54437, 0.95631, 0.69079)
    expected_ave <- c(26.70042, 26.5699, 28.52404, 24.94016, 26.82584, 28.36693)
    
    header <- c(
        rep(1,3), rep(2,3), rep(3,3), 
        rep(4,3), rep(5,3), rep(6,3),
        rep(7,3), rep(8,3), rep(9,3))
    
    out <- calculateWelch(small_df, header, c(4,5))
    
    expect_that(all.equal(expected_p, round(out[["P"]], 5)), is_true())
    expect_that(all.equal(expected_fdr, round(out[["FDR"]], 5)), is_true())
    expect_that(all.equal(expected_fold, round(out[["Fold"]], 5)), is_true())
    expect_that(all.equal(expected_ave, round(out[["Ave"]], 5)), is_true())
})

test_that("calculateWelch_limited", {
    
    small_df <- data.frame(
        "a1" = c(1, 1, 1),
        "a2" = c(1.1, 1, 2),
        "a3" = c(0.9, 1, 3),
        "b1" = c(2, 1, 4),
        "b2" = c(2.1, 1, 5),
        "b3" = c(1.9, 1, 6)
    )
    
    expected_p <- c(0.00026, NA, 0.02131)
    expected_fdr <- c(0.00051, NA, 0.02131)
    expected_fold <- c(-1, 0, -3)
    expected_ave <- c(1.5, 1.0, 3.5)
    
    out <- calculateWelch(small_df, c(1,1,1,2,2,2), c(1,2))
    
    expect_that(all.equal(expected_p, round(out[["P"]], 5)), is_true())
    expect_that(all.equal(expected_fdr, round(out[["FDR"]], 5)), is_true())
    expect_that(all.equal(expected_fold, round(out[["Fold"]], 5)), is_true())
    expect_that(all.equal(expected_ave, round(out[["Ave"]], 5)), is_true())
})

test_that("calculateLimmaContrast_data_test", {
    
    small_df <- log2(head(example_data_only_values))[, seq(10, 15)]

    # small_df content
    # 
    # s_500amol_1 s_500amol_2 s_500amol_3 s_2500amol_1 s_2500amol_2 s_2500amol_3
    # [1,]    26.54878    26.52997    26.63352     26.89629     26.79226     26.74493
    # [2,]    26.62507    26.39669    26.75476     26.56381     26.36598     26.17868
    # [3,]    28.42273    28.45091    28.58193           NA     28.56563     28.66707
    # [4,]    25.53961    25.90364    25.38250     23.43634     24.55913     24.19715
    # [5,]    27.36890    27.29489    27.37707     26.60043     26.34719     26.22433
    # [6,]    28.90688    28.83352    28.85170     28.12027     28.31122     28.08823
    
    header <- c(rep(4, 3), rep(5, 3))
    levels <- c(4,5)
    
    Variable <- as.factor(header)
    model <- ~0+Variable
    limmaDesign <- stats::model.matrix(model)
    limmaFit <- limma::lmFit(small_df, limmaDesign)
    
    out <- calculateLimmaContrast(
        small_df, 
        limmaDesign, 
        limmaFit, 
        levels)
    
    expected_p <- c(0.0169, 0.14704, 0.2116, 0.00154, 8e-05, 0.00012)
    expected_fdr <- c(0.02535, 0.17645, 0.2116, 0.00308, 0.00037, 0.00037)
    expected_ave <- c(26.69096, 26.48083, 28.53765, 24.8364, 26.8688, 28.51864)
    expected_fold <- c(-0.2404, 0.22268, -0.13116, 1.54437, 0.95631, 0.69079)
    
    expect_that(all.equal(expected_p, round(out[["P"]], 5)), is_true())
    expect_that(all.equal(expected_fdr, round(out[["FDR"]], 5)), is_true())
    expect_that(all.equal(expected_fold, round(out[["Fold"]], 5)), is_true())
    expect_that(all.equal(expected_ave, round(out[["Ave"]], 5)), is_true())
})

test_that("calculateLimmaContrast_data_batch_test", {
    
    small_df <- log2(head(example_data_only_values))[, seq(10, 15)]
    
    # small_df content
    # 
    # s_500amol_1 s_500amol_2 s_500amol_3 s_2500amol_1 s_2500amol_2 s_2500amol_3
    # [1,]    26.54878    26.52997    26.63352     26.89629     26.79226     26.74493
    # [2,]    26.62507    26.39669    26.75476     26.56381     26.36598     26.17868
    # [3,]    28.42273    28.45091    28.58193           NA     28.56563     28.66707
    # [4,]    25.53961    25.90364    25.38250     23.43634     24.55913     24.19715
    # [5,]    27.36890    27.29489    27.37707     26.60043     26.34719     26.22433
    # [6,]    28.90688    28.83352    28.85170     28.12027     28.31122     28.08823
    
    header <- c(rep(4, 3), rep(5, 3))
    batch <- c(rep(c(1,2), 3))
    levels <- c(4,5)
    
    Variable <- as.factor(header)
    Batch <- as.factor(batch)
    model <- ~0+Variable+Batch
    limmaDesign <- stats::model.matrix(model)
    limmaFit <- limma::lmFit(small_df, limmaDesign)
    
    out <- calculateLimmaContrast(
        small_df, 
        limmaDesign, 
        limmaFit, 
        levels)
    
    expected_p <- c(0.04056, 0.28436, 0.28509, 0.00571, 0.00053, 0.00047)
    expected_fdr <- c(0.06083, 0.28509, 0.28509, 0.01142, 0.00158, 0.00158)
    expected_ave <- c(26.69096, 26.48083, 28.53765, 24.8364, 26.8688, 28.51864)
    expected_fold <- c(-0.24587, 0.17468, -0.12881, 1.49441, 0.95415, 0.64867)
    
    expect_that(all.equal(expected_p, round(out[["P"]], 5)), is_true())
    expect_that(all.equal(expected_fdr, round(out[["FDR"]], 5)), is_true())
    expect_that(all.equal(expected_fold, round(out[["Fold"]], 5)), is_true())
    expect_that(all.equal(expected_ave, round(out[["Ave"]], 5)), is_true())
})



