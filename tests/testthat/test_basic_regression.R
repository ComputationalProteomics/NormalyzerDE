context("Basic run including DE evaluation")

library(NormalyzerDE)

deName <- "RegressionTestDE"
outBase <- "output/BasicRunWithDE"
outEvalBase <- outBase
outDEBase <- paste(outBase, deName, sep="/")

designPath <- "data/basic_regression/test_design.tsv"
dataPath <- "data/basic_regression/test_data_w_val.tsv"

baselineEvalDir <- "data/basic_regression/BasicRunWithDE/RegressionTest"
baselineDEDir <- "data/basic_regression/BasicRunWithDE/RegressionTestDE"

normalyzer(jobName="RegressionTest",
           designPath=designPath,
           dataPath=dataPath,
           outputDir=outEvalBase,
           quiet=TRUE)

normalizations <- c(
    "AI-G"     = "AI-G-normalized.txt", 
    "Loess-G"  = "Loess-G-normalized.txt", 
    "Log-2"    = "Log2-normalized.txt", 
    "MeanI-G"  = "MeanI-G-normalized.txt", 
    "MedI-G"   = "MedI-G-normalized.txt", 
    "Quantile" = "Quantile-normalized.txt", 
    "submitted_rawdata" = "submitted_rawdata.txt", 
    "VSN-G"    = "VSN-G-normalized.txt"
)

test_that("Evaluation matrices are identical", {
    
    for (norm in normalizations) {
        
        baseline_md5 <- tools::md5sum(paste(baselineEvalDir, norm, sep="/"))
        current_md5 <- tools::md5sum(paste(outEvalBase, "RegressionTest", norm, sep="/"))
        expect_true(baseline_md5 == current_md5, info = paste("Testing for normalization:", norm))
    }
    
})

for (norm_name in names(normalizations)) {

    norm <- normalizations[norm_name]
    
    allComps <- c("2-3", "2-4", "3-4")
    
    normalyzerDE(
        jobName=paste(deName, norm_name, sep="_"), 
        designPath=designPath, 
        dataPath=paste(baselineEvalDir, norm, sep="/"),
        outputDir=outDEBase,
        comparisons=allComps,
        type="limma",
        quiet=TRUE    
    )
    
    test_that("DE runs are identical", {
        
        currDeName <- paste(deName, norm_name, sep="_")
        
        baseline_base <- paste0(baselineDEDir, "/", currDeName)
        curr_base <- paste0(outDEBase, "/", currDeName)
        
        baseline_md5 <- tools::md5sum(paste0(baseline_base, "/", currDeName, "_stats.tsv"))
        current_md5 <- tools::md5sum(paste0(curr_base, "/", currDeName, "_stats.tsv"))

        expect_true(baseline_md5 == current_md5, info = paste("Testing for normalization:", norm))
    })
}



