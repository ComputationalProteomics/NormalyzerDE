context("Global check of evaluation and DE by MD5sum comparisons to previous runs")

library(NormalyzerDE)

outBase <- "output/BasicRunWithDE"
outEvalBase <- outBase
outDEBase <- paste(outBase, "RegressionTestDE", sep="/")

designPath <- "data/regression_cases/web_dataset/test_design.tsv"
dataPath <- "data/regression_cases/web_dataset/test_data.tsv"

baselineEvalDir <- "data/regression_cases/web_dataset/output/eval"
baselineDEDir <- "data/regression_cases/web_dataset/output/de"

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
    
    # baseline_base <- "data/regression_cases/web_dataset/output/eval"
    # curr_base <- "RegressionTest/"
    
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
        jobName=paste("RegressionTestDE", norm_name, sep="_"), 
        designPath=designPath, 
        dataPath=paste(baselineEvalDir, norm, sep="/"),
        outputDir=outDEBase,
        comparisons=allComps,
        type="limma",
        quiet=TRUE
    )
    
    test_that("DE runs are identical", {
        
        baseline_base <- paste0(baselineDEDir, "/", norm_name, "/", norm_name)
        curr_base <- paste0(outDEBase, "/", "RegressionTestDE_", norm_name)
        
        baseline_md5 <- tools::md5sum(paste0(baseline_base, "/", norm_name, "_stats.tsv"))
        current_md5 <- tools::md5sum(paste0(curr_base, "/", "RegressionTestDE_", norm_name, "_stats.tsv"))

        expect_true(baseline_md5 == current_md5, info = paste("Testing for normalization:", norm))
    })
}



