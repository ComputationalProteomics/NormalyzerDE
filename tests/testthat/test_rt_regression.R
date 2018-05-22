context("Global check of evaluation and DE by MD5sum comparisons to previous runs")

library(NormalyzerDE)

runName <- "TrickyRTRun"
runNameDE <- "TrickyRTRunDE"

outBase <- "output/TrickyRTRun"
# outEvalBase <- "TrickyRTRun"
# outDEBase <- "TrickyRTRun"

designPath <- "data/rt_regression/tricky_data_design.tsv"
dataPath <- "data/rt_regression/tricky_data_data_200.tsv"

baselineEvalDir <- "data/rt_regression/TrickyRTRun"
baselineDEDir <- "data/rt_regression/TrickyRTRunDE"
allComps <- c("8-9")

normalyzer(jobName=runName,
           designPath=designPath,
           dataPath=dataPath,
           outputDir=outBase,
           quiet=TRUE)

normalizations <- c(
    "RT-Loess"   = "RT-Loess-normalized.txt", 
    "RT-Mean"    = "RT-mean-normalized.txt", 
    "RT-Med"     = "RT-med-normalized.txt"
)

test_that("Evaluation matrices are identical", {
    
    for (norm in normalizations) {
        
        baseline_path <- paste(baselineEvalDir, norm, sep="/")
        eval_path <- paste(outBase, runName, norm, sep="/")
        
        baseline_md5 <- tools::md5sum(baseline_path)
        current_md5 <- tools::md5sum(eval_path)
        expect_true(baseline_md5 == current_md5, info = paste("Testing for normalization:", norm))
    }
})

for (norm_name in names(normalizations)) {

    norm <- normalizations[norm_name]
    
    normalyzerDE(
        jobName=paste(runNameDE, norm_name, sep="_"), 
        designPath=designPath, 
        dataPath=paste(baselineEvalDir, norm, sep="/"),
        outputDir=outBase,
        comparisons=allComps,
        type="limma",
        quiet=TRUE
    )
    
    test_that("DE runs are identical", {
        
        baseline_base <- paste0(baselineDEDir, "/", norm_name)
        curr_base <- paste0(outBase, "/", runNameDE, "_", norm_name)
        
        baseline_md5 <- tools::md5sum(paste0(baseline_base, "/", norm_name, "_stats.tsv"))
        current_md5 <- tools::md5sum(paste0(curr_base, "/", runNameDE, "_", norm_name, "_stats.tsv"))

        expect_true(baseline_md5 == current_md5, info = paste("Testing for normalization:", norm))
    })
}



