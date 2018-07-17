context("Basic run including DE evaluation")

library(NormalyzerDE)

evalName <- "RegressionTest"
deName <- "RegressionTestDE"
outBase <- "output/BasicRunWithDE"
outEvalBase <- outBase
outDEBase <- paste(outBase, deName, sep="/")

designPath <- "data/basic_regression/test_design.tsv"
dataPath <- "data/basic_regression/test_data_w_val.tsv"

normalyzer(jobName=evalName,
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

# Read prepared MD5-sums into map
normalizationsMd5Fp <- "data/basic_regression/md5sums/normalizations.md5"
normalizationsMd5Df <- read.csv(normalizationsMd5Fp, sep="", col.names=c("md5", "name"), header=FALSE)
normalizationsMd5Map <- as.vector(normalizationsMd5Df[["md5"]])
names(normalizationsMd5Map) <- normalizationsMd5Df[["name"]]

statisticsMd5Fp <- "data/basic_regression/md5sums/statistics.md5"
statisticsMd5Df <- read.csv(statisticsMd5Fp, sep="", col.names=c("md5", "name"), header=FALSE)
statisticsMd5Map <- as.vector(statisticsMd5Df[["md5"]])
names(statisticsMd5Map) <- statisticsMd5Df[["name"]]

# Perform tests
test_that("Evaluation matrices are identical", {
    
    for (norm in normalizations) {
        
        baselineMd5 <- normalizationsMd5Map[[norm]]
        currentMd5 <- tools::md5sum(paste(outEvalBase, "RegressionTest", norm, sep="/"))
        expect_true(baselineMd5 == currentMd5, info = paste("Testing for normalization:", norm))
    }
    
})

for (normName in names(normalizations)) {

    norm <- normalizations[normName]
    
    allComps <- c("2-3", "2-4", "3-4")
    
    normalyzerDE(
        jobName=paste(deName, normName, sep="_"), 
        designPath=designPath, 
        dataPath=paste(outEvalBase, evalName, norm, sep="/"),
        outputDir=outDEBase,
        comparisons=allComps,
        type="limma",
        quiet=TRUE    
    )
    
    test_that("DE runs are identical", {
        
        currDeName <- paste(deName, normName, sep="_")
        currDeFileName <- paste0(currDeName, "_stats.tsv")
        currBase <- paste0(outDEBase, "/", currDeName)
        baselineMd5 <- statisticsMd5Map[[currDeFileName]]
        currentMd5 <- tools::md5sum(paste0(currBase, "/", currDeFileName))

        expect_true(baselineMd5 == currentMd5, info = paste("Testing for normalization:", norm))
    })
}



