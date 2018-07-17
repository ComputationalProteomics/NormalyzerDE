context("RT run evaluation")

library(NormalyzerDE)

runName <- "TrickyRTRun"
runNameDE <- "TrickyRTRunDE"

outBase <- "output/TrickyRTRun"

designPath <- "data/rt_regression/tricky_data_design.tsv"
dataPath <- "data/rt_regression/tricky_data_data_200.tsv"

allComps <- c("8-9")

normalyzer(jobName=runName,
           designPath=designPath,
           dataPath=dataPath,
           outputDir=outBase,
           quiet=TRUE)

normalizations <- c(
    "RT-Loess"   = "RT-Loess-normalized.txt", 
    "RT-mean"    = "RT-mean-normalized.txt", 
    "RT-med"     = "RT-med-normalized.txt"
)

# Read prepared MD5-sums into map
normalizationsMd5Fp <- "data/rt_regression/md5sums/normalizations.md5"
normalizationsMd5Df <- read.csv(normalizationsMd5Fp, sep="", col.names=c("md5", "name"), header=FALSE)
normalizationsMd5Map <- as.vector(normalizationsMd5Df[["md5"]])
names(normalizationsMd5Map) <- normalizationsMd5Df[["name"]]

statisticsMd5Fp <- "data/rt_regression/md5sums/statistics.md5"
statisticsMd5Df <- read.csv(statisticsMd5Fp, sep="", col.names=c("md5", "name"), header=FALSE)
statisticsMd5Map <- as.vector(statisticsMd5Df[["md5"]])
names(statisticsMd5Map) <- statisticsMd5Df[["name"]]

# Perform tests
test_that("Evaluation matrices are identical", {
    
    for (normName in normalizations) {
        
        baselineMd5 <- normalizationsMd5Map[[normName]]
        evalPath <- paste(outBase, runName, normName, sep="/")
        currentMd5 <- tools::md5sum(evalPath)
        expect_true(baselineMd5 == currentMd5, info = paste("Testing for normalization:", normName))
    }
})

for (normNameRaw in names(normalizations)) {

    normName <- normalizations[normNameRaw]
    
    normalyzerDE(
        jobName=paste(runNameDE, normName, sep="_"), 
        designPath=designPath, 
        dataPath=paste(outBase, runName, normName, sep="/"),
        outputDir=outBase,
        comparisons=allComps,
        type="limma",
        quiet=TRUE,
        batchCol="batch"
    )
    
    test_that("DE runs are identical", {
        
        currBase <- paste0(outBase, "/", runNameDE, "_", normName)
        baselineMd5 <- statisticsMd5Map[[paste0(normNameRaw, "_stats.tsv")]]
        currentMd5 <- tools::md5sum(paste0(currBase, "/", runNameDE, "_", normName, "_stats.tsv"))
        # print(paste0(currBase,  "/", runNameDE, "_", normName))
        # print(paste("MD5s", baselineMd5, currentMd5))
        
        expect_true(baselineMd5 == currentMd5, info = paste("Testing for normalization:", normName))
    })
}



