print('Normalyzer sourced')

normalyzer <- function(datafile, jobName, outputDir=NULL) {
    
    print('start')
    
    require(Rcmdr)
    require(PerformanceAnalytics)
    require(vsn)
    require(preprocessCore)
    require(limma)
    require(MASS)
    require(abind)
    require(e1071)
    require(ape)
    require(raster)
    require(car)
    require(gridExtra)
    require(ggplot2)

    require(grid)
    
    source("generatePlots.R")
    source("normfinder-pipeline.R")
    source("normMethods.R")
    source("printMeta.R")
    source("printPlots.R")

    source("NormalyzerDataset.R")
    source("NormalyzerResults.R")
    source("NormalizationEvaluationResults.R")
    
    source("utils.R")
    source("inputVerification.R")
    source("analyzeResults.R")
    
    normObj <- getVerifiedNormalyzerObjectFromFile(datafile, jobName)
    jobDir <- setupJobDir(jobName, outputDir)
        
    print("Normalizing data....")
    normalyzerResultsObject <- normMethods(normObj, jobName, jobDir)
    print("Finished Normalization")
    
    print("Analyzing results...")
    normalyzerResultsObject <- analyzeNormalizations(normalyzerResultsObject, jobName)
    print("Finished analysing results")
    
    print("Generating plots...")
    generatePlots(normalyzerResultsObject, jobDir)
    
    print(paste("Done! Results are stored in ", jobDir))
}


