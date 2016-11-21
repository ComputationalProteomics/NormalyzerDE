
print('Normalyzer sourced')

setupJobDir <- function(jobName, outputDir) {
    
    print("DEBUG: Setup job dir called")
    
    varprint(jobName)
    varprint(outputDir)
    
    if (is.null(outputDir)) {
        jobDir <- paste(getwd(), "/", jobName[1], sep="")
    }
    else {
        jobDir <- paste(outputDir, "/", jobName[1], sep="")
    }
    createDirectory(jobDir)
    
    varprint(jobDir)
    
    jobDir
}

createDirectory <- function(targetPath) {
    
    if (file.exists(targetPath)) {
        abc <- "Directory already exists"
        class(abc) <- "try-error"
        if (inherits(abc, "try-error")) {
            return(abc)
        }
        stop("Directory already exists")
    } 
    else {
        dir.create(targetPath)
    }
}

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
    
    source("analyzeAndPlot.R")
    source("normfinder-pipeline.R")
    source("normMethods.R")
    source("printMeta.R")
    source("printPlots.R")

    source("NormalyzerDataset.R")
    source("NormalyzerResults.R")
    
    source("utils.R")
    source("inputVerification.R")
    
    # stopFunction()
    
    normObj <- getVerifiedNormalyzerObjectFromFile(datafile, jobName)

    jobDir <- setupJobDir(jobName, outputDir)
        
    print("Normalizing data....")
    # try.result <- try(normalizeddata <- normMethods(datafile, getjob))
    
    # normalizeddata <- normMethods(datafile, getjob, outputDir=outputDir)
    # normalizeddata <- normMethods(normObj, jobName, jobDir)
    normalyzerResultsObject <- normMethods(normObj, jobName, jobDir)
    
    # if (inherits(try.result, "try-error")) {
    #     return(try.result)
    # }
    
    print("Finished Normalization")
    print("Analyzing data....")
    
    # try.result <- try(analyzeAndPlot(normalizeddata, getjob))
    # analyzeAndPlot(normalizeddata, jobName, outputDir=outputDir)
    analyzeAndPlot(normalyzerResultsObject, jobName, outputDir=outputDir)
    
    
        
    # if(inherits(try.result, "try-error")){
    #     return(try.result)
    # }
    
    print("Done! Results are stored in the working directory")
}


