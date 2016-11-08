
print('sourced')

normalyzer <- function(datafile, getjob){
    
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

    source("NormalyzerObject.R")

    # stopFunction()
    
    print("Normalizing data....")
    # try.result <- try(normalizeddata <- normMethods(datafile, getjob))
    normalizeddata <- normMethods(datafile, getjob)
    
    # if (inherits(try.result, "try-error")) {
    #     return(try.result)
    # }
    
    print("Finished Normalization")
    print("Analyzing data....")
    
    # try.result <- try(analyzeAndPlot(normalizeddata, getjob))
    analyzeAndPlot(normalizeddata, getjob)
    
    # if(inherits(try.result, "try-error")){
    #     return(try.result)
    # }
    
    print("Done! Results are stored in the working directory")
}

myprint <- function(..., debug=T) {
    if (debug) {
        print(paste("DEBUG: ", ...))
    }
    else {
        print(paste(...))
    }
}

varprint <- function(variable, debug=T) {
    varName <- deparse(substitute(variable))
    
    if (debug) {
        cat(varName, variable, sprintf("[%s]", typeof(variable)), sprintf("[%s elements]", length(variable)), "\n")
    }
    else {
        cat("DEBUG: ", varName, variable, sprintf("[%s]", typeof(variable)), sprintf("[%s elements]", length(variable)), "\n")
    }
}

