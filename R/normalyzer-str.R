#' Normalyzer pipeline entry point
#' 
#' @param inputPath CSV delimited input file containing raw counts and 
#'  replicate information in header.
#' @param jobName Give the current run a name.
#' @param outputDir Specify an output directory for generated files.
#'  Defaults to current working directory.
#' @export
#' @examples
#' normalyzer("path/to/input.csv", "my_test_run")
#' normalyzer("path/to/input.csv", "my_test_run", outputDir="my_output")

normalyzer <- function(inputPath, jobName, outputDir=NULL) {
    
    print('start')
    
    # RcmdrMisc::numSummary used to summarize statistics
    require(Rcmdr)
    
    
    # Used for justvsn
    require(vsn)
    
    # Used for quantile normalization
    require(preprocessCore)

    # Used for normalizeCyclicLoess
    require(limma)
    
    # Used for rlm
    require(MASS)
    
    # ape::as.phylo
    require(ape)
    
    # Used for raster::cv - Coefficient of variation
    require(raster)

    # Used for car::showLabels
    require(car)
    
    # gridExtra::arrange
    require(gridExtra)
    
    # Extensively used for plotting
    require(ggplot2)
    
    # grid::grid.layout, and much more for plotting
    require(grid)

    # Biobase::rowMedians ??

    # require(PerformanceAnalytics)
    # require(abind)
    # require(e1071)
    
        
        
    source("generatePlots.R")
    source("normfinder-pipeline.R")
    source("normMethods.R")
    source("printMeta.R")
    source("printPlots.R")

    source("NormalyzerDataset.R")
    source("NormalizationEvaluationResults.R")
    source("NormalyzerResults.R")

    source("utils.R")
    source("inputVerification.R")
    source("analyzeResults.R")
    
    normObj <- getVerifiedNormalyzerObjectFromFile(inputPath, jobName)
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


