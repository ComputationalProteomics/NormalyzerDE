# require("grDevices")
# require("graphics")
# require("methods")
# require("stats")
# require("utils")



#' Normalyzer pipeline entry point
#' 
#' @param inputPath CSV delimited input file containing raw counts and 
#'  replicate information in header.
#' @param jobName Give the current run a name.
#' @param outputDir Specify an output directory for generated files.
#'  Defaults to current working directory.
#' @return None
#' @export
#' @import MASS Rcmdr limma preprocessCore methods
normalyzer <- function(inputPath, jobName, outputDir=NULL) {
    
    print('start')
    
    # require(Rcmdr)  # RcmdrMisc::numSummary, CRAN
    # require(vsn)    # Used for justvsn, BioConductor
    # require(preprocessCore)     # Used for quantile normalization, BC
    # require(limma)  # Used for normalizeCyclicLoess, BioConductor
    # require(MASS)   # Used for rlm, BioConductor
    # require(ape)    # ape::as.phylo
    # require(raster) # Used for raster::cv - Coefficient of variation
    # require(car)    # Used for car::showLabels, BioConductor
    # require(gridExtra)  # gridExtra::arrange, CRAN
    # require(ggplot2)  # Extensively used for plotting, CRAN
    # require(grid)     # grid::grid.layout, and much more for plotting, base
    
    
    # Biobase::rowMedians ??
    

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
    normalyzerResultsObject <- analyzeNormalizations(normalyzerResultsObject, 
                                                     jobName)
    print("Finished analysing results")
    
    print("Generating plots...")
    generatePlots(normalyzerResultsObject, jobDir)
    
    print(paste("Done! Results are stored in ", jobDir))
}


