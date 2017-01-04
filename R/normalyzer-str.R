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
#' @param forceAllMethods Debugging function. Run all normalizations even if
#'  they aren't in the recommended range of number of values
#' @param omitLowAbundSamples Automatically remove samples with fewer non-NA
#'  values compared to threshold given by sampleAbundThres. 
#'  Will otherwise stop with error message if such sample is encountered.
#' @param sampleAbundThres Threshold for omitting low-abundant
#'  samples. Is by default set to 15.
#' @return None
#' @export
#' @import MASS Rcmdr limma preprocessCore methods
normalyzer <- function(inputPath, 
                       jobName, 
                       outputDir=NULL,
                       forceAllMethods=FALSE,
                       omitLowAbundSamples=FALSE,
                       sampleAbundThres=15,
                       requireReplicates=TRUE) {
    
    print('start')
    startTime <- Sys.time()
    
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
    
    normObj <- getVerifiedNormalyzerObject(inputPath, 
                                           jobName,
                                           threshold=sampleAbundThres,
                                           omitSamples=omitLowAbundSamples,
                                           requireReplicates=requireReplicates)
    jobDir <- setupJobDir(jobName, outputDir)
    
    
    print("Normalizing data....")
    normalyzerResultsObject <- normMethods(normObj, 
                                           jobName, 
                                           jobDir, 
                                           forceAll=forceAllMethods)
    print("Finished Normalization")
    
    print("Analyzing results...")
    normalyzerResultsObject <- analyzeNormalizations(normalyzerResultsObject, 
                                                     jobName)
    print("Finished analysing results")
    
    print("Generating plots...")
    generatePlots(normalyzerResultsObject, jobDir)
    
    print(paste("Done! Results are stored in ", jobDir))
    
    endTime <- Sys.time()
    elapsedTime <- elapsedSecondsBetweenSystimes(startTime, endTime)
    print(paste("Processing took", elapsedTime, "seconds"))
}


