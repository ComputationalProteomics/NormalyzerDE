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
#' @import MASS limma preprocessCore methods RcmdrMisc
normalyzer <- function(inputPath, 
                       jobName, 
                       designMatrix=NULL,
                       outputDir=NULL,
                       forceAllMethods=FALSE,
                       omitLowAbundSamples=FALSE,
                       sampleAbundThres=15,
                       requireReplicates=TRUE,
                       normalizeRetentionTime=FALSE,
                       retentionTimeWindow=1,
                       include_pvals=FALSE,
                       pairwise_comparisons=NULL,
                       include_cv_col=FALSE,
                       categorical_anova=TRUE,
                       include_anova_p=FALSE,
                       var_filter_frac=NULL,
                       plot_rows=3,
                       plot_cols=4,
                       source_files=FALSE,
                       zeroToNA=FALSE,
                       inputFormat="default") {
    
    print('start')
    startTime <- Sys.time()
    
    # require(RcmdrMisc)    # RcmdrMisc::numSummary, CRAN (could maybe be omitted after re-implementation of its 'numSummary' function?)
    # require(vsn)          # Used for justvsn, BioConductor
    # require(preprocessCore)     # Used for quantile normalization, BC
    # require(limma)        # Used for normalizeCyclicLoess, BioConductor
    # require(MASS)         # Used for rlm, BioConductor
    # require(ape)          # ape::as.phylo
    # require(raster)       # Used for raster::cv - Coefficient of variation
    # require(car)          # Used for car::showLabels, BioConductor
    # require(gridExtra)    # gridExtra::arrange, CRAN
    # require(ggplot2)      # Extensively used for plotting, CRAN
    # require(grid)         # grid::grid.layout, and much more for plotting, base
    # require(hexbin)       # Needed for the meanSdPlot
    # require(raster)       # CRAN package    
    
    # Biobase::rowMedians ??
    
    
    if (source_files) {
        source("generatePlots.R")
        source("normfinder-pipeline.R")
        source("normMethods.R")
        source("higherOrderNormMethods.R")
        source("printMeta.R")
        source("printPlots.R")

        source("NormalyzerDataset.R")
        source("NormalizationEvaluationResults.R")
        source("NormalyzerResults.R")

        source("utils.R")
        source("inputVerification.R")
        source("analyzeResults.R")
        source("preparsers.R")
        source("outputUtils.R")
    }

    
    
    # source("evaluationMain.R")
    
    normObj <- getVerifiedNormalyzerObject(inputPath,
                                           jobName,
                                           designMatrix,
                                           threshold=sampleAbundThres,
                                           omitSamples=omitLowAbundSamples,
                                           requireReplicates=requireReplicates,
                                           zeroToNA=zeroToNA,
                                           inputFormat=inputFormat)
    jobDir <- setupJobDir(jobName, outputDir)
    
    print("Normalizing data...")
    normalyzerResultsObject <- normMethods(normObj,
                                           jobName,
                                           forceAll=forceAllMethods,
                                           normalizeRetentionTime=normalizeRetentionTime,
                                           retentionTimeWindow=retentionTimeWindow)
    print("Finished Normalization")
    
    print("Analyzing results...")
    normalyzerResultsObject <- analyzeNormalizations(normalyzerResultsObject, 
                                                     jobName, 
                                                     comparisons=pairwise_comparisons,
                                                     categorical_anova=categorical_anova,
                                                     var_filter_frac=var_filter_frac)
    print("Finished analysing results")

    print("Writing matrices to file")
    writeNormalizedDatasets(normalyzerResultsObject, 
                            jobDir, 
                            include_pvals=include_pvals, 
                            include_pairwise_comparisons=!is.null(pairwise_comparisons),
                            include_cv_col=include_cv_col,
                            include_anova_p=include_anova_p)

    print("Generating plots...")
    generatePlots(normalyzerResultsObject, jobDir, plot_rows=plot_rows, plot_cols=plot_cols)
    
    print(paste("Done! Results are stored in ", jobDir))
    
    endTime <- Sys.time()
    print(difftime(endTime, startTime))
    # elapsedTime <- elapsedSecondsBetweenSystimes(startTime, endTime)
    # print(paste("Processing took", elapsedTime, "seconds"))
}


