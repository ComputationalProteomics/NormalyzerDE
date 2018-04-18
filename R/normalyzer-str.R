#' Normalyzer pipeline entry point
#' 
#' @param inputPath CSV delimited input file containing raw counts and 
#'  replicate information in header.
#' @param jobName Give the current run a name.
#' @param designMatrix Path to file containing design matrix.
#' @param outputDir Specify an output directory for generated files.
#'  Defaults to current working directory.
#' @param forceAllMethods Debugging function. Run all normalizations even if
#'  they aren't in the recommended range of number of values
#' @param omitLowAbundSamples Automatically remove samples with fewer non-NA
#'  values compared to threshold given by sampleAbundThres. 
#'  Will otherwise stop with error message if such sample is encountered.
#' @param sampleAbundThres Threshold for omitting low-abundant
#'  samples. Is by default set to 15.
#' @param requireReplicates Require multiple samples per condition to pass input validation.
#' @param normalizeRetentionTime Perform normalizations over retention time.
#' @param retentionTimeWindow Retention time normalization window size.
#' @param include_pvals Include p-values in output (is this used?)
#' @param pairwise_comparisons Vector of pairwise sample comparisons.
#' @param include_cv_col Include CV in output matrices.
#' @param categorical_anova Perform categorical ANOVA (alternatively numerical).
#' @param include_anova_p Include ANOVA p-value in output.
#' @param var_filter_frac Perform variance filtering of highly variant features.
#' @param plot_rows Number of plot-rows in output documentation.
#' @param plot_cols Number of plot-columns in output documentation.
#' @param source_files Source code files (for development purposes only).
#' @param source_base Base to source code files from (for development purposes only).
#' @param zeroToNA Convert zero values to NA.
#' @param sampleColName Column name in design matrix containing sample IDs.
#' @param groupColName Column name in design matrix containing condition IDs.
#' @param inputFormat Type of input format.
#' @param skipAnalysis Only perform normalization steps.
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
                       pairwise_comparisons=NULL,
                       include_cv_col=FALSE,
                       categorical_anova=TRUE,
                       include_anova_p=FALSE,
                       var_filter_frac=NULL,
                       plot_rows=3,
                       plot_cols=4,
                       source_files=FALSE,
                       source_base=NULL,
                       zeroToNA=FALSE,
                       sampleColName="sample",
                       groupColName="group",
                       inputFormat="default",
                       skipAnalysis=FALSE) {

    startTime <- Sys.time()
    
    if (source_files) {
        
        source(paste(source_base, "generatePlots.R", sep="/"))
        source(paste(source_base, "normfinder-pipeline.R", sep="/"))
        source(paste(source_base, "normMethods.R", sep="/"))
        source(paste(source_base, "higherOrderNormMethods.R", sep="/"))
        source(paste(source_base, "printMeta.R", sep="/"))
        source(paste(source_base, "printPlots.R", sep="/"))

        source(paste(source_base, "NormalyzerDataset.R", sep="/"))
        source(paste(source_base, "NormalizationEvaluationResults.R", sep="/"))
        source(paste(source_base, "NormalyzerResults.R", sep="/"))

        source(paste(source_base, "utils.R", sep="/"))
        source(paste(source_base, "inputVerification.R", sep="/"))
        source(paste(source_base, "analyzeResults.R", sep="/"))
        source(paste(source_base, "preparsers.R", sep="/"))
        source(paste(source_base, "outputUtils.R", sep="/"))
    }

    print("[Step 1/5] Verifying input")
    normObj <- getVerifiedNormalyzerObject(inputPath,
                                           jobName,
                                           designMatrix,
                                           threshold=sampleAbundThres,
                                           omitSamples=omitLowAbundSamples,
                                           requireReplicates=requireReplicates,
                                           zeroToNA=zeroToNA,
                                           inputFormat=inputFormat,
                                           sampleCol=sampleColName,
                                           groupCol=groupColName)
    jobDir <- setupJobDir(jobName, outputDir)
    print(paste("[Step 1/5] Input verified, job directory prepared at:", jobDir))
    
    print("[Step 2/5] Performing normalizations")
    normalyzerResultsObject <- normMethods(normObj,
                                           forceAll=forceAllMethods,
                                           normalizeRetentionTime=normalizeRetentionTime,
                                           retentionTimeWindow=retentionTimeWindow)
    print("[Step 2/5] Done!")
    
    if (!skipAnalysis) {
        print("[Step 3/5] Generating evaluation measures...")
        normalyzerResultsObject <- analyzeNormalizations(normalyzerResultsObject, 
                                                         comparisons=pairwise_comparisons,
                                                         categorical_anova=categorical_anova,
                                                         var_filter_frac=var_filter_frac)
        print("[Step 3/5] Done!")
    }
    else {
        "[Step 3/5] skipAnalysis flag set so no analysis performed"
    }

    print("[Step 4/5] Writing matrices to file")
    writeNormalizedDatasets(normalyzerResultsObject, 
                            jobDir, 
                            include_pairwise_comparisons=!is.null(pairwise_comparisons),
                            include_cv_col=include_cv_col,
                            include_anova_p=include_anova_p)
    print("[Step 4/5] Matrices successfully written")
    
    if (!skipAnalysis) {
        print("[Step 5/5] Generating plots...")
        generatePlots(normalyzerResultsObject, jobDir, plot_rows=plot_rows, plot_cols=plot_cols)
        print("[Step 5/5] Plots successfully generated")
    }
    else {
        print("[Step 5/5] skipAnalysis flag set so no analysis performed - skipping evaluation plots")
    }
    
    endTime <- Sys.time()
    totTime <- difftime(endTime, startTime)
    print(paste0("All done! Results are stored in: ", jobDir, ", processing time was ", round(totTime, 1), " seconds"))
    
}


