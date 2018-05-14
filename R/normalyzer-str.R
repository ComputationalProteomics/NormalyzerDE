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
#' @importFrom raster cv
#' @examples
#' normalyzer("data.tsv", "my_jobname", designMatrix="design.tsv", outputDir="path/to/output")
#' normalyzer("data.tsv", "my_jobname", designMatrix="design.tsv", outputDir="path/to/output", normalizeRetentionTime=TRUE, retentionTimeWindow=2)
#' normalyzer("data.tsv", "my_jobname", designMatrix="design.tsv", outputDir="path/to/output", inputFormat="maxquantprot")
normalyzer <- function(inputPath, 
                       jobName, 
                       designMatrix=NULL,
                       outputDir=NULL,
                       forceAllMethods=FALSE,
                       omitLowAbundSamples=FALSE,
                       sampleAbundThres=15,
                       requireReplicates=TRUE,
                       normalizeRetentionTime=TRUE,
                       retentionTimeWindow=1,
                       pairwiseComparisons=NULL,
                       includeCvCol=FALSE,
                       categoricalAnova=TRUE,
                       includeAnovaP=FALSE,
                       varFilterFrac=NULL,
                       plotRows=3,
                       plotCols=4,
                       sourceFiles=FALSE,
                       sourceBase=NULL,
                       zeroToNA=FALSE,
                       sampleColName="sample",
                       groupColName="group",
                       inputFormat="default",
                       skipAnalysis=FALSE) {

    startTime <- Sys.time()
    
    if (sourceFiles) {
        sourceScripts(sourceBase)
    }

    print("[Step 1/6] Verifying input")
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
    print(paste("[Step 1/6] Input verified, job directory prepared at:", jobDir))
    
    print("[Step 2/6] Performing normalizations")
    normalyzerResultsObject <- normMethods(normObj,
                                           forceAll=forceAllMethods,
                                           normalizeRetentionTime=normalizeRetentionTime,
                                           retentionTimeWindow=retentionTimeWindow)
    print("[Step 2/6] Done!")
    
    if (!skipAnalysis) {
        print("[Step 3/6] Generating evaluation measures...")
        normalyzerResultsObject <- analyzeNormalizations(normalyzerResultsObject, 
                                                         comparisons=pairwiseComparisons,
                                                         categoricalAnova=categoricalAnova,
                                                         varFilterFrac=varFilterFrac)
        print("[Step 3/6] Done!")
    }
    else {
        "[Step 3/6] skipAnalysis flag set so no analysis performed"
    }
    
    if (!skipAnalysis) {
        print("[Step 4/6] Calculating statistics...")
        # normalyzerResultsObject <- calculateStatistics(normalyzerResultsObject, )
        print("[Step 4/6] Plots successfully generated")
    }
    else {
        "[Step 4/6] skipAnalysis flag set so no statistics calculated"
    }

    print("[Step 5/6] Writing matrices to file")
    writeNormalizedDatasets(normalyzerResultsObject, 
                            jobDir, 
                            includePairwiseComparisons=!is.null(pairwiseComparisons),
                            includeCvCol=includeCvCol,
                            includeAnovaP=includeAnovaP)
    print("[Step 5/6] Matrices successfully written")
    
    if (!skipAnalysis) {
        print("[Step 6/6] Generating plots...")
        generatePlots(normalyzerResultsObject, jobDir, plot_rows=plotRows, plot_cols=plotCols)
        print("[Step 6/6] Plots successfully generated")
    }
    else {
        print("[Step 6/6] skipAnalysis flag set so no plots generated")
    }
    
    endTime <- Sys.time()
    totTime <- difftime(endTime, startTime, units="mins")
    print(paste0("All done! Results are stored in: ", jobDir, ", processing time was ", round(totTime, 1), " minutes"))
}

#' Normalyzer differential expression
#' 
#' @param dataFp File path to normalized matrix
#' @param designFp File path to design matrix
#' @param jobName Name of job
#' @param comparisons Character vector containing target contrasts. If comparing condA with condB, then the
#' vector would be c("condA-condB")
#' @param outputDir Path to output directory
#' @param logTrans Log transform the input (needed if providing non-logged input)
#' @param robustLimma Perform robust Limma estimate in the eBayes fit step
#' @param type Type of statistical comparison, "limma" or "welch"
#' @param batchCol Provide an optional column for inclusion of possible batch variance in the model
#' @return None
#' @export
#' @examples
#' normalyzerDE("results/normalized.tsv", "design.tsv", "my_jobname", c("1-2", "1-4"),  outputDir="path/to/output")
#' @export
normalyzerDE <- function(dataFp, designFp, jobName, comparisons, outputDir=NULL, logTrans=FALSE, 
                         robustLimma=FALSE, type="limma", sampleCol="sample", condCol="group", 
                         batchCol=NULL, techRepCol=NULL, leastRepCount=1, sourceFiles=FALSE,
                         sourceBase=NULL) {

    if (sourceFiles) {
        sourceScripts(sourceBase)
    }
    
    startTime <- Sys.time()
    jobDir <- setupJobDir(jobName, outputDir)

    print("Setting up statistics object")
    nst <- setupStatisticsObject(dataFp, designFp, comparisons, logTrans=logTrans, leastRepCount=leastRepCount)
    
    if (!is.null(techRepCol)) {
        print("Reducing technical replicates")
        nst@dataMat <- reduce_technical_replicates(nst@dataMat, nst@designDf[, techRepCol])
        nst@designDf <- reduce_design(nst@designDf, nst@designDf[, techRepCol])
    }
    
    print("Calculating statistical contrasts...")
    nst <- calculateContrasts(nst, comparisons, condCol="group", type=type, batchCol=batchCol)
    print("Contrast calculations done!")
    
    annotDf <- generateAnnotatedMatrix(nst)
    outPath <- paste0(jobDir, "/", jobName, "_stats.tsv")
    
    print(paste("Writing", nrow(annotDf), "annotated rows to", outPath))
    write.table(annotDf, file=outPath, sep="\t", row.names = F)
    print(paste("Writing statistics report"))
    generateStatsReport(nst, jobName, jobDir)
    
    endTime <- Sys.time()
    totTime <- difftime(endTime, startTime, units="mins")
    print(paste0("All done! Results are stored in: ", jobDir, ", processing time was ", round(totTime, 1), " minutes"))
}

sourceScripts <- function(sourceBase) {
    source(paste(sourceBase, "generatePlots.R", sep="/"))
    source(paste(sourceBase, "normfinder-pipeline.R", sep="/"))
    source(paste(sourceBase, "normMethods.R", sep="/"))
    source(paste(sourceBase, "higherOrderNormMethods.R", sep="/"))
    source(paste(sourceBase, "printMeta.R", sep="/"))
    source(paste(sourceBase, "printPlots.R", sep="/"))
    source(paste(sourceBase, "calculateStatistics.R", sep="/"))
    
    source(paste(sourceBase, "NormalyzerDataset.R", sep="/"))
    source(paste(sourceBase, "NormalizationEvaluationResults.R", sep="/"))
    source(paste(sourceBase, "NormalyzerResults.R", sep="/"))
    source(paste(sourceBase, "NormalyzerStatistics.R", sep="/"))
    
    source(paste(sourceBase, "utils.R", sep="/"))
    source(paste(sourceBase, "inputVerification.R", sep="/"))
    source(paste(sourceBase, "analyzeResults.R", sep="/"))
    source(paste(sourceBase, "preparsers.R", sep="/"))
    source(paste(sourceBase, "outputUtils.R", sep="/"))
}















