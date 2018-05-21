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
#' @param pairwiseComparisons Vector of pairwise sample comparisons.
#' @param includeCvCol Include CV in output matrices.
#' @param categoricalAnova Perform categorical ANOVA (alternatively numerical).
#' @param includeAnovaP Include ANOVA p-value in output.
#' @param varFilterFrac Perform variance filtering of highly variant features.
#' @param plotRows Number of plot-rows in output documentation.
#' @param plotCols Number of plot-columns in output documentation.
#' @param sourceFiles Source code files (for development purposes only).
#' @param sourceBase Base to source code files from (for development purposes only).
#' @param zeroToNA Convert zero values to NA.
#' @param sampleColName Column name in design matrix containing sample IDs.
#' @param groupColName Column name in design matrix containing condition IDs.
#' @param inputFormat Type of input format.
#' @param skipAnalysis Only perform normalization steps.
#' @param quiet Omit status messages printed during run
#' @return None
#' @export
#' @import MASS limma preprocessCore methods RcmdrMisc
#' @importFrom raster cv
#' @examples
#' normalyzer("data.tsv", "my_jobname", designMatrix="design.tsv", outputDir="path/to/output")
#' normalyzer("data.tsv", "my_jobname", designMatrix="design.tsv", outputDir="path/to/output", normalizeRetentionTime=TRUE, retentionTimeWindow=2)
#' normalyzer("data.tsv", "my_jobname", designMatrix="design.tsv", outputDir="path/to/output", inputFormat="maxquantprot")
normalyzer <- function(jobName,
                       designPath, 
                       dataPath,
                       outputDir=".",
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
                       sourceBase=".",
                       zeroToNA=FALSE,
                       sampleColName="sample",
                       groupColName="group",
                       inputFormat="default",
                       skipAnalysis=FALSE,
                       quiet=FALSE) {

    startTime <- Sys.time()
    
    if (sourceFiles) {
        sourceScripts(sourceBase)
    }

    if (!quiet) print("[Step 1/5] Verifying input")
    normObj <- getVerifiedNormalyzerObject(jobName=jobName,
                                           designPath=designPath,
                                           dataPath=dataPath,
                                           threshold=sampleAbundThres,
                                           omitSamples=omitLowAbundSamples,
                                           requireReplicates=requireReplicates,
                                           zeroToNA=zeroToNA,
                                           inputFormat=inputFormat,
                                           sampleCol=sampleColName,
                                           groupCol=groupColName)
    jobDir <- setupJobDir(jobName, outputDir)
    if (!quiet) print(paste("[Step 1/5] Input verified, job directory prepared at:", jobDir))
    
    if (!quiet) print("[Step 2/5] Performing normalizations")
    normalyzerResultsObject <- normMethods(normObj,
                                           forceAll=forceAllMethods,
                                           normalizeRetentionTime=normalizeRetentionTime,
                                           retentionTimeWindow=retentionTimeWindow)
    if (!quiet) print("[Step 2/5] Done!")
    
    if (!skipAnalysis) {
        if (!quiet) print("[Step 3/5] Generating evaluation measures...")
        normalyzerResultsObject <- analyzeNormalizations(normalyzerResultsObject, 
                                                         comparisons=pairwiseComparisons,
                                                         categoricalAnova=categoricalAnova,
                                                         varFilterFrac=varFilterFrac)
        if (!quiet) print("[Step 3/5] Done!")
    }
    else {
        if (!quiet) print("[Step 3/5] skipAnalysis flag set so no analysis performed")
    }
    
    if (!quiet) print("[Step 4/5] Writing matrices to file")
    writeNormalizedDatasets(normalyzerResultsObject, 
                            jobDir, 
                            includePairwiseComparisons=!is.null(pairwiseComparisons),
                            includeCvCol=includeCvCol,
                            includeAnovaP=includeAnovaP)
    if (!quiet) print("[Step 4/5] Matrices successfully written")
    
    if (!skipAnalysis) {
        if (!quiet) print("[Step 5/5] Generating plots...")
        generatePlots(normalyzerResultsObject, jobDir, plotRows=plotRows, plotCols=plotCols)
        if (!quiet) print("[Step 5/5] Plots successfully generated")
    }
    else {
        if (!quiet) print("[Step 5/5] skipAnalysis flag set so no plots generated")
    }
    
    endTime <- Sys.time()
    totTime <- difftime(endTime, startTime, units="mins")
    if (!quiet) print(paste0("All done! Results are stored in: ", jobDir, ", processing time was ", round(totTime, 1), " minutes"))
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
#' @param quiet Omit status messages printed during run
#' @return None
#' @export
#' @examples
#' normalyzerDE("results/normalized.tsv", "design.tsv", "my_jobname", c("1-2", "1-4"),  outputDir="path/to/output")
#' @export
normalyzerDE <- function(jobName, designPath, dataPath, comparisons, outputDir=".", logTrans=FALSE, 
                         robustLimma=FALSE, type="limma", sampleCol="sample", condCol="group", 
                         batchCol=NULL, techRepCol=NULL, leastRepCount=1, sourceFiles=FALSE,
                         sourceBase=NULL, quiet=FALSE) {

    if (sourceFiles) {
        sourceScripts(sourceBase)
    }
    
    startTime <- Sys.time()
    jobDir <- setupJobDir(jobName, outputDir)

    if (!quiet) print("Setting up statistics object")
    nst <- setupStatisticsObject(designPath, dataPath, comparisons, logTrans=logTrans, leastRepCount=leastRepCount)
    
    if (!is.null(techRepCol)) {
        if (!quiet) print("Reducing technical replicates")
        nst@dataMat <- reduceTechnicalReplicates(nst@dataMat, nst@designDf[, techRepCol])
        nst@designDf <- reduceDesign(nst@designDf, nst@designDf[, techRepCol])
    }
    
    if (!quiet) print("Calculating statistical contrasts...")
    nst <- calculateContrasts(nst, comparisons, condCol="group", type=type, batchCol=batchCol)
    if (!quiet) print("Contrast calculations done!")
    
    annotDf <- generateAnnotatedMatrix(nst)
    outPath <- paste0(jobDir, "/", jobName, "_stats.tsv")
    
    if (!quiet) print(paste("Writing", nrow(annotDf), "annotated rows to", outPath))
    write.table(annotDf, file=outPath, sep="\t", row.names = F)
    if (!quiet) print(paste("Writing statistics report"))
    generateStatsReport(nst, jobName, jobDir)
    
    endTime <- Sys.time()
    totTime <- difftime(endTime, startTime, units="mins")
    if (!quiet) print(paste0("All done! Results are stored in: ", jobDir, ", processing time was ", round(totTime, 1), " minutes"))
}

sourceScripts <- function(sourceBase) {
    source(paste(sourceBase, "generatePlots.R", sep="/"))
    source(paste(sourceBase, "normMethods.R", sep="/"))
    source(paste(sourceBase, "higherOrderNormMethods.R", sep="/"))
    source(paste(sourceBase, "printMeta.R", sep="/"))
    source(paste(sourceBase, "printPlots.R", sep="/"))
    source(paste(sourceBase, "calculateStatistics.R", sep="/"))
    
    source(paste(sourceBase, "NormalyzerDataset.R", sep="/"))
    source(paste(sourceBase, "NormalyzerEvaluationResults.R", sep="/"))
    source(paste(sourceBase, "NormalyzerResults.R", sep="/"))
    source(paste(sourceBase, "NormalyzerStatistics.R", sep="/"))
    
    source(paste(sourceBase, "utils.R", sep="/"))
    source(paste(sourceBase, "inputVerification.R", sep="/"))
    source(paste(sourceBase, "analyzeResults.R", sep="/"))
    source(paste(sourceBase, "preparsers.R", sep="/"))
    source(paste(sourceBase, "outputUtils.R", sep="/"))
}

