#' NormalyzerDE pipeline entry point
#' 
#' This function is the main execution point for the normalization part of
#' the NormalyzerDE analysis pipeline. When executed it performs the following
#' steps:
#' 
#' 1: Loads the data matrix containing expression values and optional
#' annotations, as well as the design matrix containing the experimental setup
#' 2: Performs input data verification to validate that the data is in correct
#' format. This step captures many common formatting errors. It returns an
#' instance of the NormalyzerDataset class representing the unprocessed data.
#' 3: Calculate a range of normalizations for the dataset. The result is
#' provided as a NormalyzerResults object containing the resulting data matrices
#' from each normalization.
#' 4: Analyze the normalizations and generate performance measures for each
#' of the normalized datasets. This result is provided as a 
#' NormalyzerEvaluationResults object.
#' 5: Output the matrices containing the normalized datasets to files.
#' 6: Generate visualizations overviewing the performance measures and
#' write them to a PDF report.
#' 
#' @param jobName Give the current run a name.
#' @param designPath Path to file containing design matrix.
#' @param dataPath Specify an output directory for generated files.
#'  Defaults to current working directory.
#' @param experimentObj SummarizedExperiment object, can be provided as input
#'  as alternative to 'designPath' and 'dataPath'
#' @param outputDir Directory where results folder is created.
#' @param forceAllMethods Debugging function. Run all normalizations even if
#'  they aren't in the recommended range of number of values
#' @param omitLowAbundSamples Automatically remove samples with fewer non-NA
#'  values compared to threshold given by sampleAbundThres. 
#'  Will otherwise stop with error message if such sample is encountered.
#' @param sampleAbundThres Threshold for omitting low-abundant
#'  samples. Is by default set to 15.
#' @param tinyRunThres If total number of features is less than this, a limited
#'  run is performed.
#' @param requireReplicates Require multiple samples per condition to pass input 
#'  validation.
#' @param normalizeRetentionTime Perform normalizations over retention time.
#' @param plotRows Number of plot-rows in output documentation.
#' @param plotCols Number of plot-columns in output documentation.
#' @param zeroToNA Convert zero values to NA.
#' @param sampleColName Column name in design matrix containing sample IDs.
#' @param groupColName Column name in design matrix containing condition IDs.
#' @param inputFormat Type of input format.
#' @param skipAnalysis Only perform normalization steps.
#' @param quiet Omit status messages printed during run.
#' @param noLogTransform Don't log-transform the input.
#' 
#' @param rtStepSizeMinutes Retention time normalization window size.
#' @param rtWindowMinCount Minimum number of datapoints in each retention-time
#'   segment.
#' @param rtWindowShifts Number of layered retention time normalized windows.
#' @param rtWindowMergeMethod Merge approach for layered retention time windows.
#' 
#' @return None
#' @export
#' @import MASS limma preprocessCore methods RcmdrMisc
#' @importFrom raster cv
#' @examples
#' \dontrun{
#' data_path <- system.file(package="NormalyzerDE", "extdata", "tiny_data.tsv")
#' design_path <- system.file(package="NormalyzerDE", "extdata", "tiny_design.tsv")
#' out_dir <- tempdir()
#' normalyzer(
#'     jobName="my_jobname", 
#'     designPath=design_path, 
#'     dataPath=data_path, 
#'     outputDir=out_dir)
#' normalyzer(
#'     "my_jobname", 
#'     designMatrix="design.tsv", 
#'     "data.tsv", 
#'     outputDir="path/to/output", 
#'     normalizeRetentionTime=TRUE, 
#'     retentionTimeWindow=2)
#' normalyzer(
#'     "my_jobname", 
#'     designMatrix="design.tsv", 
#'     "data.tsv", 
#'     outputDir="path/to/output", 
#'     inputFormat="maxquantprot")
#' }
normalyzer <- function(
        jobName,
        designPath=NULL, 
        dataPath=NULL,
        experimentObj=NULL,
        outputDir=".",
        forceAllMethods=FALSE,
        omitLowAbundSamples=FALSE,
        sampleAbundThres=5,
        tinyRunThres=50,
        requireReplicates=TRUE,
        normalizeRetentionTime=TRUE,
        plotRows=3,
        plotCols=4,
        zeroToNA=FALSE,
        sampleColName="sample",
        groupColName="group",
        inputFormat="default",
        skipAnalysis=FALSE,
        quiet=FALSE,
        noLogTransform=FALSE,
        
        rtStepSizeMinutes=1,
        rtWindowMinCount=100,
        rtWindowShifts=1,
        rtWindowMergeMethod="mean"
    ) {
    
    if (!quiet) message("You are running version ", utils::packageVersion("NormalyzerDE"), " of NormalyzerDE")

    if (is.null(experimentObj) && (is.null(designPath) || is.null(dataPath))) {
        stop("Either options 'designPath' plus 'dataPath' or 'summarizedExp' need to be provided")
    }
    
    startTime <- Sys.time()
    
    if (!quiet) message("[Step 1/5] Load data and verify input")

    if (is.null(experimentObj)) {
        experimentObj <- setupRawDataObject(dataPath, designPath, inputFormat, 
                                            zeroToNA, sampleColName, groupColName)
    }
    else {
        verifySummarizedExperiment(experimentObj, sampleColName)
        SummarizedExperiment::colData(experimentObj)[[sampleColName]] <- 
            as.character(SummarizedExperiment::colData(experimentObj)[[sampleColName]])
        SummarizedExperiment::metadata(experimentObj) <- list(
            sample=sampleColName, 
            group=groupColName)
    }

    normObj <- getVerifiedNormalyzerObject(
        jobName=jobName,
        summarizedExp=experimentObj,
        threshold=sampleAbundThres,
        omitSamples=omitLowAbundSamples,
        requireReplicates=requireReplicates,
        quiet=quiet,
        noLogTransform=noLogTransform,
        tinyRunThres=tinyRunThres
    )
        
    jobDir <- setupJobDir(jobName, outputDir)
    if (!quiet) message(
        "[Step 1/5] Input verified, job directory prepared at:", 
        jobDir
    )
    
    if (!quiet) message("[Step 2/5] Performing normalizations")
    normalyzerResultsObject <- normMethods(
        normObj,
        forceAll=forceAllMethods,
        normalizeRetentionTime=normalizeRetentionTime,
        rtStepSizeMinutes=rtStepSizeMinutes,
        rtWindowMinCount=rtWindowMinCount,
        rtWindowShifts=rtWindowShifts,
        rtWindowMergeMethod=rtWindowMergeMethod,
        quiet=quiet,
        noLogTransform=noLogTransform
    )
    if (!quiet) message("[Step 2/5] Done!")
    
    if (!skipAnalysis) {
        if (!quiet) message("[Step 3/5] Generating evaluation measures...")
        normalyzerResultsObject <- analyzeNormalizations(normalyzerResultsObject)
        if (!quiet) message("[Step 3/5] Done!")
    }
    else {
        if (!quiet) message("[Step 3/5] skipAnalysis flag set so no analysis 
                          performed")
    }

    if (!quiet) message("[Step 4/5] Writing matrices to file")
    writeNormalizedDatasets(normalyzerResultsObject, jobDir)
    if (!quiet) message("[Step 4/5] Matrices successfully written")
    
    if (!skipAnalysis) {
        if (!quiet) message("[Step 5/5] Generating plots...")
        generatePlots(normalyzerResultsObject, jobDir, plotRows=plotRows, 
                      plotCols=plotCols)
        if (!quiet) message("[Step 5/5] Plots successfully generated")
    }
    else {
        if (!quiet) message("[Step 5/5] skipAnalysis flag set so no plots 
                          generated")
    }
    
    endTime <- Sys.time()
    totTime <- difftime(endTime, startTime, units="mins")
    if (!quiet) message("All done! Results are stored in: ", jobDir, 
                        ", processing time was ", round(totTime, 1), 
                        " minutes")
}

#' NormalyzerDE differential expression
#' 
#' Performs differential expression analysis on a normalization matrix.
#' This command executes a pipeline processing the data and generates an
#' annotated normalization matrix and a report containing p-value histograms
#' for each of the performed comparisons.
#' 
#' When executed, it performs the following steps:
#' 
#' 1: Read the data and the design matrices into dataframes.
#' 2: Generate an instance of the NormalyzerStatistics class representing the
#' data and their statistical comparisons.
#' 3: Optionally reduce technical replicates in both the data matrix and the 
#' design matrix
#' 4: Calculate statistical contrats between supplied groups
#' 5: Generate an annotated version of the original dataframe where columns
#' containing statistical key measures have been added
#' 6: Write the table to file
#' 7: Generate a PDF report displaying p-value histograms for each calculated
#' contrast
#' 
#' @param jobName Name of job
#' @param designPath File path to design matrix
#' @param dataPath File path to normalized matrix
#' @param experimentObj SummarizedExperiment object, can be provided as input
#'  as alternative to 'designPath' and 'dataPath'
#' @param comparisons Character vector containing target contrasts. 
#'   If comparing condA with condB, then the vector would be c("condA-condB")
#' @param outputDir Path to output directory
#' @param logTrans Log transform the input (needed if providing non-logged 
#'   input)
#' @param type Type of statistical comparison, "limma" or "welch"
#' @param sampleCol Design matrix column header for column containing sample IDs
#' @param condCol Design matrix column header for column containing sample 
#'   conditions
#' @param batchCol Provide an optional column for inclusion of possible batch 
#'   variance in the model
#' @param techRepCol Design matrix column header for column containing technical 
#'   replicates
#' @param leastRepCount Minimum required replicate count
#' @param quiet Omit status messages printed during run
#' 
#' @param sigThres Significance threshold use for illustrating significant hits
#'   in diagnostic plots
#' @param sigThresType Type of significance threshold, "fdr" or "p". "fdr" is
#'   strongly recommended (Benjamini-Hochberg corrected p-values)
#' @param log2FoldThres Fold-size cutoff for being considered significant in
#'   diagnostic plots
#' 
#' @return None
#' @export
#' @examples 
#' data_path <- system.file(package="NormalyzerDE", "extdata", "tiny_data.tsv")
#' design_path <- system.file(package="NormalyzerDE", "extdata", "tiny_design.tsv")
#' out_dir <- tempdir()
#' normalyzerDE(
#'   jobName="my_jobname", 
#'   comparisons=c("4-5"), 
#'   designPath=design_path, 
#'   dataPath=data_path,
#'   outputDir=out_dir,
#'   condCol="group")
#' @export
normalyzerDE <- function(jobName, comparisons, designPath=NULL, dataPath=NULL, experimentObj=NULL, 
                         outputDir=".", logTrans=FALSE, type="limma", sampleCol="sample", condCol="group", 
                         batchCol=NULL, techRepCol=NULL, leastRepCount=1, quiet=FALSE, 
                         sigThres=0.1, sigThresType="fdr", log2FoldThres=0) {

    if (!quiet) message("You are running version ", utils::packageVersion("NormalyzerDE"), " of NormalyzerDE")
    
    if (is.null(experimentObj) && (is.null(designPath) || is.null(dataPath))) {
        stop("Either options 'designPath' plus 'dataPath' or 'summarizedExp' need to be provided")
    }
    
    if (is.null(comparisons)) {
        stop("Argument 'comparisons' must be provided. Specify one or more comparisons as a vector.\n",
             "Example, one comparison between group 1 and 2: c('1-2')\n",
             "Example, two comparisons between groups 1 and 2, and groups 2 and 3: c('1-2', '2-3')")
    }
    
    startTime <- Sys.time()
    jobDir <- setupJobDir(jobName, outputDir)

    if (!quiet) print("Setting up statistics object")
    if (is.null(experimentObj)) {
        experimentObj <- setupRawContrastObject(dataPath, designPath, sampleCol)
    }
    else {
        verifySummarizedExperiment(experimentObj, sampleCol)
    }
    
    if (!is.null(techRepCol)) {
        if (!quiet) print("Reducing technical replicates")
        experimentObj <- reduceTechnicalReplicates(experimentObj, techRepCol, sampleCol)
    }

    nst <- NormalyzerStatistics(
        experimentObj, 
        logTrans=logTrans
    )
        
    if (!quiet) print("Calculating statistical contrasts...")
    nst <- calculateContrasts(nst, comparisons, type=type, condCol=condCol, batchCol=batchCol, leastRepCount=leastRepCount)
    if (!quiet) print("Contrast calculations done!")
    
    annotDf <- generateAnnotatedMatrix(nst)
    outPath <- paste0(jobDir, "/", jobName, "_stats.tsv")

    if (!quiet) print(paste("Writing", nrow(annotDf), "annotated rows to", outPath))
    utils::write.table(annotDf, file=outPath, sep="\t", row.names = FALSE, quote=FALSE)
    if (!quiet) print(paste("Writing statistics report"))
    generateStatsReport(nst, jobName, jobDir, sigThres, sigThresType, log2FoldThres)
    
    endTime <- Sys.time()
    totTime <- difftime(endTime, startTime, units="mins")
    if (!quiet) print(paste0("All done! Results are stored in: ", jobDir, ", processing time was ", round(totTime, 1), " minutes"))
}

