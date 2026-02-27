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
#' @param dataPath Path to file containing data matrix.
#' @param experimentObj SummarizedExperiment object, can be provided as input
#'  as alternative to 'designPath' and 'dataPath'
#' @param outputDir Directory where results folder is created.
#' @param forceAllMethods Debugging function. Run all normalizations even if
#'  they aren't in the recommended range of number of values
#' @param omitLowAbundSamples Automatically remove samples with fewer non-NA
#'  values compared to threshold given by sampleAbundThres.
#'  Will otherwise stop with error message if such sample is encountered.
#' @param sampleAbundThres Threshold for omitting low-abundant
#'  samples. Is by default set to 5.
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
#' @param inputOptions Optional list of input-reader options. For DIA-NN, use
#'   \code{\link{diannInputOptions}}. For delimited inputs, use
#'   \code{\link{defaultInputOptions}}, \code{\link{proteiosInputOptions}}, or
#'   \code{\link{maxQuantInputOptions}} to configure the delimiter.
#' @param skipAnalysis Only perform normalization steps.
#' @param quiet Omit status messages printed during run.
#' @param noLogTransform Don't log-transform the input.
#' @param writeReportAsPngs Output the evaluation report as PNG files instead of
#'  a single PDF
#'
#' @param rtStepSizeMinutes Retention time normalization window size.
#' @param rtWindowMinCount Minimum number of datapoints in each retention-time
#'   segment.
#' @param rtWindowShifts Number of layered retention time normalized windows.
#' @param rtWindowMergeMethod Merge approach for layered retention time windows.
#'
#' @return None
#' @export
#' @import MASS limma methods
#' @examples
#' data_path <- system.file(package="NormalyzerDE", "extdata", "tiny_data.tsv")
#' design_path <- system.file(package="NormalyzerDE", "extdata", "tiny_design.tsv")
#' out_dir <- tempdir()
#' job_name <- basename(tempfile("job_"))
#' normalyzer(
#'     jobName=job_name,
#'     designPath=design_path,
#'     dataPath=data_path,
#'     outputDir=out_dir,
#'     skipAnalysis=TRUE,
#'     quiet=TRUE)
normalyzer <- function(
  jobName,
  designPath = NULL,
  dataPath = NULL,
  experimentObj = NULL,
  outputDir = ".",
  forceAllMethods = FALSE,
  omitLowAbundSamples = FALSE,
  sampleAbundThres = 5,
  tinyRunThres = 50,
  requireReplicates = TRUE,
  normalizeRetentionTime = TRUE,
  plotRows = 3,
  plotCols = 4,
  zeroToNA = FALSE,
  sampleColName = "sample",
  groupColName = "group",
  inputFormat = "default",
  inputOptions = NULL,
  skipAnalysis = FALSE,
  quiet = FALSE,
  noLogTransform = FALSE,
  writeReportAsPngs = FALSE,

  rtStepSizeMinutes = 1,
  rtWindowMinCount = 100,
  rtWindowShifts = 1,
  rtWindowMergeMethod = "mean"
) {
  if (!quiet) {
    message(
      "You are running version ",
      utils::packageVersion("NormalyzerDE"),
      " of NormalyzerDE"
    )
  }

  if (is.null(experimentObj) && (is.null(designPath) || is.null(dataPath))) {
    stop(
      "Either options 'designPath' plus 'dataPath' or 'summarizedExp' need to be provided"
    )
  }

  startTime <- Sys.time()

  if (!quiet) {
    message("[Step 1/5] Load data and verify input")
  }

  if (is.null(experimentObj)) {
    experimentObj <- setupRawDataObject(
      dataPath = dataPath,
      designPath = designPath,
      inputFormat = inputFormat,
      zeroToNA = zeroToNA,
      sampleColName = sampleColName,
      groupColName = groupColName,
      inputOptions = inputOptions
    )
  } else {
    verifySummarizedExperiment(experimentObj, sampleColName)
    SummarizedExperiment::colData(experimentObj)[[sampleColName]] <-
      as.character(SummarizedExperiment::colData(experimentObj)[[
        sampleColName
      ]])
    S4Vectors::metadata(experimentObj) <- list(
      sample = sampleColName,
      group = groupColName
    )
  }

  normObj <- getVerifiedNormalyzerObject(
    jobName = jobName,
    summarizedExp = experimentObj,
    threshold = sampleAbundThres,
    omitSamples = omitLowAbundSamples,
    requireReplicates = requireReplicates,
    quiet = quiet,
    noLogTransform = noLogTransform,
    tinyRunThres = tinyRunThres
  )

  jobDir <- setupJobDir(jobName, outputDir)
  if (!quiet) {
    message(
      "[Step 1/5] Input verified, job directory prepared at:",
      jobDir
    )
  }

  if (!quiet) {
    message("[Step 2/5] Performing normalizations")
  }
  normalyzerResultsObject <- normMethods(
    normObj,
    forceAll = forceAllMethods,
    normalizeRetentionTime = normalizeRetentionTime,
    rtStepSizeMinutes = rtStepSizeMinutes,
    rtWindowMinCount = rtWindowMinCount,
    rtWindowShifts = rtWindowShifts,
    rtWindowMergeMethod = rtWindowMergeMethod,
    quiet = quiet,
    noLogTransform = noLogTransform
  )
  if (!quiet) {
    message("[Step 2/5] Done!")
  }

  if (!skipAnalysis) {
    if (!quiet) {
      message("[Step 3/5] Generating evaluation measures...")
    }
    normalyzerResultsObject <- analyzeNormalizations(normalyzerResultsObject)
    if (!quiet) message("[Step 3/5] Done!")
  } else {
    if (!quiet) {
      message(
        "[Step 3/5] skipAnalysis flag set so no analysis 
                          performed"
      )
    }
  }

  if (!quiet) {
    message("[Step 4/5] Writing matrices to file")
  }
  writeNormalizedDatasets(normalyzerResultsObject, jobDir)
  if (!quiet) {
    message("[Step 4/5] Matrices successfully written")
  }

  if (!skipAnalysis) {
    if (!quiet) {
      message("[Step 5/5] Generating plots...")
    }
    generatePlots(
      normalyzerResultsObject,
      jobDir,
      plotRows = plotRows,
      plotCols = plotCols,
      writeAsPngs = writeReportAsPngs
    )
    if (!quiet) message("[Step 5/5] Plots successfully generated")
  } else {
    if (!quiet) {
      message(
        "[Step 5/5] skipAnalysis flag set so no plots 
                          generated"
      )
    }
  }

  endTime <- Sys.time()
  totTime <- difftime(endTime, startTime, units = "mins")
  if (!quiet) {
    message(
      "All done! Results are stored in: ",
      jobDir,
      ", processing time was ",
      round(totTime, 1),
      " minutes"
    )
  }
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
#' @details
#' For \code{type="limpa"}, NormalyzerDE uses the Bioconductor package
#' \pkg{limpa} to model intensity-dependent missing values via a detection
#' probability curve (DPC) and to propagate quantification uncertainty into the
#' differential expression analysis. The input should be on the log2 scale with
#' missing values encoded as \code{NA}. Between-sample normalization (if desired)
#' can be performed upstream (for example by the quantification tool or via
#' \code{\link{normalyzer}}) or after \code{limpa::dpcQuant()} using
#' \code{limpaPostQuantNorm}. Avoid applying multiple normalizations
#' unintentionally.
#'
#' For PTM-level data (e.g., phosphoproteomics) where each row corresponds to a
#' modified site, set \code{limpaByRow=TRUE} (or \code{limpaProteinIdCol=NULL})
#' to keep each row separate rather than summarizing to protein-level.
#'
#' By default, NormalyzerDE uses a fixed DPC slope
#' (\code{limpaQuantArgs$dpc.slope}, default \code{0.8}) and lets limpa estimate
#' the intercept. To estimate both DPC parameters from your data, set
#' \code{limpaDpcMethod="dpc"} (or \code{"dpcCN"} for a complete-normal estimate,
#' which can be more robust for datasets with very large fold-changes).
#'
#' @param jobName Name of job
#' @param designPath File path to design matrix
#' @param dataPath File path to normalized matrix
#' @param experimentObj SummarizedExperiment object, can be provided as input
#'  as alternative to 'designPath' and 'dataPath'
#' @param comparisons Character vector containing target contrasts.
#'   If comparing condA with condB, then the vector would be c("condA-condB").
#'   Ignored if \code{oneVsRest=TRUE}.
#' @param outputDir Path to output directory
#' @param logTrans Log transform the input (needed if providing non-logged
#'   input)
#' @param type Type of statistical comparison, "limma", "limma_intensity" or
#'  "welch" or "limpa", where "limma_intensity" allows the prior to be fit
#'  according to intensity rather than using a flat prior. "limpa" uses the
#'  optional Bioconductor package \pkg{limpa} to handle missing values via a
#'  detection probability curve (DPC) model.
#' @param sampleCol Design matrix column header for column containing sample IDs
#' @param condCol Design matrix column header for column containing sample
#'   conditions
#' @param batchCol Provide an optional column for inclusion of possible batch
#'   variance in the model
#' @param techRepCol Design matrix column header for column containing technical
#'   replicates
#' @param leastRepCount Minimum required replicate count. For \code{type="limpa"},
#'   a feature is retained if at least one group has \code{leastRepCount}
#'   observed samples (features entirely missing across all samples are removed).
#' @param impute Whether to impute values (ignored for \code{type="limpa"}).
#' @param imputeMinFraction Minimum fraction non-NA values for an analyte in
#'   any group to impute in other groups (ignored for \code{type="limpa"}).
#' @param subsetByComparison If TRUE, subset data and design to each comparison
#'   before NA-filtering, imputation and model fitting.
#' @param oneVsRest If TRUE, compute one-vs-rest contrasts for each group in
#'   \code{condCol} (or the subset in \code{oneVsRestGroups}).
#' @param oneVsRestGroups Optional character vector specifying which groups in
#'   \code{condCol} to compare against all other samples.
#' @param limpaProteinIdCol For \code{type="limpa"}, optionally summarize
#'   peptide/precursor rows to protein-level using \code{limpa::dpcQuant()}.
#'   Set to a column name in the row annotation (for example \code{"Protein.Group"})
#'   to use as the protein identifier. Use \code{"auto"} (default) to try common
#'   identifiers. If the chosen column contains duplicate identifiers, the data
#'   are summarized once across all samples and the output rows correspond to
#'   proteins. Set to \code{NULL} to disable protein summarization and treat each
#'   row as one protein (recommended for PTM-level data such as phosphoproteomics
#'   where each row corresponds to a modified site).
#' @param limpaByRow For \code{type="limpa"}, treat each input row as a separate
#'   protein and always use \code{limpa::dpcQuantByRow()} instead of summarizing
#'   via \code{limpa::dpcQuant()}. This is recommended for PTM-level matrices
#'   (e.g., phosphosites). Equivalent to setting \code{limpaProteinIdCol=NULL}.
#' @param limpaDpc For \code{type="limpa"}, optional DPC parameters to pass to
#'   \code{limpa::dpcQuant()} / \code{limpa::dpcQuantByRow()}. Can be a list as
#'   returned by \code{limpa::dpc()}, or a numeric vector \code{c(beta0, beta1)}.
#' @param limpaDpcMethod For \code{type="limpa"}, optional method to estimate the
#'   DPC parameters from the data when \code{limpaDpc} is not supplied.
#'   \code{"none"} (default) uses a fixed slope (\code{limpaQuantArgs$dpc.slope},
#'   default \code{0.8}) and lets limpa estimate the intercept internally.
#'   \code{"dpc"} estimates both DPC parameters from the observed-normal model
#'   via \code{limpa::dpc()}. \code{"dpcCN"} estimates the DPC from the
#'   complete-normal model via \code{limpa::dpcCN()}, which can be more robust
#'   for datasets with very large fold-changes.
#' @param limpaDpcArgs For \code{type="limpa"}, optional named list of additional
#'   arguments forwarded to \code{limpa::dpc()} or \code{limpa::dpcCN()} when
#'   \code{limpaDpcMethod} is not \code{"none"}. Argument \code{y} is ignored.
#'   For \code{limpaDpcMethod="dpcCN"}, \code{dpc.slope.start} defaults to
#'   \code{limpaQuantArgs$dpc.slope}.
#' @param limpaQuantArgs For \code{type="limpa"}, optional named list of
#'   additional arguments forwarded to \code{limpa::dpcQuant()} /
#'   \code{limpa::dpcQuantByRow()}. Use this to set \code{dpc.slope} (default
#'   \code{0.8}), \code{chunk} (default \code{1000L}), and \code{verbose} (default
#'   \code{FALSE}), plus any additional \code{...} arguments supported by limpa.
#'   Arguments \code{y}, \code{protein.id}, and \code{dpc} are ignored.
#' @param limpaPostQuantNorm For \code{type="limpa"}, optional between-sample
#'   normalization applied to the quantified expression matrix after
#'   \code{limpa::dpcQuant()} / \code{limpa::dpcQuantByRow()} and before
#'   \code{limpa::dpcDE()}. One of \code{"none"} (default) or \code{"quantile"}.
#'   Avoid double-normalization if your input was already normalized upstream.
#' @param limpaDEArgs For \code{type="limpa"}, optional named list of additional
#'   arguments forwarded to \code{limpa::dpcDE()} (and then to
#'   \code{limpa::voomaLmFitWithImputation()}). To enable limma sample weights,
#'   set \code{limpaDEArgs = list(sample.weights = TRUE)}.
#' @param limpaKeep For \code{type="limpa"}, optionally store intermediate limpa
#'   objects in \code{backendData(nst)$limpa} for reuse/debugging. Set to
#'   \code{"elist"} to keep the \code{EList} object (completed expression matrix
#'   plus uncertainty estimates), \code{"fit"} to keep the fitted \code{MArrayLM}
#'   object(s), or \code{"all"} to keep both. Default is \code{"none"}.
#' @param quiet Omit status messages printed during run
#'
#' @param sigThres Significance threshold use for illustrating significant hits
#'   in diagnostic plots
#' @param sigThresType Type of significance threshold, "fdr" or "p". "fdr" is
#'   strongly recommended (Benjamini-Hochberg corrected p-values)
#' @param log2FoldThres Fold-size cutoff for being considered significant in
#'   diagnostic plots
#' @param writeReportAsPngs Output report as separate PNG files instead of a
#'   single PDF
#' @param inputFormat Type of input format for \code{dataPath} when reading from
#'   files. Supports \code{"default"} and \code{"diann"}.
#' @param inputOptions Optional list of input-reader options. For DIA-NN, use
#'   \code{\link{diannInputOptions}}. For delimited inputs, use
#'   \code{\link{defaultInputOptions}}, \code{\link{proteiosInputOptions}}, or
#'   \code{\link{maxQuantInputOptions}} to configure the delimiter.
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
normalyzerDE <- function(
  jobName,
  comparisons = NULL,
  designPath = NULL,
  dataPath = NULL,
  experimentObj = NULL,
  outputDir = ".",
  logTrans = FALSE,
  type = "limma",
  sampleCol = "sample",
  condCol = "group",
  batchCol = NULL,
  techRepCol = NULL,
  leastRepCount = 1,
  impute = FALSE,
  imputeMinFraction = 0.75,
  subsetByComparison = FALSE,
  oneVsRest = FALSE,
  oneVsRestGroups = NULL,
  quiet = FALSE,
  sigThres = 0.1,
  sigThresType = "fdr",
  log2FoldThres = 0,
  writeReportAsPngs = FALSE,
  limpaProteinIdCol = "auto",
  limpaDpc = NULL,
  limpaDpcMethod = c("none", "dpc", "dpcCN"),
  limpaDpcArgs = NULL,
  limpaQuantArgs = NULL,
  limpaDEArgs = NULL,
  limpaKeep = c("none", "elist", "fit", "all"),
  inputFormat = "default",
  inputOptions = NULL,
  limpaByRow = FALSE,
  limpaPostQuantNorm = c("none", "quantile")
) {
  if (!quiet) {
    message(
      "You are running version ",
      utils::packageVersion("NormalyzerDE"),
      " of NormalyzerDE"
    )
  }

  if (is.null(experimentObj) && (is.null(designPath) || is.null(dataPath))) {
    stop(
      "Either options 'designPath' plus 'dataPath' or 'summarizedExp' need to be provided"
    )
  }

  if (!oneVsRest && is.null(comparisons)) {
    stop(
      "Argument 'comparisons' must be provided. Specify one or more comparisons as a vector.\n",
      "Example, one comparison between group 1 and 2: c('1-2')\n",
      "Example, two comparisons between groups 1 and 2, and groups 2 and 3: c('1-2', '2-3')"
    )
  }

  startTime <- Sys.time()
  jobDir <- setupJobDir(jobName, outputDir)
  safeJobName <- basename(jobDir)

  if (!quiet) {
    message("Setting up statistics object")
  }
  if (is.null(experimentObj)) {
    experimentObj <- setupRawContrastObject(
      dataPath,
      designPath,
      sampleCol,
      inputFormat = inputFormat,
      inputOptions = inputOptions
    )
  } else {
    verifySummarizedExperiment(experimentObj, sampleCol)
  }

  if (!is.null(techRepCol)) {
    if (!quiet) {
      message("Reducing technical replicates")
    }
    experimentObj <- reduceTechnicalReplicates(
      experimentObj,
      techRepCol,
      sampleCol
    )
  }

  nst <- NormalyzerStatistics(
    experimentObj,
    logTrans = logTrans
  )

  if (!quiet) {
    message("Calculating statistical contrasts...")
  }
  nst <- calculateContrasts(
    nst,
    comparisons,
    type = type,
    condCol = condCol,
    batchCol = batchCol,
    leastRepCount = leastRepCount,
    impute = impute,
    imputeMinFraction = imputeMinFraction,
    subsetByComparison = subsetByComparison,
    oneVsRest = oneVsRest,
    oneVsRestGroups = oneVsRestGroups,
    limpaProteinIdCol = limpaProteinIdCol,
    limpaByRow = limpaByRow,
    limpaDpc = limpaDpc,
    limpaDpcMethod = limpaDpcMethod,
    limpaDpcArgs = limpaDpcArgs,
    limpaQuantArgs = limpaQuantArgs,
    limpaPostQuantNorm = limpaPostQuantNorm,
    limpaDEArgs = limpaDEArgs,
    limpaKeep = limpaKeep
  )
  if (!quiet) {
    message("Contrast calculations done!")
  }

  annotDf <- generateAnnotatedMatrix(nst)
  outPath <- paste0(jobDir, "/", safeJobName, "_stats.tsv")

  if (!quiet) {
    message("Writing ", nrow(annotDf), " annotated rows to ", outPath)
  }
  utils::write.table(
    annotDf,
    file = outPath,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  if (!quiet) {
    message("Writing statistics report")
  }
  generateStatsReport(
    nst,
    safeJobName,
    jobDir,
    sigThres,
    sigThresType,
    log2FoldThres,
    writeAsPngs = writeReportAsPngs
  )

  endTime <- Sys.time()
  totTime <- difftime(endTime, startTime, units = "mins")
  if (!quiet) {
    message(
      "All done! Results are stored in: ",
      jobDir,
      ", processing time was ",
      round(totTime, 1),
      " minutes"
    )
  }
}
