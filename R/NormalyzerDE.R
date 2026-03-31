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
#' @param reuseOutputDir Reuse an existing non-empty output directory. By default,
#'   NormalyzerDE errors to avoid mixing outputs from different runs.
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
#' @param noLogTransform Don't log2-transform the input.
#' @param writeReportAsPngs Write the evaluation report as separate PNG files
#'   instead of a single PDF
#'
#' @param rtStepSizeMinutes Retention time normalization window size.
#' @param rtWindowMinCount Minimum number of datapoints in each retention-time
#'   segment.
#' @param rtWindowShifts Number of layered retention time normalized windows.
#' @param rtWindowMergeMethod Merge approach for layered retention time windows.
#'
#' @param limpaOptions Optional helper created by \code{\link{limpaOptions}}.
#'   For \code{preQuant = "limpa"}, use this to configure \pkg{limpa}
#'   quantification settings such as \code{proteinIdCol}, \code{byRow},
#'   \code{dpc}, \code{dpcMethod}, \code{dpcArgs}, and \code{quantArgs}.
#'   See \code{\link{limpaOptions}} for the available fields.
#' @param preQuant Optional pre-quantification step applied before running the
#'   Normalyzer normalization evaluation. Use \code{"limpa"} to first complete
#'   the data matrix with \code{limpa::dpcQuant()} / \code{limpa::dpcQuantByRow()}
#'   and then evaluate Normalyzer normalizations on the post-quant log2 matrix.
#'   This follows the recommended way to use limpa (quantify first, then
#'   normalize). For \code{preQuant="limpa"} with \code{inputFormat="diann"},
#'   ambiguous DIA-NN reports default to precursor-level input unless
#'   \code{inputOptions} explicitly sets the level or columns. When enabled, the
#'   quantified \code{EList} is saved as
#'   an RDS file \code{<jobDir>/<basename(jobDir)>_limpa_quantified.rds} for reuse with
#'   \code{\link{normalyzerDE}} via \code{limpaOptions(quantifiedRds = ...)}.
#'   Reuse is explicit: \code{\link{normalyzerDE}} does not pick up nearby
#'   quantified caches automatically. The saved \code{EList} can also be
#'   filtered after quantification using
#'   \code{y$other$n.observations} and normalized in a study-specific way before
#'   differential testing.
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
  reuseOutputDir = FALSE,
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
  rtWindowMergeMethod = "mean",

  preQuant = c("none", "limpa"),
  limpaOptions = NULL
) {
  if (!quiet) {
    version <- utils::packageVersion("NormalyzerDE")
    cli::cli_alert_info("You are running version {version} of NormalyzerDE")
  }

  if (is.null(experimentObj) && (is.null(designPath) || is.null(dataPath))) {
    cli::cli_abort(
      "Provide {.arg designPath} + {.arg dataPath}, or provide {.arg experimentObj}.",
      class = "normalyzerde_error",
      call = NULL
    )
  }

  startTime <- Sys.time()
  preQuantUse <- match.arg(preQuant)
  totalSteps <- if (identical(preQuantUse, "limpa")) 6L else 5L

  stepTag <- function(step) {
    sprintf("[Step %d/%d]", step, totalSteps)
  }

  if (!quiet) {
    cli::cli_alert_info("{.strong {stepTag(1)}} Load data and verify input")
  }

  if (is.null(experimentObj)) {
    inputOptionsUse <- inputOptions
    if (identical(preQuantUse, "limpa") && identical(inputFormat, "diann")) {
      inputOptionsResolved <- diannPreferPrecursorInputOptionsForLimpa(
        inputOptionsUse
      )
      if (
        !quiet &&
          diannWasDefaultedToPrecursorForLimpa(
            inputOptionsUse,
            inputOptionsResolved
          )
      ) {
        cli::cli_alert_info(
          "For {.arg inputFormat}={.val diann} with {.arg preQuant}={.val limpa}, defaulting to precursor-level DIA-NN columns because {.arg inputOptions} did not specify the level or columns explicitly."
        )
      }
      inputOptionsUse <- inputOptionsResolved
      oldWarn <- getOption("NormalyzerDE.warnDiannAutoAmbiguous")
      options(NormalyzerDE.warnDiannAutoAmbiguous = TRUE)
      on.exit(
        options(NormalyzerDE.warnDiannAutoAmbiguous = oldWarn),
        add = TRUE
      )
    }
    experimentObj <- setupRawDataObject(
      dataPath = dataPath,
      designPath = designPath,
      inputFormat = inputFormat,
      zeroToNA = zeroToNA,
      sampleColName = sampleColName,
      groupColName = groupColName,
      inputOptions = inputOptionsUse
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

  jobDir <- setupJobDir(jobName, outputDir, reuseOutputDir = reuseOutputDir)
  if (!quiet) {
    cli::cli_alert_success(
      "{.strong {stepTag(1)}} Input verified; output directory prepared at {.path {jobDir}}"
    )
  }

  noLogTransformUse <- noLogTransform
  if (identical(preQuantUse, "limpa")) {
    limpaOptions <- resolveLimpaOptions(limpaOptions)
    limpaProteinIdCol <- limpaOptions$proteinIdCol
    limpaByRow <- limpaOptions$byRow
    limpaDpc <- limpaOptions$dpc
    limpaDpcMethod <- limpaOptions$dpcMethod
    limpaDpcArgs <- limpaOptions$dpcArgs
    limpaQuantArgs <- limpaOptions$quantArgs

    if (!quiet) {
      cli::cli_alert_info(
        "{.strong {stepTag(2)}} Running limpa pre-quantification"
      )
    }

    requireLimpaPackageInternal("preQuant='limpa'")
    limpaQuantByRow <- resolveLimpaQuantByRowFn()

    limpaDpcMethod <- match.arg(
      limpaDpcMethod,
      c("none", "dpc", "dpcON", "dpcCN")
    )
    validateLimpaDpc(limpaDpc)

    dpcArgsUse <- sanitizeLimpaDpcArgs(limpaDpcArgs)
    quantArgsUse <- applyLimpaQuantDefaultsAndValidate(
      sanitizeLimpaQuantArgs(limpaQuantArgs)
    )

    byRowConfig <- normalizeLimpaByRowConfig(
      limpaByRow = limpaByRow,
      limpaProteinIdCol = limpaProteinIdCol
    )
    limpaByRowUse <- byRowConfig$limpaByRow
    limpaProteinIdCol <- byRowConfig$limpaProteinIdCol

    verboseUse <- isTRUE(quantArgsUse[["verbose"]])
    dpcSlopeUse <- as.numeric(quantArgsUse[["dpc.slope"]])

    log2WithNonFiniteAsNA <- function(mat) {
      wasMissing <- is.na(mat)
      out <- log2(mat)
      nonFinite <- !is.finite(out) & !wasMissing
      if (any(nonFinite)) {
        cli::cli_warn(
          "Non-finite values produced by log2 transform (e.g. zeros or negative values) were treated as missing (set to NA).",
          class = "normalyzerde_warning",
          call = NULL
        )
        out[nonFinite] <- NA_real_
      }
      out
    }

    log2Mat <- filterrawdata(normObj)
    if (!noLogTransform) {
      log2Mat <- log2WithNonFiniteAsNA(log2Mat)
    }

    annotMatRaw <- annotationValues(normObj)
    keepRows <- rowSums(!is.na(log2Mat)) > 0
    log2Mat <- log2Mat[keepRows, , drop = FALSE]
    annotMatRaw <- annotMatRaw[keepRows, , drop = FALSE]

    genesInputAll <- as.data.frame(
      annotMatRaw,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    limpaProteinIdColUsed <- inferLimpaProteinIdCol(
      annotMatRaw,
      limpaProteinIdCol
    )

    limpaDpcUse <- if (!is.null(limpaDpc)) {
      if (!quiet && limpaDpcMethod != "none" && isTRUE(verboseUse)) {
        cli::cli_alert_info(
          "{.arg dpc} was supplied; ignoring {.arg dpcMethod}={.val {limpaDpcMethod}}."
        )
      }
      limpaDpc
    } else {
      estimateLimpaDpcFromData(
        dataMat = log2Mat,
        limpaDpcMethod = limpaDpcMethod,
        limpaDpcArgs = dpcArgsUse,
        dpcSlope = dpcSlopeUse,
        verbose = verboseUse
      )
    }

    quantifyByRow <- function(mat, genesDf) {
      quantifyLimpaByRowInternal(
        dataMat = mat,
        genesDf = genesDf,
        dpc = limpaDpcUse,
        quantArgs = quantArgsUse,
        limpaQuantByRow = limpaQuantByRow
      )
    }

    quantifyByProtein <- function(mat, genesDf, proteinIdCol) {
      quantifyLimpaByProteinInternal(
        dataMat = mat,
        genesDf = genesDf,
        proteinIdCol = proteinIdCol,
        dpc = limpaDpcUse,
        quantArgs = quantArgsUse
      )
    }

    yQuant <- NULL
    if (!is.null(limpaProteinIdColUsed)) {
      proteinIdVec <- as.character(genesInputAll[[limpaProteinIdColUsed]])
      if (anyDuplicated(proteinIdVec) > 0) {
        yQuant <- quantifyByProtein(
          log2Mat,
          genesDf = genesInputAll,
          proteinIdCol = limpaProteinIdColUsed
        )
      } else {
        if (!quiet) {
          cli::cli_alert_info(
            "{.arg proteinIdCol} {.val {limpaProteinIdColUsed}} contains no duplicated identifiers; skipping peptide/precursor-to-protein summarization."
          )
        }
        yQuant <- quantifyByRow(log2Mat, genesDf = genesInputAll)
      }
    } else {
      yQuant <- quantifyByRow(log2Mat, genesDf = genesInputAll)
    }

    yQuant <- normalizeLimpaQuantifiedEList(yQuant)

    safeJobName <- basename(jobDir)
    quantifiedRds <- file.path(
      jobDir,
      paste0(safeJobName, "_limpa_quantified.rds")
    )
    saveRDS(yQuant, file = quantifiedRds)

    if (!quiet) {
      cli::cli_alert_info(
        "{.strong {stepTag(2)}} Saved quantified EList to {.path {quantifiedRds}}"
      )
    }

    noLogTransformUse <- TRUE
    designDf <- designMatrix(normObj)
    sampleColUsed <- sampleNameCol(normObj)
    groupColUsed <- groupNameCol(normObj)

    sampleNamesUse <- as.character(designDf[[sampleColUsed]])
    postQuantMat <- as.matrix(yQuant$E)[, sampleNamesUse, drop = FALSE]
    annotMatQuant <- if (!is.null(yQuant$genes)) {
      as.matrix(yQuant$genes)
    } else {
      matrix(nrow = nrow(postQuantMat), ncol = 0)
    }

    normObj <- NormalyzerDataset(
      jobName = jobName,
      designMatrix = designDf,
      rawData = postQuantMat,
      annotationData = annotMatQuant,
      sampleNameCol = sampleColUsed,
      groupNameCol = groupColUsed,
      tinyRunThres = tinyRunThres,
      quiet = quiet
    )

    if (!quiet) {
      cli::cli_alert_success(
        "{.strong {stepTag(2)}} limpa quantification completed"
      )
    }
  }

  stepLabel <- if (identical(preQuantUse, "limpa")) 3 else 2
  if (!quiet) {
    cli::cli_alert_info(
      "{.strong {stepTag(stepLabel)}} Performing normalizations"
    )
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
    noLogTransform = noLogTransformUse
  )
  if (!quiet) {
    cli::cli_alert_success(
      "{.strong {stepTag(stepLabel)}} Normalizations completed"
    )
  }

  analysisStepLabel <- if (identical(preQuantUse, "limpa")) 4 else 3
  if (!skipAnalysis) {
    if (!quiet) {
      cli::cli_alert_info(
        "{.strong {stepTag(analysisStepLabel)}} Generating evaluation measures..."
      )
    }
    normalyzerResultsObject <- analyzeNormalizations(normalyzerResultsObject)
    if (!quiet) {
      cli::cli_alert_success(
        "{.strong {stepTag(analysisStepLabel)}} Evaluation measures generated"
      )
    }
  } else {
    if (!quiet) {
      cli::cli_alert_success(
        "{.strong {stepTag(analysisStepLabel)}} Skipped evaluation measures (skipAnalysis=TRUE)"
      )
    }
  }

  writeStepLabel <- if (identical(preQuantUse, "limpa")) 5 else 4
  if (!quiet) {
    cli::cli_alert_info(
      "{.strong {stepTag(writeStepLabel)}} Writing matrices to file"
    )
  }
  writeNormalizedDatasets(normalyzerResultsObject, jobDir)
  if (!quiet) {
    cli::cli_alert_success(
      "{.strong {stepTag(writeStepLabel)}} Matrices successfully written"
    )
  }

  plotStepLabel <- if (identical(preQuantUse, "limpa")) 6 else 5
  if (!skipAnalysis) {
    if (!quiet) {
      cli::cli_alert_info(
        "{.strong {stepTag(plotStepLabel)}} Generating plots..."
      )
    }
    generatePlots(
      normalyzerResultsObject,
      jobDir,
      plotRows = plotRows,
      plotCols = plotCols,
      writeAsPngs = writeReportAsPngs
    )
    if (!quiet) {
      cli::cli_alert_success(
        "{.strong {stepTag(plotStepLabel)}} Plots successfully generated"
      )
    }
  } else {
    if (!quiet) {
      cli::cli_alert_success(
        "{.strong {stepTag(plotStepLabel)}} Skipped plot generation (skipAnalysis=TRUE)"
      )
    }
  }

  endTime <- Sys.time()
  totTime <- difftime(endTime, startTime, units = "mins")
  if (!quiet) {
    cli::cli_alert_success(
      "All done! Results are saved in {.path {jobDir}}; processing time was {round(totTime, 1)} minutes"
    )
  }

  invisible(NULL)
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
#' \code{limpaOptions(postQuantNorm = ...)}. Whether normalization is appropriate is
#' dataset-dependent; for studies with expected global shifts between groups,
#' extra normalization can remove the biological signal of interest. Avoid
#' applying multiple normalizations unintentionally. For \code{type="limpa"}
#' with \code{inputFormat="diann"},
#' ambiguous DIA-NN reports default to precursor-level input unless
#' \code{inputOptions} explicitly sets the level or columns.
#' In \code{limpaOptions()}, \code{postQuantNorm} accepts \code{"none"},
#' \code{"GI"}, \code{"median"}, \code{"mean"}, \code{"Quantile"} (or
#' \code{"quantile"}), \code{"CycLoess"}, and \code{"RLR"}.
#'
#' A common \pkg{limpa} workflow is to keep all samples through
#' \code{dpcQuant()}, then filter sparse proteins using
#' \code{y$other$n.observations} before differential testing (or, within
#' NormalyzerDE, use \code{leastRepCount} or pass a filtered quantified object
#' via \code{limpaOptions(quantifiedRds = ...)}). Reuse is explicit:
#' \code{\link{normalyzerDE}} does not auto-detect nearby quantified caches.
#' Outlier samples are often
#' better assessed with sample-specific
#' weights and QC/MDS plots than with density plots alone. To estimate sample
#' weights in NormalyzerDE, set
#' \code{limpaOptions(deArgs = list(sample.weights = TRUE))} and extract them with
#' \code{\link{getLimpaSampleWeights}} or from the sample-weights TSV written by
#' \code{\link{normalyzerDE}}.
#'
#' For PTM-level data (e.g., phosphoproteomics) where each row corresponds to a
#' modified site, set \code{limpaOptions(byRow = TRUE)} (or
#' \code{limpaOptions(proteinIdCol = NULL)}) to keep each row separate rather
#' than summarizing to protein-level.
#'
#' By default, NormalyzerDE uses a fixed DPC slope
#' (\code{limpaOptions(quantArgs = list(dpc.slope = 0.8))}) and lets limpa
#' estimate the intercept. To estimate both DPC parameters from your data, set
#' \code{limpaOptions(dpcMethod = "dpc")} (or \code{"dpcON"} /
#' \code{"dpcCN"} for more robust estimates, especially in datasets with very
#' large fold-changes).
#'
#' @param jobName Name of job
#' @param designPath File path to design matrix
#' @param dataPath File path to normalized matrix or completed log2 matrix.
#' @param experimentObj SummarizedExperiment object, can be provided as input
#'  as alternative to 'designPath' and 'dataPath'
#' @param comparisons Character vector containing target contrasts.
#'   If comparing condA with condB, then the vector would be c("condA-condB").
#'   Ignored if \code{oneVsRest=TRUE}.
#' @param outputDir Path to output directory
#' @param reuseOutputDir Reuse an existing non-empty output directory. By default,
#'   NormalyzerDE errors to avoid mixing outputs from different runs.
#' @param logTrans Log2-transform the input (needed if providing non-log2
#'   input).
#' @param type Type of statistical comparison, "limma", "limma_intensity" or
#'  "welch" or "limpa", where "limma_intensity" allows the prior to be fit
#'  according to intensity rather than using a flat prior. "limpa" uses the
#'  Bioconductor package \pkg{limpa} to handle missing values via a
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
#'   before NA-filtering, imputation and model fitting. Use this when filtering
#'   or imputation should depend only on the samples in each comparison, not the
#'   full dataset.
#' @param oneVsRest If TRUE, create one comparison per selected group against
#'   all remaining samples (for all groups in \code{condCol}, or the subset in
#'   \code{oneVsRestGroups}).
#' @param oneVsRestGroups Optional character vector specifying which groups in
#'   \code{condCol} to compare against all other samples.
#' @param limpaOptions Optional helper created by \code{\link{limpaOptions}}.
#'   For \code{type = "limpa"}, use this to configure protein-level
#'   summarization (\code{proteinIdCol}, \code{byRow}), DPC estimation
#'   (\code{dpc}, \code{dpcMethod}, \code{dpcArgs}), quantification
#'   (\code{quantArgs}, \code{quantifiedRds}), differential expression
#'   (\code{deArgs}), optional post-quantification normalization
#'   (\code{postQuantNorm}), and retained backend objects (\code{keep}).
#'   See \code{\link{limpaOptions}} for the available fields.
#' @param quiet Omit status messages printed during run
#'
#' @param sigThres Significance threshold use for illustrating significant hits
#'   in diagnostic plots
#' @param sigThresType Type of significance threshold, "fdr" or "p". "fdr" is
#'   strongly recommended (Benjamini-Hochberg corrected p-values)
#' @param log2FoldThres Fold-size cutoff for being considered significant in
#'   diagnostic plots
#' @param writeReportAsPngs Write the report as separate PNG files instead of a
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
  reuseOutputDir = FALSE,
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
  limpaOptions = NULL,
  inputFormat = "default",
  inputOptions = NULL
) {
  if (!quiet) {
    version <- utils::packageVersion("NormalyzerDE")
    cli::cli_alert_info("You are running version {version} of NormalyzerDE")
  }

  if (is.null(experimentObj) && (is.null(designPath) || is.null(dataPath))) {
    cli::cli_abort(
      "Provide {.arg designPath} + {.arg dataPath}, or provide {.arg experimentObj}.",
      class = "normalyzerde_error",
      call = NULL
    )
  }

  isLimpa <- identical(as.character(type)[1], "limpa")

  limpaOptions <- if (isLimpa) {
    resolveLimpaOptions(limpaOptions)
  } else {
    limpaOptions
  }
  limpaProteinIdCol <- if (isLimpa) {
    limpaOptions$proteinIdCol
  } else {
    NULL
  }
  limpaByRow <- if (isLimpa) {
    limpaOptions$byRow
  } else {
    FALSE
  }
  limpaDpc <- if (isLimpa) {
    limpaOptions$dpc
  } else {
    NULL
  }
  limpaDpcMethod <- if (isLimpa) {
    limpaOptions$dpcMethod
  } else {
    "none"
  }
  limpaDpcArgs <- if (isLimpa) {
    limpaOptions$dpcArgs
  } else {
    list()
  }
  limpaQuantArgs <- if (isLimpa) {
    limpaOptions$quantArgs
  } else {
    list()
  }
  limpaQuantifiedRds <- if (isLimpa) {
    limpaOptions$quantifiedRds
  } else {
    NULL
  }
  limpaDEArgs <- if (isLimpa) {
    limpaOptions$deArgs
  } else {
    list()
  }
  limpaKeep <- if (isLimpa) {
    limpaOptions$keep
  } else {
    "none"
  }
  limpaPostQuantNorm <- if (isLimpa) {
    limpaOptions$postQuantNorm
  } else {
    "none"
  }

  if (!oneVsRest && is.null(comparisons)) {
    cli::cli_abort(
      c(
        "Provide {.arg comparisons}, or set {.arg oneVsRest}={.val TRUE}.",
        i = "Example (one comparison): {.code comparisons = c('1-2')}",
        i = "Example (two comparisons): {.code comparisons = c('1-2', '2-3')}"
      ),
      class = "normalyzerde_error",
      call = NULL
    )
  }

  if (isLimpa) {
    limpaOptions$quantifiedRds <- limpaQuantifiedRds
  }

  if (isLimpa) {
    limpaPostQuantNormUse <- match.arg(
      limpaPostQuantNorm,
      c(
        "none",
        "GI",
        "median",
        "mean",
        "Quantile",
        "CycLoess",
        "RLR",
        "quantile"
      )
    )

    if (
      !identical(limpaPostQuantNormUse, "none") &&
        !is.null(dataPath) &&
        is.character(dataPath) &&
        length(dataPath) == 1 &&
        !is.na(dataPath) &&
        nzchar(dataPath)
    ) {
      base <- basename(dataPath)
      normPattern <- "-normalized\\.(txt|tsv)$"
      if (grepl(normPattern, base, ignore.case = TRUE)) {
        normPrefix <- sub(normPattern, "", base, ignore.case = TRUE)
        prefixLower <- tolower(normPrefix)

        if (!identical(prefixLower, "log2")) {
          normLower <- tolower(limpaPostQuantNormUse)
          if (identical(prefixLower, normLower)) {
            cli::cli_abort(
              c(
                "Input file {.path {base}} appears to be already normalized ({.val {normPrefix}}). Applying {.arg postQuantNorm}={.val {limpaPostQuantNormUse}} would repeat the same normalization step.",
                i = "Use {.arg postQuantNorm}={.val none} when supplying a pre-normalized matrix.",
                i = "Or use the completed log2 matrix ({.file log2-normalized.txt}) and set {.arg postQuantNorm} to the desired method."
              ),
              class = "normalyzerde_error",
              call = NULL
            )
          }

          cli::cli_warn(
            c(
              "Input file {.path {base}} appears to be already normalized ({.val {normPrefix}}). Applying {.arg postQuantNorm}={.val {limpaPostQuantNormUse}} would add another normalization step.",
              i = "This may be unintended double-normalization."
            ),
            class = "normalyzerde_warning",
            call = NULL
          )
        }
      }
    }
  }

  startTime <- Sys.time()
  totalSteps <- 6L

  stepTag <- function(step) {
    sprintf("[Step %d/%d]", step, totalSteps)
  }

  if (!quiet) {
    cli::cli_alert_info("{.strong {stepTag(1)}} Load data and verify input")
  }
  jobDir <- setupJobDir(jobName, outputDir, reuseOutputDir = reuseOutputDir)
  safeJobName <- basename(jobDir)

  if (is.null(experimentObj)) {
    inputOptionsUse <- inputOptions
    if (isLimpa && identical(inputFormat, "diann")) {
      inputOptionsResolved <- diannPreferPrecursorInputOptionsForLimpa(
        inputOptionsUse
      )
      if (
        !quiet &&
          diannWasDefaultedToPrecursorForLimpa(
            inputOptionsUse,
            inputOptionsResolved
          )
      ) {
        cli::cli_alert_info(
          "For {.arg inputFormat}={.val diann} with {.arg type}={.val limpa}, defaulting to precursor-level DIA-NN columns because {.arg inputOptions} did not specify the level or columns explicitly."
        )
      }
      inputOptionsUse <- inputOptionsResolved
      oldWarn <- getOption("NormalyzerDE.warnDiannAutoAmbiguous")
      options(NormalyzerDE.warnDiannAutoAmbiguous = TRUE)
      on.exit(
        options(NormalyzerDE.warnDiannAutoAmbiguous = oldWarn),
        add = TRUE
      )
    }
    experimentObj <- setupRawContrastObject(
      dataPath,
      designPath,
      sampleCol,
      inputFormat = inputFormat,
      inputOptions = inputOptionsUse
    )
  } else {
    verifySummarizedExperiment(experimentObj, sampleCol)
  }
  if (!quiet) {
    cli::cli_alert_success(
      "{.strong {stepTag(1)}} Input verified; output directory prepared at {.path {jobDir}}"
    )
  }

  if (is.null(techRepCol)) {
    if (!quiet) {
      cli::cli_alert_success(
        "{.strong {stepTag(2)}} Skipped technical replicate reduction (techRepCol is NULL)"
      )
    }
  } else {
    if (!quiet) {
      cli::cli_alert_info(
        "{.strong {stepTag(2)}} Reducing technical replicates"
      )
    }
    experimentObj <- reduceTechnicalReplicates(
      experimentObj,
      techRepCol,
      sampleCol
    )
    if (!quiet) {
      cli::cli_alert_success(
        "{.strong {stepTag(2)}} Technical replicates reduced"
      )
    }
  }

  if (!quiet) {
    cli::cli_alert_info("{.strong {stepTag(3)}} Setting up statistics object")
  }
  nst <- NormalyzerStatistics(
    experimentObj,
    logTrans = logTrans
  )
  if (!quiet) {
    cli::cli_alert_success(
      "{.strong {stepTag(3)}} Statistics object prepared"
    )
  }

  if (!quiet) {
    cli::cli_alert_info(
      "{.strong {stepTag(4)}} Calculating statistical contrasts..."
    )
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
    limpaOptions = limpaOptions
  )
  if (!quiet) {
    cli::cli_alert_success(
      "{.strong {stepTag(4)}} Contrast calculations done!"
    )
  }

  annotDf <- generateAnnotatedMatrix(nst)
  outPath <- paste0(jobDir, "/", safeJobName, "_stats.tsv")

  if (!quiet) {
    cli::cli_alert_info(
      "{.strong {stepTag(5)}} Writing {nrow(annotDf)} annotated rows to {.path {outPath}}"
    )
  }
  utils::write.table(
    annotDf,
    file = outPath,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  if (!quiet) {
    cli::cli_alert_success("{.strong {stepTag(5)}} Annotated matrix written")
  }

  sampleWeightsDf <- collectLimpaSampleWeights(nst)
  if (!is.null(sampleWeightsDf) && nrow(sampleWeightsDf) > 0) {
    sampleWeightsPath <- paste0(jobDir, "/", safeJobName, "_sample_weights.tsv")
    if (!quiet) {
      cli::cli_alert_info(
        "{.strong {stepTag(5)}} Writing {nrow(sampleWeightsDf)} limpa sample-weight rows to {.path {sampleWeightsPath}}"
      )
    }
    writeLimpaSampleWeights(sampleWeightsDf, sampleWeightsPath)
    if (!quiet) {
      cli::cli_alert_success(
        "{.strong {stepTag(5)}} Sample weights written"
      )
    }
  }

  if (!quiet) {
    cli::cli_alert_info("{.strong {stepTag(6)}} Writing statistics report")
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
    cli::cli_alert_success("{.strong {stepTag(6)}} Statistics report written")
    cli::cli_alert_success(
      "All done! Results are saved in {.path {jobDir}}; processing time was {round(totTime, 1)} minutes"
    )
  }

  invisible(NULL)
}
