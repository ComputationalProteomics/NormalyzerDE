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
#' @param preQuant Optional pre-quantification step applied before running the
#'   Normalyzer normalization evaluation. Use \code{"limpa"} to first complete
#'   the data matrix with \code{limpa::dpcQuant()} / \code{limpa::dpcQuantByRow()}
#'   and then evaluate Normalyzer normalizations on the post-quant log2 matrix.
#'   This follows the recommended way to use limpa (quantify first, then
#'   normalize). When enabled, the quantified \code{EList} is saved as
#'   an RDS file \code{<jobDir>/<basename(jobDir)>_limpa_quantified.rds} for reuse with
#'   \code{\link{normalyzerDE}} via \code{limpaQuantifiedRds}.
#' @param limpaProteinIdCol For \code{preQuant="limpa"}, optionally summarize
#'   peptide/precursor rows to protein-level using \code{limpa::dpcQuant()}.
#'   Set to a column name in the row annotation (for example
#'   \code{"Protein.Group"}) to use as the protein identifier. Use
#'   \code{"auto"} (default) to try common identifiers. If the chosen column
#'   contains duplicate identifiers, the data are summarized once across all
#'   samples and the output rows correspond to proteins. Set to \code{NULL} to
#'   disable protein summarization and quantify each row separately as one feature via
#'   \code{limpa::dpcQuantByRow()} (recommended for PTM-level data such as
#'   phosphoproteomics where each row corresponds to a modified site).
#' @param limpaByRow For \code{preQuant="limpa"}, treat each input row as a
#'   separate feature and always use \code{limpa::dpcQuantByRow()} instead of
#'   summarizing via \code{limpa::dpcQuant()}. This is recommended for PTM-level
#'   matrices (e.g., phosphosites). Equivalent to setting
#'   \code{limpaProteinIdCol=NULL}.
#' @param limpaDpc For \code{preQuant="limpa"}, optional DPC parameters to pass
#'   to \code{limpa::dpcQuant()} / \code{limpa::dpcQuantByRow()}. Can be a list as
#'   returned by \code{limpa::dpc()}/\code{limpa::dpcON()}/\code{limpa::dpcCN()},
#'   or a numeric vector \code{c(beta0, beta1)}.
#' @param limpaDpcMethod For \code{preQuant="limpa"}, optional method to estimate
#'   the DPC parameters from the data when \code{limpaDpc} is not supplied.
#'   \code{"none"} (default) uses a fixed slope (\code{limpaQuantArgs$dpc.slope},
#'   default \code{0.8}) and lets limpa estimate the intercept internally.
#'   \code{"dpc"} estimates both DPC parameters from the observed-normal model
#'   via \code{limpa::dpc()}. \code{"dpcON"} estimates the DPC from the
#'   observed-normal model via the newer \code{limpa::dpcON()} (for a robust fit,
#'   set \code{limpaDpcArgs=list(robust=TRUE)}). \code{"dpcCN"} estimates the DPC from the
#'   complete-normal model via \code{limpa::dpcCN()}, which can be more robust
#'   for datasets with very large fold-changes.
#' @param limpaDpcArgs For \code{preQuant="limpa"}, optional named list of
#'   additional arguments forwarded to \code{limpa::dpc()}, \code{limpa::dpcON()},
#'   or \code{limpa::dpcCN()}
#'   when \code{limpaDpcMethod} is not \code{"none"}. Argument \code{y} is ignored.
#'   For \code{limpaDpcMethod="dpcON"} and \code{"dpcCN"}, \code{dpc.slope.start} defaults to
#'   \code{limpaQuantArgs$dpc.slope}.
#' @param limpaQuantArgs For \code{preQuant="limpa"}, optional named list of
#'   additional arguments forwarded to \code{limpa::dpcQuant()} /
#'   \code{limpa::dpcQuantByRow()}. Use this to set \code{dpc.slope} (default
#'   \code{0.8}), \code{chunk} (default \code{1000L}), and \code{verbose} (default
#'   \code{FALSE}), plus any additional \code{...} arguments supported by limpa.
#'   Arguments \code{y}, \code{protein.id}, and \code{dpc} are ignored.
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
  rtWindowMergeMethod = "mean",

  preQuant = c("none", "limpa"),
  limpaProteinIdCol = "auto",
  limpaByRow = FALSE,
  limpaDpc = NULL,
  limpaDpcMethod = c("none", "dpc", "dpcON", "dpcCN"),
  limpaDpcArgs = NULL,
  limpaQuantArgs = NULL
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
  preQuantUse <- match.arg(preQuant)
  totalSteps <- if (identical(preQuantUse, "limpa")) 6 else 5

  if (!quiet) {
    message("[Step 1/", totalSteps, "] Load data and verify input")
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
      "[Step 1/",
      totalSteps,
      "] Input verified, job directory prepared at:",
      jobDir
    )
  }

  noLogTransformUse <- noLogTransform
  if (identical(preQuantUse, "limpa")) {
    if (!quiet) {
      message("[Step 2/", totalSteps, "] Running limpa pre-quantification")
    }

    requireLimpaPackageInternal("preQuant='limpa'")

    limpaDpcMethod <- match.arg(limpaDpcMethod)
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
        warning(
          "Non-finite values produced by log2 transform (e.g. zeros or negative values) ",
          "were treated as missing (set to NA)."
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
      if (limpaDpcMethod != "none" && isTRUE(verboseUse)) {
        message(
          "limpaDpc was supplied; ignoring limpaDpcMethod='",
          limpaDpcMethod,
          "'."
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
      y <- methods::new("EList", list(E = mat, genes = genesDf))
      do.call(
        limpa::dpcQuantByRow,
        c(list(y = y, dpc = limpaDpcUse), quantArgsUse)
      )
    }

    quantifyByProtein <- function(mat, genesDf, proteinIdCol) {
      proteinId <- genesDf[[proteinIdCol]]
      proteinId <- as.character(proteinId)

      if (length(proteinId) != nrow(mat)) {
        stop(
          "Row annotation column '",
          proteinIdCol,
          "' does not match the number of rows in the data matrix."
        )
      }

      if (anyNA(proteinId) || any(proteinId == "")) {
        stop(
          "Row annotation column '",
          proteinIdCol,
          "' contains missing or empty protein identifiers. ",
          "Remove these rows or choose another column."
        )
      }

      stableCols <- inferStableProteinAnnotationCols(
        genesDf = genesDf,
        proteinId = proteinId,
        proteinIdCol = proteinIdCol
      )
      genesForQuant <- genesDf[,
        unique(c(proteinIdCol, stableCols)),
        drop = FALSE
      ]

      yPeptide <- methods::new("EList", list(E = mat, genes = genesForQuant))
      yProtein <- do.call(
        limpa::dpcQuant,
        c(
          list(
            y = yPeptide,
            protein.id = proteinIdCol,
            dpc = limpaDpcUse
          ),
          quantArgsUse
        )
      )

      proteinIds <- rownames(yProtein$E)
      rowIds <- as.character(seq_len(nrow(yProtein$E)))

      rownames(yProtein$E) <- rowIds
      if (!is.null(yProtein$other$n.observations)) {
        rownames(yProtein$other$n.observations) <- rowIds
      }
      if (!is.null(yProtein$other$standard.error)) {
        rownames(yProtein$other$standard.error) <- rowIds
      }

      genes <- if (!is.null(yProtein$genes)) {
        as.data.frame(yProtein$genes, check.names = FALSE)
      } else {
        data.frame(check.names = FALSE)
      }
      if (!(proteinIdCol %in% colnames(genes))) {
        genes[[proteinIdCol]] <- proteinIds
      }
      genes <- genes[,
        c(proteinIdCol, setdiff(names(genes), proteinIdCol)),
        drop = FALSE
      ]
      rownames(genes) <- rowIds
      yProtein$genes <- genes

      yProtein
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
          message(
            "limpaProteinIdCol '",
            limpaProteinIdColUsed,
            "' contains no duplicated identifiers, so no peptide/precursor-to-protein summarization ",
            "will be performed."
          )
        }
        yQuant <- quantifyByRow(log2Mat, genesDf = genesInputAll)
      }
    } else {
      yQuant <- quantifyByRow(log2Mat, genesDf = genesInputAll)
    }

    rowIds <- as.character(seq_len(nrow(yQuant$E)))
    rownames(yQuant$E) <- rowIds
    if (!is.null(yQuant$other$n.observations)) {
      rownames(yQuant$other$n.observations) <- rowIds
    }
    if (!is.null(yQuant$other$standard.error)) {
      rownames(yQuant$other$standard.error) <- rowIds
    }
    if (!is.null(yQuant$genes)) {
      rownames(yQuant$genes) <- rowIds
    }

    safeJobName <- basename(jobDir)
    quantifiedRds <- file.path(
      jobDir,
      paste0(safeJobName, "_limpa_quantified.rds")
    )
    saveRDS(yQuant, file = quantifiedRds)

    if (!quiet) {
      message(
        "[Step 2/",
        totalSteps,
        "] Saved quantified EList to: ",
        quantifiedRds
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
      message("[Step 2/", totalSteps, "] Done!")
    }
  }

  if (!quiet) {
    stepLabel <- if (identical(preQuantUse, "limpa")) 3 else 2
    message("[Step ", stepLabel, "/", totalSteps, "] Performing normalizations")
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
    stepLabel <- if (identical(preQuantUse, "limpa")) 3 else 2
    message("[Step ", stepLabel, "/", totalSteps, "] Done!")
  }

  if (!skipAnalysis) {
    if (!quiet) {
      stepLabel <- if (identical(preQuantUse, "limpa")) 4 else 3
      message(
        "[Step ",
        stepLabel,
        "/",
        totalSteps,
        "] Generating evaluation measures..."
      )
    }
    normalyzerResultsObject <- analyzeNormalizations(normalyzerResultsObject)
    if (!quiet) {
      stepLabel <- if (identical(preQuantUse, "limpa")) 4 else 3
      message("[Step ", stepLabel, "/", totalSteps, "] Done!")
    }
  } else {
    if (!quiet) {
      message(
        "[Step ",
        if (identical(preQuantUse, "limpa")) 4 else 3,
        "/",
        totalSteps,
        "] skipAnalysis flag set so no analysis performed"
      )
    }
  }

  if (!quiet) {
    stepLabel <- if (identical(preQuantUse, "limpa")) 5 else 4
    message("[Step ", stepLabel, "/", totalSteps, "] Writing matrices to file")
  }
  writeNormalizedDatasets(normalyzerResultsObject, jobDir)
  if (!quiet) {
    stepLabel <- if (identical(preQuantUse, "limpa")) 5 else 4
    message(
      "[Step ",
      stepLabel,
      "/",
      totalSteps,
      "] Matrices successfully written"
    )
  }

  if (!skipAnalysis) {
    if (!quiet) {
      stepLabel <- if (identical(preQuantUse, "limpa")) 6 else 5
      message("[Step ", stepLabel, "/", totalSteps, "] Generating plots...")
    }
    generatePlots(
      normalyzerResultsObject,
      jobDir,
      plotRows = plotRows,
      plotCols = plotCols,
      writeAsPngs = writeReportAsPngs
    )
    if (!quiet) {
      stepLabel <- if (identical(preQuantUse, "limpa")) 6 else 5
      message(
        "[Step ",
        stepLabel,
        "/",
        totalSteps,
        "] Plots successfully generated"
      )
    }
  } else {
    if (!quiet) {
      message(
        "[Step ",
        if (identical(preQuantUse, "limpa")) 6 else 5,
        "/",
        totalSteps,
        "] skipAnalysis flag set so no plots generated"
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
#' \code{limpaDpcMethod="dpc"} (or \code{"dpcON"} / \code{"dpcCN"} for more robust estimates,
#' especially in datasets with very large fold-changes).
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
#'   row as one feature (recommended for PTM-level data such as phosphoproteomics
#'   where each row corresponds to a modified site).
#' @param limpaByRow For \code{type="limpa"}, treat each input row as a separate
#'   feature and always use \code{limpa::dpcQuantByRow()} instead of summarizing
#'   via \code{limpa::dpcQuant()}. This is recommended for PTM-level matrices
#'   (e.g., phosphosites). Equivalent to setting \code{limpaProteinIdCol=NULL}.
#' @param limpaDpc For \code{type="limpa"}, optional DPC parameters to pass to
#'   \code{limpa::dpcQuant()} / \code{limpa::dpcQuantByRow()}. Can be a list as
#'   returned by \code{limpa::dpc()}/\code{limpa::dpcON()}/\code{limpa::dpcCN()},
#'   or a numeric vector \code{c(beta0, beta1)}.
#' @param limpaDpcMethod For \code{type="limpa"}, optional method to estimate the
#'   DPC parameters from the data when \code{limpaDpc} is not supplied.
#'   \code{"none"} (default) uses a fixed slope (\code{limpaQuantArgs$dpc.slope},
#'   default \code{0.8}) and lets limpa estimate the intercept internally.
#'   \code{"dpc"} estimates both DPC parameters from the observed-normal model
#'   via \code{limpa::dpc()}. \code{"dpcON"} estimates the DPC from the
#'   observed-normal model via the newer \code{limpa::dpcON()} (for a robust fit,
#'   set \code{limpaDpcArgs=list(robust=TRUE)}). \code{"dpcCN"} estimates the DPC from the
#'   complete-normal model via \code{limpa::dpcCN()}, which can be more robust
#'   for datasets with very large fold-changes.
#' @param limpaDpcArgs For \code{type="limpa"}, optional named list of additional
#'   arguments forwarded to \code{limpa::dpc()}, \code{limpa::dpcON()}, or \code{limpa::dpcCN()} when
#'   \code{limpaDpcMethod} is not \code{"none"}. Argument \code{y} is ignored.
#'   For \code{limpaDpcMethod="dpcON"} and \code{"dpcCN"}, \code{dpc.slope.start} defaults to
#'   \code{limpaQuantArgs$dpc.slope}.
#' @param limpaQuantArgs For \code{type="limpa"}, optional named list of
#'   additional arguments forwarded to \code{limpa::dpcQuant()} /
#'   \code{limpa::dpcQuantByRow()}. Use this to set \code{dpc.slope} (default
#'   \code{0.8}), \code{chunk} (default \code{1000L}), and \code{verbose} (default
#'   \code{FALSE}), plus any additional \code{...} arguments supported by limpa.
#'   Arguments \code{y}, \code{protein.id}, and \code{dpc} are ignored.
#' @param limpaQuantifiedRds For \code{type="limpa"}, optional path to an RDS
#'   file containing a quantified \code{EList} object (as returned by
#'   \code{limpa::dpcQuant()} or \code{limpa::dpcQuantByRow()}). When provided,
#'   NormalyzerDE skips the internal \code{dpcQuant*()} step and reuses the
#'   cached quantification (including \code{standard.error} and
#'   \code{n.observations}) for differential expression. This is required to
#'   preserve quantification uncertainty when you first ran \code{dpcQuant*()}
#'   upstream (for example via \code{normalyzer(preQuant=\"limpa\")}).
#'   If \code{NULL}, NormalyzerDE first looks for a canonical cache file named
#'   \code{"<basename(dataDir)>_limpa_quantified.rds"} next to \code{dataPath}.
#'   If that canonical file is not present and multiple cache candidates match
#'   \code{"*_limpa_quantified.rds"}, specify \code{limpaQuantifiedRds}
#'   explicitly.
#' @param limpaPostQuantNorm For \code{type="limpa"}, optional between-sample
#'   normalization applied to the quantified expression matrix after
#'   \code{limpa::dpcQuant()} / \code{limpa::dpcQuantByRow()} and before
#'   \code{limpa::dpcDE()}. One of \code{"none"} (default), \code{"GI"},
#'   \code{"median"}, \code{"mean"}, \code{"Quantile"} (or \code{"quantile"}),
#'   \code{"CycLoess"}, or \code{"RLR"}.
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
  limpaDpcMethod = c("none", "dpc", "dpcON", "dpcCN"),
  limpaDpcArgs = NULL,
  limpaQuantArgs = NULL,
  limpaQuantifiedRds = NULL,
  limpaDEArgs = NULL,
  limpaKeep = c("none", "elist", "fit", "all"),
  inputFormat = "default",
  inputOptions = NULL,
  limpaByRow = FALSE,
  limpaPostQuantNorm = c(
    "none",
    "GI",
    "median",
    "mean",
    "Quantile",
    "CycLoess",
    "RLR",
    "quantile"
  )
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

  if (
    identical(as.character(type)[1], "limpa") && is.null(limpaQuantifiedRds)
  ) {
    limpaQuantifiedRds <- autoDetectLimpaQuantifiedRds(dataPath)
    if (!is.null(limpaQuantifiedRds) && !quiet) {
      message("Auto-detected limpaQuantifiedRds: ", limpaQuantifiedRds)
    }
  }

  if (identical(as.character(type)[1], "limpa")) {
    limpaPostQuantNormUse <- match.arg(limpaPostQuantNorm)

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
            stop(
              "The input file '",
              base,
              "' appears to already be normalized with '",
              normPrefix,
              "', but `limpaPostQuantNorm='",
              limpaPostQuantNormUse,
              "'` would apply the same normalization again.\n",
              "Use `limpaPostQuantNorm='none'` when supplying a pre-normalized matrix, ",
              "or use the 'log2-normalized' matrix and set `limpaPostQuantNorm` to the desired method."
            )
          }

          warning(
            "The input file '",
            base,
            "' appears to already be normalized with '",
            normPrefix,
            "', but `limpaPostQuantNorm='",
            limpaPostQuantNormUse,
            "'` will apply an additional normalization. ",
            "This may be unintended double-normalization.",
            call. = FALSE
          )
        }
      }
    }
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
    limpaQuantifiedRds = limpaQuantifiedRds,
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
