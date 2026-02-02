# NormalyzerDE pipeline entry point

This function is the main execution point for the normalization part of
the NormalyzerDE analysis pipeline. When executed it performs the
following steps:

## Usage

``` r
normalyzer(
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
  skipAnalysis = FALSE,
  quiet = FALSE,
  noLogTransform = FALSE,
  writeReportAsPngs = FALSE,
  rtStepSizeMinutes = 1,
  rtWindowMinCount = 100,
  rtWindowShifts = 1,
  rtWindowMergeMethod = "mean"
)
```

## Arguments

- jobName:

  Give the current run a name.

- designPath:

  Path to file containing design matrix.

- dataPath:

  Specify an output directory for generated files. Defaults to current
  working directory.

- experimentObj:

  SummarizedExperiment object, can be provided as input as alternative
  to 'designPath' and 'dataPath'

- outputDir:

  Directory where results folder is created.

- forceAllMethods:

  Debugging function. Run all normalizations even if they aren't in the
  recommended range of number of values

- omitLowAbundSamples:

  Automatically remove samples with fewer non-NA values compared to
  threshold given by sampleAbundThres. Will otherwise stop with error
  message if such sample is encountered.

- sampleAbundThres:

  Threshold for omitting low-abundant samples. Is by default set to 15.

- tinyRunThres:

  If total number of features is less than this, a limited run is
  performed.

- requireReplicates:

  Require multiple samples per condition to pass input validation.

- normalizeRetentionTime:

  Perform normalizations over retention time.

- plotRows:

  Number of plot-rows in output documentation.

- plotCols:

  Number of plot-columns in output documentation.

- zeroToNA:

  Convert zero values to NA.

- sampleColName:

  Column name in design matrix containing sample IDs.

- groupColName:

  Column name in design matrix containing condition IDs.

- inputFormat:

  Type of input format.

- skipAnalysis:

  Only perform normalization steps.

- quiet:

  Omit status messages printed during run.

- noLogTransform:

  Don't log-transform the input.

- writeReportAsPngs:

  Output the evaluation report as PNG files instead of a single PDF

- rtStepSizeMinutes:

  Retention time normalization window size.

- rtWindowMinCount:

  Minimum number of datapoints in each retention-time segment.

- rtWindowShifts:

  Number of layered retention time normalized windows.

- rtWindowMergeMethod:

  Merge approach for layered retention time windows.

## Value

None

## Details

1: Loads the data matrix containing expression values and optional
annotations, as well as the design matrix containing the experimental
setup 2: Performs input data verification to validate that the data is
in correct format. This step captures many common formatting errors. It
returns an instance of the NormalyzerDataset class representing the
unprocessed data. 3: Calculate a range of normalizations for the
dataset. The result is provided as a NormalyzerResults object containing
the resulting data matrices from each normalization. 4: Analyze the
normalizations and generate performance measures for each of the
normalized datasets. This result is provided as a
NormalyzerEvaluationResults object. 5: Output the matrices containing
the normalized datasets to files. 6: Generate visualizations overviewing
the performance measures and write them to a PDF report.

## Examples

``` r
if (FALSE) { # \dontrun{
data_path <- system.file(package="NormalyzerDE", "extdata", "tiny_data.tsv")
design_path <- system.file(package="NormalyzerDE", "extdata", "tiny_design.tsv")
out_dir <- tempdir()
normalyzer(
    jobName="my_jobname", 
    designPath=design_path, 
    dataPath=data_path, 
    outputDir=out_dir)
normalyzer(
    "my_jobname", 
    designMatrix="design.tsv", 
    "data.tsv", 
    outputDir="path/to/output", 
    normalizeRetentionTime=TRUE, 
    retentionTimeWindow=2)
normalyzer(
    "my_jobname", 
    designMatrix="design.tsv", 
    "data.tsv", 
    outputDir="path/to/output", 
    inputFormat="maxquantprot")
} # }
```
