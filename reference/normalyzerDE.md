# NormalyzerDE differential expression

Performs differential expression analysis on a normalization matrix.
This command executes a pipeline processing the data and generates an
annotated normalization matrix and a report containing p-value
histograms for each of the performed comparisons.

## Usage

``` r
normalyzerDE(
  jobName,
  comparisons,
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
  quiet = FALSE,
  sigThres = 0.1,
  sigThresType = "fdr",
  log2FoldThres = 0,
  writeReportAsPngs = FALSE
)
```

## Arguments

- jobName:

  Name of job

- comparisons:

  Character vector containing target contrasts. If comparing condA with
  condB, then the vector would be c("condA-condB")

- designPath:

  File path to design matrix

- dataPath:

  File path to normalized matrix

- experimentObj:

  SummarizedExperiment object, can be provided as input as alternative
  to 'designPath' and 'dataPath'

- outputDir:

  Path to output directory

- logTrans:

  Log transform the input (needed if providing non-logged input)

- type:

  Type of statistical comparison, "limma", "limma_intensity" or "welch",
  where "limma_intensity" allows the prior to be fit according to
  intensity rather than using a flat prior

- sampleCol:

  Design matrix column header for column containing sample IDs

- condCol:

  Design matrix column header for column containing sample conditions

- batchCol:

  Provide an optional column for inclusion of possible batch variance in
  the model

- techRepCol:

  Design matrix column header for column containing technical replicates

- leastRepCount:

  Minimum required replicate count

- quiet:

  Omit status messages printed during run

- sigThres:

  Significance threshold use for illustrating significant hits in
  diagnostic plots

- sigThresType:

  Type of significance threshold, "fdr" or "p". "fdr" is strongly
  recommended (Benjamini-Hochberg corrected p-values)

- log2FoldThres:

  Fold-size cutoff for being considered significant in diagnostic plots

- writeReportAsPngs:

  Output report as separate PNG files instead of a single PDF

## Value

None

## Details

When executed, it performs the following steps:

1: Read the data and the design matrices into dataframes. 2: Generate an
instance of the NormalyzerStatistics class representing the data and
their statistical comparisons. 3: Optionally reduce technical replicates
in both the data matrix and the design matrix 4: Calculate statistical
contrats between supplied groups 5: Generate an annotated version of the
original dataframe where columns containing statistical key measures
have been added 6: Write the table to file 7: Generate a PDF report
displaying p-value histograms for each calculated contrast

## Examples

``` r
data_path <- system.file(package="NormalyzerDE", "extdata", "tiny_data.tsv")
design_path <- system.file(package="NormalyzerDE", "extdata", "tiny_design.tsv")
out_dir <- tempdir()
normalyzerDE(
  jobName="my_jobname", 
  comparisons=c("4-5"), 
  designPath=design_path, 
  dataPath=data_path,
  outputDir=out_dir,
  condCol="group")
#> You are running version 1.23.5 of NormalyzerDE
#> [1] "Setting up statistics object"
#> [1] "Calculating statistical contrasts..."
#> [1] "Contrast calculations done!"
#> [1] "Writing 100 annotated rows to /tmp/RtmpVapksj/my_jobname/my_jobname_stats.tsv"
#> [1] "Writing statistics report"
#> [1] "All done! Results are stored in: /tmp/RtmpVapksj/my_jobname, processing time was 0 minutes"
```
