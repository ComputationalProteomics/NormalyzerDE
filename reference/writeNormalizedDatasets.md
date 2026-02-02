# Write normalization matrices to file

Outputs each of the normalized datasets to the specified directory.

## Usage

``` r
writeNormalizedDatasets(
  nr,
  jobdir,
  includePairwiseComparisons = FALSE,
  includeCvCol = FALSE,
  includeAnovaP = FALSE,
  normSuffix = "-normalized.txt",
  rawdataName = "submitted_rawdata.txt"
)
```

## Arguments

- nr:

  Results object.

- jobdir:

  Path to output directory.

- includePairwiseComparisons:

  Include p-values for pairwise comparisons.

- includeCvCol:

  Include CV column in output.

- includeAnovaP:

  Include ANOVA p-value in output.

- normSuffix:

  String used to name output together with normalization names.

- rawdataName:

  Name of output raw data file.

## Value

None

## Examples

``` r
data(example_summarized_experiment)
normObj <- getVerifiedNormalyzerObject("job_name", example_summarized_experiment)
#> Input data checked. All fields are valid.
#> Sample check: More than one sample group found
#> Sample replication check: All samples have replicates
#> RT annotation column found (5)
normResults <- normMethods(normObj)
normResultsWithEval <- analyzeNormalizations(normResults)
outputDir <- tempdir()
writeNormalizedDatasets(normResultsWithEval, outputDir)
```
