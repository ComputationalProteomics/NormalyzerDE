# Verify that input data is in correct format, and if so, return a generated NormalyzerDE data object from that input data

This function performs a number of checks on the input data and provides
informative error messages if the data isn't fulfilling the required
format. Checks include verifying that the design matrix matches to the
data matrix, that the data matrix contains valid numbers and that
samples have enough values for analysis

## Usage

``` r
getVerifiedNormalyzerObject(
  jobName,
  summarizedExp,
  threshold = 15,
  omitSamples = FALSE,
  requireReplicates = TRUE,
  quiet = FALSE,
  noLogTransform = FALSE,
  tinyRunThres = 50
)
```

## Arguments

- jobName:

  Name of ongoing run.

- summarizedExp:

  Summarized experiment input object

- threshold:

  Minimum number of features.

- omitSamples:

  Automatically omit invalid samples from analysis.

- requireReplicates:

  Require there to be at least to samples per condition

- quiet:

  Don't print output messages during processing

- noLogTransform:

  Don't log-transform the provided data

- tinyRunThres:

  If less features in run, a limited run is performed

## Value

Normalyzer data object representing verified input data.

## Examples

``` r
data(example_summarized_experiment)
normObj <- getVerifiedNormalyzerObject("job_name", example_summarized_experiment)
#> Input data checked. All fields are valid.
#> Sample check: More than one sample group found
#> Sample replication check: All samples have replicates
#> RT annotation column found (5)
```
