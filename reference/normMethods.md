# Perform normalizations on Normalyzer dataset

Perform normalizations on Normalyzer dataset

## Usage

``` r
normMethods(
  nds,
  forceAll = FALSE,
  normalizeRetentionTime = TRUE,
  quiet = FALSE,
  rtStepSizeMinutes = 1,
  rtWindowMinCount = 100,
  rtWindowShifts = 1,
  rtWindowMergeMethod = "mean",
  noLogTransform = FALSE
)
```

## Arguments

- nds:

  Normalyzer dataset object.

- forceAll:

  Force all methods to run despite not qualifying for thresholds.

- normalizeRetentionTime:

  Perform retention time based normalization methods.

- quiet:

  Prevent diagnostic output

- rtStepSizeMinutes:

  Retention time normalization window size.

- rtWindowMinCount:

  Minimum number of datapoints in each retention-time segment.

- rtWindowShifts:

  Number of layered retention time normalized windows.

- rtWindowMergeMethod:

  Merge approach for layered retention time windows.

- noLogTransform:

  Per default NormalyzerDE performs a log-transformation on the input
  data. If not needed, specify this option

## Value

Returns Normalyzer results object with performed analyzes assigned as
attributes

## Examples

``` r
data(example_summarized_experiment)
normObj <- getVerifiedNormalyzerObject("job_name", example_summarized_experiment)
#> Input data checked. All fields are valid.
#> Sample check: More than one sample group found
#> Sample replication check: All samples have replicates
#> RT annotation column found (5)
normResults <- normMethods(normObj)
```
