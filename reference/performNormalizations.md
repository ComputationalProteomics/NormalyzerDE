# Main function for executing normalizations

Main function for executing normalizations

## Usage

``` r
performNormalizations(
  nr,
  forceAll = FALSE,
  rtNorm = FALSE,
  rtStepSizeMinutes = 1,
  rtWindowMinCount = 100,
  rtWindowShifts = 1,
  rtWindowMergeMethod = "median",
  noLogTransform = FALSE,
  quiet = FALSE
)

# S4 method for class 'NormalyzerResults'
performNormalizations(
  nr,
  forceAll = FALSE,
  rtNorm = FALSE,
  rtStepSizeMinutes = 1,
  rtWindowMinCount = 100,
  rtWindowShifts = 1,
  rtWindowMergeMethod = "median",
  noLogTransform = FALSE,
  quiet = FALSE
)
```

## Arguments

- nr:

  Normalyzer results object.

- forceAll:

  Ignore dataset size limits and run all normalizations (only meant for
  testing purposes)

- rtNorm:

  Perform retention time based normalizations

- rtStepSizeMinutes:

  Retention time normalization window size.

- rtWindowMinCount:

  Minimum number of datapoints in each retention-time segment.

- rtWindowShifts:

  Number of layered retention time normalized windows.

- rtWindowMergeMethod:

  Merge approach for layered retention time windows.

- noLogTransform:

  Prevent log-transforming input

- quiet:

  Don't show regular output messages

## Value

nr NormalyzerDE results object
