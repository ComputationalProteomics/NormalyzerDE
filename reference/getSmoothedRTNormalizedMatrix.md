# Generate multiple RT time-window normalized matrices where one is shifted. Merge them using a specified method (mean or median) and return the result.

Uses the function getRTNormalizedMatrix to generate multiple normalized
matrices which are shifted respective to each other and finally merged
into a single matrix. This could potentially reduce effect of
fluctuations within individual windows.

## Usage

``` r
getSmoothedRTNormalizedMatrix(
  rawMatrix,
  retentionTimes,
  normMethod,
  stepSizeMinutes,
  windowShifts = 2,
  windowMinCount = 100,
  mergeMethod = "mean",
  noLogTransform = FALSE
)
```

## Arguments

- rawMatrix:

  Target matrix to be normalized

- retentionTimes:

  Vector of retention times corresponding to rawMatrix

- normMethod:

  The normalization method to apply to the time windows

- stepSizeMinutes:

  Size of windows to be normalized

- windowShifts:

  Number of frame shifts.

- windowMinCount:

  Minimum number of features within window.

- mergeMethod:

  Layer merging approach. Mean or median.

- noLogTransform:

  Don't log transform the input

## Value

Normalized matrix

## Examples

``` r
data(example_data_small)
data(example_data_only_values)
data(example_design_small)
retentionTimes <- as.numeric(example_data[, "Average.RT"])
dataMat <- example_data_only_values
performCyclicLoessNormalization <- function(rawMatrix) {
    log2Matrix <- log2(rawMatrix)
    normMatrix <- limma::normalizeCyclicLoess(log2Matrix, method="fast")
    colnames(normMatrix) <- colnames(rawMatrix)
    normMatrix
}
rtNormMat <- getSmoothedRTNormalizedMatrix(dataMat, retentionTimes, 
    performCyclicLoessNormalization, stepSizeMinutes=1, windowMinCount=100, 
    windowShifts=2, mergeMethod="median")
```
