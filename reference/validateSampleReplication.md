# Check whether all samples have replicates

Check whether all samples have replicates

## Usage

``` r
validateSampleReplication(
  dataMatrix,
  groups,
  requireReplicates = TRUE,
  quiet = FALSE
)
```

## Arguments

- dataMatrix:

  Prepared matrix containing expression data.

- groups:

  Vector containing condition levels

- requireReplicates:

  By default stops processing if not all samples have replicates

## Value

None
