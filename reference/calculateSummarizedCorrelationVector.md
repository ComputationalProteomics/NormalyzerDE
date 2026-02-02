# Calculates correlation values between replicates for each condition matrix. Finally returns a matrix containing the results for all dataset

Calculates correlation values between replicates for each condition
matrix. Finally returns a matrix containing the results for all dataset

## Usage

``` r
calculateSummarizedCorrelationVector(
  methodlist,
  allReplicateGroups,
  sampleGroupsWithReplicates,
  corrType
)
```

## Arguments

- methodlist:

  List containing normalized matrices for each normalization method

- allReplicateGroups:

  Vector with condition groups matching the columns found in the
  normalization methods

- sampleGroupsWithReplicates:

  Unique vector with condition groups present in two or more samples

- corrType:

  Type of correlation (Pearson or Spearman)

## Value

avgCorSum Matrix with column corresponding to normalization approaches
and rows corresponding to replicate group
