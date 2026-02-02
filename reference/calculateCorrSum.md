# Calculates internal correlations for each condition having at least two samples and returns a vector with correlation values corresponding to each condition

Calculates internal correlations for each condition having at least two
samples and returns a vector with correlation values corresponding to
each condition

## Usage

``` r
calculateCorrSum(
  methodData,
  allReplicateGroups,
  sampleGroupsWithReplicates,
  corrType
)
```

## Arguments

- methodData:

  Expression data matrix

- allReplicateGroups:

  Full condition header corresponding to data tables columns

- sampleGroupsWithReplicates:

  Unique conditions where number of replicates exceeds one

- corrType:

  Type of correlation (Pearson or Spearman)

## Value

corSums
