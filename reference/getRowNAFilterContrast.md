# Get contrast vector (TRUE/FALSE-values) indicating whether both at least half values are present, and each sample has at least one non-NA value

Get contrast vector (TRUE/FALSE-values) indicating whether both at least
half values are present, and each sample has at least one non-NA value

## Usage

``` r
getRowNAFilterContrast(dataMatrix, replicateHeader, minCount = 1)
```

## Arguments

- dataMatrix:

  Matrix with expression values for entities in replicate samples.

- replicateHeader:

  Header showing how samples in matrix are replicated.

- minCount:

  Minimum number of required values present in samples.

## Value

Contrast vector
