# Cyclic Loess normalization

Log2 transformed data is normalized by Loess method using the function
"normalizeCyclicLoess". Further information is available for the
function "normalizeCyclicLoess" in the Limma package.

## Usage

``` r
performCyclicLoessNormalization(rawMatrix, noLogTransform = FALSE)
```

## Arguments

- rawMatrix:

  Target matrix to be normalized

- noLogTransform:

  Assumes no need for log transformation

## Value

Normalized matrix

## Examples

``` r
data(example_data_only_values_small)
normMatrix <- performCyclicLoessNormalization(example_data_only_values)
```
