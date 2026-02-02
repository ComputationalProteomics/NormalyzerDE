# Quantile normalization is performed by the function "normalize.quantiles" from the package preprocessCore.

It makes the assumption that the data in different samples should
originate from an identical distribution. It does this by generating a
reference distribution and then scaling the other samples accordingly.

## Usage

``` r
performQuantileNormalization(rawMatrix, noLogTransform = FALSE)
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
normMatrix <- performQuantileNormalization(example_data_only_values)
```
