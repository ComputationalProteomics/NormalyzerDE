# Global linear regression normalization

Log2 transformed data is normalized by robust linear regression using
the function "rlm" from the MASS package.

## Usage

``` r
performGlobalRLRNormalization(rawMatrix, noLogTransform = FALSE)
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
normMatrix <- performGlobalRLRNormalization(example_data_only_values)
```
