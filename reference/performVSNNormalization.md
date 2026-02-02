# Log2 transformed data is normalized using the function "justvsn" from the VSN package.

The VSN (Variance Stabilizing Normalization) attempts to transform the
data in such a way that the variance remains nearly constant over the
intensity spectrum

## Usage

``` r
performVSNNormalization(rawMatrix)
```

## Arguments

- rawMatrix:

  Target matrix to be normalized

## Value

Normalized matrix

## Examples

``` r
data(example_data_only_values_small)
normMatrix <- performVSNNormalization(example_data_only_values)
```
