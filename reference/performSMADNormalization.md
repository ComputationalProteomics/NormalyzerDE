# Median absolute deviation normalization Normalization subtracts the median and divides the data by the median absolute deviation (MAD).

Median absolute deviation normalization Normalization subtracts the
median and divides the data by the median absolute deviation (MAD).

## Usage

``` r
performSMADNormalization(rawMatrix, noLogTransform = FALSE)
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
normMatrix <- performSMADNormalization(example_data_only_values)
```
