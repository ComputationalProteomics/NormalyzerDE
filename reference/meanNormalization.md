# Intensity of each variable in a given sample is divided by the mean of sum of intensities of all variables in the sample and then multiplied by the mean of sum of intensities of all variables in all samples. The normalized data is then transformed to log2.

Intensity of each variable in a given sample is divided by the mean of
sum of intensities of all variables in the sample and then multiplied by
the mean of sum of intensities of all variables in all samples. The
normalized data is then transformed to log2.

## Usage

``` r
meanNormalization(rawMatrix, noLogTransform = FALSE)
```

## Arguments

- rawMatrix:

  Target matrix to be normalized

- noLogTransform:

  Assumes no need for log transformation

## Value

Normalized and log-transformed matrix

## Examples

``` r
data(example_data_only_values_small)
normMatrix <- meanNormalization(example_data_only_values)
```
