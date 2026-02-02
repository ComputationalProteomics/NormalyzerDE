# Intensity of each variable in a given sample is divided by the median of intensities of all variables in the sample and then multiplied by the mean of median of sum of intensities of all variables in all samples. The normalized data is then log2-transformed.

Intensity of each variable in a given sample is divided by the median of
intensities of all variables in the sample and then multiplied by the
mean of median of sum of intensities of all variables in all samples.
The normalized data is then log2-transformed.

## Usage

``` r
medianNormalization(rawMatrix, noLogTransform = FALSE)
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
normMatrix <- medianNormalization(example_data_only_values)
```
