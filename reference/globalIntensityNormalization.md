# The normalization divides the intensity of each variable in a sample with the sum of intensities of all variables in the sample and multiplies with the median of sum of intensities of all variables in all samples. The normalized data is then log2-transformed.

The normalization divides the intensity of each variable in a sample
with the sum of intensities of all variables in the sample and
multiplies with the median of sum of intensities of all variables in all
samples. The normalized data is then log2-transformed.

## Usage

``` r
globalIntensityNormalization(rawMatrix, noLogTransform = FALSE)
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
normMatrix <- globalIntensityNormalization(example_data_only_values)
```
