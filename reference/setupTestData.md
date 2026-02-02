# Generate a random test dataset with features, sample values and retention times

Generate a random test dataset with features, sample values and
retention times

## Usage

``` r
setupTestData(nSamples, nFeatures, rtMin = 40, rtMax = 80, mean = 20, sd = 4)
```

## Arguments

- nSamples:

  Number of samples

- nFeatures:

  Number of features

- rtMin:

  Minimum retention time

- rtMax:

  Maximum retention time

- mean:

  Mean value for sample intensities

- sd:

  Standard deviation for sample intensities

## Value

Test dataset

## Examples

``` r
df <- setupTestData(6, 20)
df <- setupTestData(6, 20, mean=15, sd=1)
```
