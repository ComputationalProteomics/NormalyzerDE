# Uses a list of FDR-values to extract features with low variance in the log2-transformed dataset. This is then used to calculate the average CV for these 'lowly variable' features in each normalization approach

Uses a list of FDR-values to extract features with low variance in the
log2-transformed dataset. This is then used to calculate the average CV
for these 'lowly variable' features in each normalization approach

## Usage

``` r
findLowlyVariableFeaturesCVs(referenceFDR, methodList)
```

## Arguments

- referenceFDR:

  List of FDR values used as non-normalized reference

- methodList:

  List containing normalized matrices

## Value

lowVarFeaturesAverageCVs Average CV values for lowly variable features
in each normalization approach
