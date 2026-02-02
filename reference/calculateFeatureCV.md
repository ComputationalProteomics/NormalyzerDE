# Calculate CV values for each feature. Iterates through each normalization method and calculates a matrix of CV values where each column correspond to a method and each row corresponds to a feature.

Calculate CV values for each feature. Iterates through each
normalization method and calculates a matrix of CV values where each
column correspond to a method and each row corresponds to a feature.

## Usage

``` r
calculateFeatureCV(methodList)
```

## Arguments

- methodList:

  List containing normalized matrices.

- sampleReplicateGroups:

  Condition header.

## Value

methodFeatureCVMatrix Matrix with feature as rows and normalization
method as columns
