# Calculate average MAD (Median Absolute Deviation) for each feature in each condition and then calculates the average for each replicate group

Calculate average MAD (Median Absolute Deviation) for each feature in
each condition and then calculates the average for each replicate group

## Usage

``` r
calculateAvgMadMem(methodList, sampleReplicateGroups)
```

## Arguments

- methodList:

  List containing normalized matrices.

- sampleReplicateGroups:

  Condition header.

## Value

condAvgMadMat Matrix with average MAD for each biological condition.
