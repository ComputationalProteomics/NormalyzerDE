# Calculate CV per replicate group and normalization technique

Iterates through each normalization method and calculate average CV
values per replicate group.

## Usage

``` r
calculateReplicateCV(methodList, sampleReplicateGroups)
```

## Arguments

- methodList:

  List containing normalized matrices.

- sampleReplicateGroups:

  Condition header.

## Value

avgCVPerNormAndReplicates Matrix with group CVs as rows and
normalization technique as columns
