# Calculates ANOVA p-values comparing the different condition groups returning a vector of resulting p-values with NA-values where too few values were present in at least one of the groups.

Calculates ANOVA p-values comparing the different condition groups
returning a vector of resulting p-values with NA-values where too few
values were present in at least one of the groups.

## Usage

``` r
calculateANOVAPValues(methodList, sampleReplicateGroups, categoricalANOVA)
```

## Arguments

- methodList:

  List containing normalized matrices

- sampleReplicateGroups:

  Condition header

- categoricalANOVA:

  Whether the ANOVA should be calculated using ordered or categorical
  groups

## Value

avgVarianceMat Matrix with average variance for each biological
condition
