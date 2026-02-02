# Performs statistical comparisons between the supplied conditions. It uses the design matrix and data matrix in the supplied NormalyzerStatistics object. A column is supplied specifying which of the columns in the design matrix that is used for deciding the sample groups. The comparisons vector specifies which pairwise comparisons between condition levels that are to be calculated.

Optionally, a batch column can be specified allowing compensation for
covariate variation in the statistical model. This is only compatible
with a Limma-based statistical analysis.

## Usage

``` r
calculateContrasts(
  nst,
  comparisons,
  condCol,
  batchCol = NULL,
  splitter = "-",
  type = "limma",
  leastRepCount = 1
)

# S4 method for class 'NormalyzerStatistics'
calculateContrasts(
  nst,
  comparisons,
  condCol,
  batchCol = NULL,
  splitter = "-",
  type = "limma",
  leastRepCount = 1
)
```

## Arguments

- nst:

  Results evaluation object.

- comparisons:

  String with comparisons for contrasts.

- condCol:

  Column name in design matrix containing condition information.

- batchCol:

  Column name in design matrix containing batch information.

- splitter:

  Character dividing contrast conditions.

- type:

  Type of statistical test (Limma or welch).

- leastRepCount:

  Least replicates in each group to be retained for contrast
  calculations

## Value

nst Statistics object with statistical measures calculated

## Examples

``` r
data(example_stat_summarized_experiment)
nst <- NormalyzerStatistics(example_stat_summarized_experiment)
results <- calculateContrasts(nst, c("1-2", "2-3"), "group")
resultsBatch <- calculateContrasts(nst, c("1-2", "2-3"), "group", batchCol="batch")
```
