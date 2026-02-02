# Generate an annotated data frame from statistics object

Extracts key values (p-value, adjusted p-value, log2-fold change and
average expression values) from an NormalyzerStatistics instance and
appends these to the annotation- and data-matrices

## Usage

``` r
generateAnnotatedMatrix(nst, prefixSep = "_", compLabels = NULL)
```

## Arguments

- nst:

  NormalyzerDE statistics object.

- prefixSep:

  Character string for separating the prefix names from the statistics
  suffix

- compLabels:

  Vector containing strings to use as prefix for statistical comparisons

## Value

outDf Annotated statistics matrix

## Examples

``` r
data(example_stat_summarized_experiment)
statObj <- NormalyzerStatistics(example_stat_summarized_experiment)
statObj <- calculateContrasts(statObj, comparisons=c("1-2", "2-3"), condCol="group", type="limma")
annotDf <- generateAnnotatedMatrix(statObj)
```
