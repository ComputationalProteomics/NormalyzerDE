# Calculate measures for normalization results

This function prepares an NormalyzerEvaluationResults object containing
the evaluation measures CV (coefficient of variance), MAD (median
absolute deviation), average variance, significance measures (ANOVA
between condition groups) and correlation between replicates.

## Usage

``` r
analyzeNormalizations(nr, categoricalAnova = FALSE)
```

## Arguments

- nr:

  Normalyzer results object with calculated results.

- categoricalAnova:

  Whether categorical or numerical (ordered) ANOVA should be calculated.

## Value

Normalyzer results with attached evaluation results object.

## Examples

``` r
data(example_summarized_experiment)
normObj <- getVerifiedNormalyzerObject("job_name", example_summarized_experiment)
#> Input data checked. All fields are valid.
#> Sample check: More than one sample group found
#> Sample replication check: All samples have replicates
#> RT annotation column found (5)
normResults <- normMethods(normObj)
normResultsWithEval <- analyzeNormalizations(normResults)
```
