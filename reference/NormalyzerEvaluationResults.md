# Representation of evaluation results by calculating performance measures for an an NormalyzerResults instance

Contains the resulting information from the processing which
subsequently can be used to generate the quality assessment report.

## Usage

``` r
NormalyzerEvaluationResults(nr, categoricalAnova = TRUE)

NormalyzerEvaluationResults(nr, categoricalAnova = TRUE)
```

## Arguments

- nr:

  NormalyzerResults object

- categoricalAnova:

  Should ANOVA be categorical or not

## Value

nds Generated NormalyzerEvaluationResults instance

## Slots

- `avgcvmem`:

  Average coefficient of variance per method

- `avgcvmempdiff`:

  Percentage difference of mean coefficient of variance compared to
  log2-transformed data

- `featureCVPerMethod`:

  CV calculated per feature and normalization method.

- `avgmadmem`:

  Average median absolute deviation

- `avgmadmempdiff`:

  Percentage difference of median absolute deviation compared to
  log2-transformed data

- `avgvarmem`:

  Average variance per method

- `avgvarmempdiff`:

  Percentage difference of mean variance compared to log2-transformed
  data

- `lowVarFeaturesCVs`:

  List of 5 for log2-transformed data

- `lowVarFeaturesCVsPercDiff`:

  Coefficient of variance for least variable entries

- `anovaP`:

  ANOVA calculated p-values

- `repCorPear`:

  Within group Pearson correlations

- `repCorSpear`:

  Within group Spearman correlations

## Examples

``` r
data(example_summarized_experiment)
normObj <- getVerifiedNormalyzerObject("job_name", example_summarized_experiment)
#> Input data checked. All fields are valid.
#> Sample check: More than one sample group found
#> Sample replication check: All samples have replicates
#> RT annotation column found (5)
normResults <- normMethods(normObj)
normEval <- NormalyzerEvaluationResults(normResults)
```
