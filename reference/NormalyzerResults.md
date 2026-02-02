# Representation of the results from performing normalization over a dataset

It is linked to a NormalyzerDataset instance representing the raw data
which has been processed. After performing evaluation it also links to
an instance of NormalyzerEvaluationResults representing the results from
the evaluation.

## Usage

``` r
NormalyzerResults(nds)

NormalyzerResults(nds)
```

## Arguments

- nds:

  NormalyzerDataset object

## Value

nr Prepared NormalyzerResults object

## Slots

- `normalizations`:

  SummarizedExperiment object containing calculated normalization
  results

- `nds`:

  Normalyzer dataset representing run data

- `ner`:

  Normalyzer evaluation results for running extended normalizations

## Examples

``` r
data(example_summarized_experiment)
normObj <- getVerifiedNormalyzerObject("job_name", example_summarized_experiment)
#> Input data checked. All fields are valid.
#> Sample check: More than one sample group found
#> Sample replication check: All samples have replicates
#> RT annotation column found (5)
emptyNormResults <- NormalyzerResults(normObj)
```
