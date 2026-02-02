# Class representing a dataset for statistical processing in NormalyzerDE

Is initialized with an annotation matrix, a data matrix and a design
data frame. This object can subsequently be processed to generate
statistical values and in turn used to write a full matrix with
additional statistical information as well as a graphical report of the
comparisons.

## Usage

``` r
NormalyzerStatistics(experimentObj, logTrans = FALSE)

NormalyzerStatistics(experimentObj, logTrans = FALSE)
```

## Arguments

- experimentObj:

  Instance of SummarizedExperiment containing matrix and design
  information as column data

- logTrans:

  Whether the input data should be log transformed

## Value

nds Generated NormalyzerStatistics instance

## Slots

- `annotMat`:

  Matrix containing annotation information

- `dataMat`:

  Matrix containing (normalized) expression data

- `filteredDataMat`:

  Filtered matrix with low-count rows removed

- `designDf`:

  Data frame containing design conditions

- `filteringContrast`:

  Vector showing which entries are filtered (due to low count)

- `pairwiseCompsP`:

  List with P-values for pairwise comparisons

- `pairwiseCompsFdr`:

  List with FDR-values for pairwise comparisons

- `pairwiseCompsAve`:

  List with average expression values

- `pairwiseCompsFold`:

  List with log2 fold-change values for pairwise comparisons

- `contrasts`:

  Spot for saving vector of last used contrasts

- `condCol`:

  Column containing last used conditions

- `batchCol`:

  Column containing last used batch conditions

## Examples

``` r
data(example_stat_summarized_experiment)
nst <- NormalyzerStatistics(example_stat_summarized_experiment)
```
