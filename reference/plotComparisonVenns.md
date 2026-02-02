# If multiple comparisons - Show overlap in Venn diagrams

If multiple comparisons - Show overlap in Venn diagrams

## Usage

``` r
plotComparisonVenns(
  nst,
  jobName,
  currentLayout,
  pageno,
  sigThres = 0.1,
  sigThresType = "fdr",
  log2FoldThres = 0,
  maxContrasts = 4
)
```

## Arguments

- nst:

  NormalyzerDE statistics object.

- jobName:

  Name of processing run.

- currentLayout:

  Layout used for document.

- pageno:

  Current page number.

- sigThres:

  Cutoff value for significance theshold

- sigThresType:

  Type of significance cutoff

- log2FoldThres:

  Log2-fold based cutoff threshold

- maxContrasts:

  Maximum contrasts to show pairwise comparisons for

## Value

None
