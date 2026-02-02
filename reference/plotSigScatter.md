# Takes an NormalyzerStatistics instance and generates and prints a volcano plot

Takes an NormalyzerStatistics instance and generates and prints a
volcano plot

## Usage

``` r
plotSigScatter(
  nst,
  jobName,
  currentLayout,
  pageno,
  type = "Volcano",
  sigThres = 0.1,
  sigThresType = "fdr",
  log2FoldThres = 0
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

- type:

  Specify whether to plot 'Volcano' or 'MA'.

- sigThres:

  FDR threshold for DE coloring.

## Value

None
