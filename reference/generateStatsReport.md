# Generate full output report plot document. Plots p-value histograms for each contrast in the NormalyzerStatistics instance and writes these to a PDF report.

Generate full output report plot document. Plots p-value histograms for
each contrast in the NormalyzerStatistics instance and writes these to a
PDF report.

## Usage

``` r
generateStatsReport(
  nst,
  jobName,
  jobDir,
  sigThres = 0.1,
  sigThresType = "fdr",
  log2FoldThres = 0,
  plotRows = 3,
  plotCols = 4,
  writeAsPngs = FALSE
)
```

## Arguments

- nst:

  NormalyzerDE statistics object.

- jobName:

  Name of processing run.

- jobDir:

  Path to output directory.

- sigThres:

  Significance threshold for indicating as significant

- sigThresType:

  Type of significance threshold (FDR or p)

- log2FoldThres:

  log2 fold-change required for being counted as significant

- plotRows:

  Number of plot rows.

- plotCols:

  Number of plot columns.

- writeAsPngs:

  Output the report as separate PNG files instead of a single PDF file

## Value

None

## Examples

``` r
data(example_stat_summarized_experiment)
statObj <- NormalyzerStatistics(example_stat_summarized_experiment)
statObj <- calculateContrasts(statObj, comparisons=c("1-2", "2-3"), 
  condCol="group", type="limma")
outputDir <- tempdir()
generateStatsReport(statObj, "jobName", outputDir)
#> Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
#> ℹ Please use `linewidth` instead.
#> ℹ The deprecated feature was likely used in the NormalyzerDE package.
#>   Please report the issue to the authors.
#> agg_record_288826ed677d 
#>                       2 
```
