# Generates a number of visualizations for the performance measures calculated for the normalized matrices. These contain both general measures and direct comparisons for different normalization approaches.

They include:

## Usage

``` r
generatePlots(nr, jobdir, plotRows = 3, plotCols = 4, writeAsPngs = FALSE)
```

## Arguments

- nr:

  Normalyzer results object.

- jobdir:

  Path to output directory for run.

- plotRows:

  Number of plot rows.

- plotCols:

  Number of plot columns.

- writeAsPngs:

  Output the report as PNG-plots instead of a single PDF

## Value

None

## Details

"Total intensity" Barplot showing the summed intensity in each sample
for thelog2-transformed data

"Total missing" Barplot showing the number of missing values found in
each sample for the log2-tranformed data

Log2-MDS plot: MDS plot where data is reduced to two dimensions allowing
inspection of the main global changes in the data

PCV - Intragroup: Mean of intragroup CV of all replicate groups

PMAD - Intragroup: Mean of intragroup median absolute deviation across
replicate groups

PEV - Intragroup: Mean of intragroup pooled estimate of variance across
the replicate groups

Relative PCV, PMAD and PEV compared to log2: The results from PCV, PMAD
and PEV from all normalized data compared to the log2 data

Stable variables plot: 5 analysis of log2 transformed data. Thereafter,
global CV of these variables is estimated from different normalized
datasets. A plot of global CV of the stable variables from all datsets
on the y-axis and PCV-compared to log2 on the x-axis is generated.

CV vs Raw Intensity plots: For the first replicate group in each of the
normalized dataset, a plot of PCV of each variable compared to the
average intensity of the variable in the replicate group is plotted.

MA plots: Plotted using the plotMA function of the limma package. The
first sample in each dataset is plotted against the average of the
replicate group that sample belong to.

Scatterplots: The first two samples from each dataset are plotted.

Q-Q plots: QQ-plots are plotted for the first sample in each normalized
dataset.

Boxplots: Boxplots for all samples are plotted and colored according to
the replicate grouping.

Relative Log Expression (RLE) plots: Relative log expression value
plots. Ratio between the expression of the variable and the median
expression of this variable across all samples. The samples should be
aligned around zero. Any deviation would indicate discrepancies in the
data.

Density plots: Density distributions for each sample using the density
function. Can capture outliers (if single densities lies far from the
others) and see if there is batch effects in the dataset (if for
instance there is two clear collections of lines in the data).

MDS plots Multidimensional scaling plot using the cmdscale() function
from the stats package. Is often able to show whether replicates group
together, and whether there are any clear outliers in the data.

MeanSDplots Displays the standard deviation values against values
ordered according to mean. If no dependency on mean is present (as is
desired) a flat red line is shown.

Pearson and Spearman correlation Mean of intragroup Pearson and Spearman
correlation values for each method.

Dendograms Generated using the hclust function. Data is centered and
scaled prior to analysis. Coloring of replicates is done using as.phylo
from the ape package.

P-value histograms Histogram plots of p-values after calculating an
ANOVA between different condition groups. If no effect is present in the
data a flat distribution is expected. If an effect is present a flat
distribution is still expected, but with a sharp peak close to zero. If
other effects are present it might indicate that the data doesn't
support the assumptions of ANOVA, for instance if there are batch
effects present in the data.

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
outputDir <- tempdir()
generatePlots(normResultsWithEval, outputDir)
#> Warning: `aes_string()` was deprecated in ggplot2 3.0.0.
#> ℹ Please use tidy evaluation idioms with `aes()`.
#> ℹ See also `vignette("ggplot2-in-packages")` for more information.
#> ℹ The deprecated feature was likely used in the vsn package.
#>   Please report the issue to the authors.
#> agg_record_19ab43c0ae68 
#>                       2 
```
