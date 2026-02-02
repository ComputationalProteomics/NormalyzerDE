# Remove technical replicates from data and design

Collapses sample values into their average. If only one value is present
due to NA-values in other technical replicates, then that value is used.

## Usage

``` r
reduceTechnicalReplicates(se, techRepColName, sampleColName)
```

## Arguments

- se:

  Summarized experiment where the assay contains the data to be reduced,
  and the colData the data frame

- techRepColName:

  Technical replicates column name in colData

- sampleColName:

  Sample names column name in colData

## Value

reducedSe Summarized experiment with reduced data

## Details

Takes a SummarizedExperiment where the data is present as the assay and
the colData contains the design conditions. In the design conditions
there should be one column with the technical replicate groups and one
column containing the sample names

## Examples

``` r
testData <- as.matrix(data.frame(
    c(1,1,1), 
    c(1,2,1), 
    c(7,7,7), 
    c(7,9,7)))
colnames(testData) <- c("a1", "a2", "b1", "b2")
designDf <- data.frame(
    sample=c("a1", "a2", "b1", "b2"), 
    techrep=c("a", "a", "b", "b"))
se <- SummarizedExperiment::SummarizedExperiment(
    assay=testData,
    colData=designDf
)
statObj <- reduceTechnicalReplicates(se, "techrep", "sample")
```
