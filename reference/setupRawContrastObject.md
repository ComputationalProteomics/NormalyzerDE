# Prepare SummarizedExperiment object for statistics data

Prepare SummarizedExperiment object for statistics data

## Usage

``` r
setupRawContrastObject(dataPath, designPath, sampleColName)
```

## Arguments

- dataPath:

  Path to raw data matrix

- designPath:

  Path to design matrix

- sampleColName:

  Name for column in design matrix containing sample names

## Value

experimentObj Prepared instance of SummarizedExperiment

## Examples

``` r
data_path <- system.file(package="NormalyzerDE", "extdata", "tiny_data.tsv")
design_path <- system.file(package="NormalyzerDE", "extdata", "tiny_design.tsv")
sumExpObj <- setupRawContrastObject(data_path, design_path, "sample")
```
