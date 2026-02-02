# Prepare SummarizedExperiment object for raw data to be normalized containing data, design and annotation information

Prepare SummarizedExperiment object for raw data to be normalized
containing data, design and annotation information

## Usage

``` r
setupRawDataObject(
  dataPath,
  designPath,
  inputFormat = "default",
  zeroToNA = FALSE,
  sampleColName = "sample",
  groupColName = "group"
)
```

## Arguments

- dataPath:

  File path to data matrix.

- designPath:

  File path to design matrix.

- inputFormat:

  Type of matrix for data, can be either 'default', 'proteios',
  'maxquantprot' or 'maxquantpep'

- zeroToNA:

  If TRUE zeroes in the data is automatically converted to NA values

- sampleColName:

  Column name for column containing sample names

- groupColName:

  Column name for column containing condition levels

## Value

experimentObj SummarizedExperiment object loaded with the data

## Examples

``` r
data_path <- system.file(package="NormalyzerDE", "extdata", "tiny_data.tsv")
design_path <- system.file(package="NormalyzerDE", "extdata", "tiny_design.tsv")
df <- setupRawDataObject(data_path, design_path)
```
