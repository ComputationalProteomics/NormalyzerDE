# Load raw data into dataframe

General function which allows specifying different types of input data
including "proteios", "maxquantpep" (peptide output from MaxQuant) and
"maxquantprot" (protein output from MaxQuant) formats.

## Usage

``` r
loadData(dataPath, inputFormat = "default")
```

## Arguments

- dataPath:

  File path to design matrix.

- inputFormat:

  If input is given in standard NormalyzerDE format, Proteios format or
  in MaxQuant protein or peptide format

## Value

rawData Raw data loaded into data frame

## Examples

``` r
if (FALSE) { # \dontrun{
df <- loadData("data.tsv")
} # }
```
