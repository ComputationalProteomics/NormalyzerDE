# Load raw design into dataframe

Takes a design path, loads the matrix and ensures that the sample column
is in character format and that the group column is in factor format.

## Usage

``` r
loadDesign(designPath, sampleCol = "sample", groupCol = "group")
```

## Arguments

- designPath:

  File path to design matrix.

- sampleCol:

  Column name for column containing sample names.

- groupCol:

  Column name for column containing condition levels.

## Value

designMatrix Design data loaded into data frame

## Examples

``` r
if (FALSE) { # \dontrun{
df <- loadDesign("design.tsv")
} # }
```
