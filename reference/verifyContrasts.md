# Check that a given contrast string is valid given a particular design matrix. Each level tested for in the contrast should be present in the condition column for the design matrix.

Mainly meant to verify strings received during server usage.

## Usage

``` r
verifyContrasts(designLevels, contrasts)
```

## Arguments

- designLevels:

  Vector containing condition levels present in design

- contrasts:

  A string containing one or several (comma delimited) strings for which
  contrasts should be performed

## Value

None
