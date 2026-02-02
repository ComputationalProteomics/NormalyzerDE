# Verify that samples contain at least a lowest number of values

Verify that samples contain at least a lowest number of values

## Usage

``` r
getLowCountSampleFiltered(
  dataMatrix,
  groups,
  threshold = 15,
  stopIfTooFew = TRUE
)
```

## Arguments

- dataMatrix:

  Dataframe with processed input data.

- groups:

  Vector containing condition levels.

- threshold:

  Lowest number of allowed values in a column.

- stopIfTooFew:

  Abort run if lower than threshold number of values in column

## Value

None
