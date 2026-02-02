# Filter rows with lower than given number of replicates for any condition

Filter rows with lower than given number of replicates for any condition

## Usage

``` r
filterLowRep(df, groups, leastRep = 2)
```

## Arguments

- df:

  Dataframe with expression data to filter

- groups:

  Condition groups header

- leastRep:

  Minimum number of replicates in each group to retain

## Value

collDesignDf Reduced design matrix
