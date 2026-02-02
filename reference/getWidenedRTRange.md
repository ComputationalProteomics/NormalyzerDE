# Pick datapoints before and after window until a minimum number is reached Expects the start and end retention times to match actual retention times present in the data

Pick datapoints before and after window until a minimum number is
reached Expects the start and end retention times to match actual
retention times present in the data

## Usage

``` r
getWidenedRTRange(
  rtStart,
  rtEnd,
  minimumDatapoints,
  retentionTimes,
  allowTooWideData = FALSE
)
```

## Arguments

- rtStart:

  Original retention time start point

- rtEnd:

  Original retention time end point

- minimumDatapoints:

  Required number of datapoints to fulfill

- retentionTimes:

  Vector with all retention times

## Value

Vector with start and end of new RT range
