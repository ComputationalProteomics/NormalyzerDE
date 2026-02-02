# Represents raw input data together with basic annotation information

Takes a job name, a data matrix, a design matrix as well as
specification of the group and sample columns in the design matrix.
Provides the basic representation of a dataset in the NormalyzerDE
normalization part.

## Usage

``` r
NormalyzerDataset(
  jobName,
  designMatrix,
  rawData,
  annotationData,
  sampleNameCol,
  groupNameCol,
  tinyRunThres = 50,
  quiet = FALSE
)

NormalyzerDataset(
  jobName,
  designMatrix,
  rawData,
  annotationData,
  sampleNameCol,
  groupNameCol,
  tinyRunThres = 50,
  quiet = FALSE
)
```

## Arguments

- jobName:

  Name of the NormalyzerDE processing run

- designMatrix:

  Matrix containing sample conditions

- rawData:

  Matrix containing raw input data

- annotationData:

  Matrix containing annotation information for each input feature. Is
  expected to contain the same number of rows as the data but can
  contain any number of features

- sampleNameCol:

  Name of column in design matrix containing sample information

- groupNameCol:

  Name of column in design matrix containing condition information

- tinyRunThres:

  If fewer features than this is present in the input a limited run will
  be performed to avoid some steps requiring a more extensive number of
  features.

- quiet:

  If set to TRUE no information messages will be printed

## Value

nds Generated NormalyzerDataset instance

## Slots

- `jobName`:

  Name of the job represented by the dataset.

- `rawData`:

  Matrix with raw values.

- `sampleNameCol`:

  Name column for sample.

- `groupNameCol`:

  Name column for groups.

- `designMatrix`:

  Data frame containing design.

- `sampleNames`:

  Vector containing sample names.

- `filterrawdata`:

  Reduced raw data matrix where low abundance rows are removed

- `sampleReplicateGroups`:

  Vector with sample replicate information

- `samplesGroupsWithReplicates`:

  Vector with replicated sample replicate information

- `annotationValues`:

  Annotation part of original dataframe.

- `retentionTimes`:

  Vector of retention time values.

- `singleReplicateRun`:

  Conditional whether run is single replicate.
