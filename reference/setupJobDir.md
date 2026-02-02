# Create empty directory for run

Creates a directory at provided path named to the jobname.

## Usage

``` r
setupJobDir(jobName, outputDir)
```

## Arguments

- jobName:

  Name of the run.

- outputDir:

  Path to directory where to create the output directory.

## Value

Path to newly created directory.

## Examples

``` r
setupJobDir("job_name", "path/to/outdir")
#> [1] "path/to/outdir/job_name"
```
