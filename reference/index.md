# Package index

## Main workflows

- [`normalyzerDE()`](https://computationalproteomics.github.io/NormalyzerDE/reference/normalyzerDE.md)
  : NormalyzerDE differential expression
- [`normalyzer()`](https://computationalproteomics.github.io/NormalyzerDE/reference/normalyzer.md)
  : NormalyzerDE pipeline entry point
- [`analyzeNormalizations()`](https://computationalproteomics.github.io/NormalyzerDE/reference/analyzeNormalizations.md)
  : Calculate measures for normalization results
- [`generatePlots()`](https://computationalproteomics.github.io/NormalyzerDE/reference/generatePlots.md)
  : Generates a number of visualizations for the performance measures
  calculated for the normalized matrices. These contain both general
  measures and direct comparisons for different normalization
  approaches.
- [`generateStatsReport()`](https://computationalproteomics.github.io/NormalyzerDE/reference/generateStatsReport.md)
  : Generate full output report plot document. Plots p-value histograms
  for each contrast in the NormalyzerStatistics instance and writes
  these to a PDF report.
- [`writeNormalizedDatasets()`](https://computationalproteomics.github.io/NormalyzerDE/reference/writeNormalizedDatasets.md)
  : Write normalization matrices to file

## Results objects

- [`NormalyzerResults()`](https://computationalproteomics.github.io/NormalyzerDE/reference/NormalyzerResults.md)
  : Representation of the results from performing normalization over a
  dataset
- [`NormalyzerEvaluationResults()`](https://computationalproteomics.github.io/NormalyzerDE/reference/NormalyzerEvaluationResults.md)
  : Representation of evaluation results by calculating performance
  measures for an an NormalyzerResults instance
- [`NormalyzerStatistics()`](https://computationalproteomics.github.io/NormalyzerDE/reference/NormalyzerStatistics.md)
  : Class representing a dataset for statistical processing in
  NormalyzerDE

## Input and setup

- [`loadData()`](https://computationalproteomics.github.io/NormalyzerDE/reference/loadData.md)
  : Load raw data into dataframe
- [`loadDesign()`](https://computationalproteomics.github.io/NormalyzerDE/reference/loadDesign.md)
  : Load raw design into dataframe
- [`setupJobDir()`](https://computationalproteomics.github.io/NormalyzerDE/reference/setupJobDir.md)
  : Create empty directory for run
- [`setupRawDataObject()`](https://computationalproteomics.github.io/NormalyzerDE/reference/setupRawDataObject.md)
  : Prepare SummarizedExperiment object for raw data to be normalized
  containing data, design and annotation information
- [`setupRawContrastObject()`](https://computationalproteomics.github.io/NormalyzerDE/reference/setupRawContrastObject.md)
  : Prepare SummarizedExperiment object for statistics data
- [`reduceTechnicalReplicates()`](https://computationalproteomics.github.io/NormalyzerDE/reference/reduceTechnicalReplicates.md)
  : Remove technical replicates from data and design
- [`getVerifiedNormalyzerObject()`](https://computationalproteomics.github.io/NormalyzerDE/reference/getVerifiedNormalyzerObject.md)
  : Verify that input data is in correct format, and if so, return a
  generated NormalyzerDE data object from that input data
- [`generateAnnotatedMatrix()`](https://computationalproteomics.github.io/NormalyzerDE/reference/generateAnnotatedMatrix.md)
  : Generate an annotated data frame from statistics object
- [`getRTNormalizedMatrix()`](https://computationalproteomics.github.io/NormalyzerDE/reference/getRTNormalizedMatrix.md)
  : Perform RT-segmented normalization by performing the supplied
  normalization over retention-time sliced data
- [`getSmoothedRTNormalizedMatrix()`](https://computationalproteomics.github.io/NormalyzerDE/reference/getSmoothedRTNormalizedMatrix.md)
  : Generate multiple RT time-window normalized matrices where one is
  shifted. Merge them using a specified method (mean or median) and
  return the result.

## Normalization methods

- [`normMethods()`](https://computationalproteomics.github.io/NormalyzerDE/reference/normMethods.md)
  : Perform normalizations on Normalyzer dataset
- [`globalIntensityNormalization()`](https://computationalproteomics.github.io/NormalyzerDE/reference/globalIntensityNormalization.md)
  : The normalization divides the intensity of each variable in a sample
  with the sum of intensities of all variables in the sample and
  multiplies with the median of sum of intensities of all variables in
  all samples. The normalized data is then log2-transformed.
- [`meanNormalization()`](https://computationalproteomics.github.io/NormalyzerDE/reference/meanNormalization.md)
  : Intensity of each variable in a given sample is divided by the mean
  of sum of intensities of all variables in the sample and then
  multiplied by the mean of sum of intensities of all variables in all
  samples. The normalized data is then transformed to log2.
- [`medianNormalization()`](https://computationalproteomics.github.io/NormalyzerDE/reference/medianNormalization.md)
  : Intensity of each variable in a given sample is divided by the
  median of intensities of all variables in the sample and then
  multiplied by the mean of median of sum of intensities of all
  variables in all samples. The normalized data is then
  log2-transformed.
- [`performQuantileNormalization()`](https://computationalproteomics.github.io/NormalyzerDE/reference/performQuantileNormalization.md)
  : Quantile normalization is performed by the function
  "normalize.quantiles" from the package preprocessCore.
- [`performVSNNormalization()`](https://computationalproteomics.github.io/NormalyzerDE/reference/performVSNNormalization.md)
  : Log2 transformed data is normalized using the function "justvsn"
  from the VSN package.
- [`performCyclicLoessNormalization()`](https://computationalproteomics.github.io/NormalyzerDE/reference/performCyclicLoessNormalization.md)
  : Cyclic Loess normalization
- [`performSMADNormalization()`](https://computationalproteomics.github.io/NormalyzerDE/reference/performSMADNormalization.md)
  : Median absolute deviation normalization Normalization subtracts the
  median and divides the data by the median absolute deviation (MAD).
- [`performGlobalRLRNormalization()`](https://computationalproteomics.github.io/NormalyzerDE/reference/performGlobalRLRNormalization.md)
  : Global linear regression normalization

## Differential expression

- [`calculateContrasts()`](https://computationalproteomics.github.io/NormalyzerDE/reference/calculateContrasts.md)
  : Performs statistical comparisons between the supplied conditions. It
  uses the design matrix and data matrix in the supplied
  NormalyzerStatistics object. A column is supplied specifying which of
  the columns in the design matrix that is used for deciding the sample
  groups. The comparisons vector specifies which pairwise comparisons
  between condition levels that are to be calculated.
