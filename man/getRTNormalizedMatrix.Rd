% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/higherOrderNormMethods.R
\name{getRTNormalizedMatrix}
\alias{getRTNormalizedMatrix}
\title{Perform RT-segmented normalization by performing the supplied normalization
over retention-time sliced data}
\usage{
getRTNormalizedMatrix(
  rawMatrix,
  retentionTimes,
  normMethod,
  stepSizeMinutes = 1,
  windowMinCount = 100,
  offset = 0,
  noLogTransform = FALSE
)
}
\arguments{
\item{rawMatrix}{Target matrix to be normalized}

\item{retentionTimes}{Vector of retention times corresponding to rawMatrix}

\item{normMethod}{The normalization method to apply to the time windows}

\item{stepSizeMinutes}{Size of windows to be normalized}

\item{windowMinCount}{Minimum number of values for window to not be expanded.}

\item{offset}{Whether time window should shifted half step size}

\item{noLogTransform}{Don't log-transform the data}
}
\value{
Normalized matrix
}
\description{
The function orders the retention times and steps through them using the
supplied step size (in minutes). If smaller than a fixed lower boundary 
the window is expanded to ensure a minimum amount of data in each 
normalization step. An offset can be specified which can be used to perform
multiple RT-segmentations with partial overlapping windows.
}
\examples{
data(example_data_small)
data(example_design_small)
data(example_data_only_values)
dataMat <- example_data_only_values
retentionTimes <- as.numeric(example_data[, "Average.RT"])
performCyclicLoessNormalization <- function(rawMatrix) {
    log2Matrix <- log2(rawMatrix)
    normMatrix <- limma::normalizeCyclicLoess(log2Matrix, method="fast")
    colnames(normMatrix) <- colnames(rawMatrix)
    normMatrix
}
rtNormMat <- getRTNormalizedMatrix(dataMat, retentionTimes, 
performCyclicLoessNormalization, stepSizeMinutes=1, windowMinCount=100)
}
