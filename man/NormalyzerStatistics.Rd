% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NormalyzerStatistics.R
\name{NormalyzerStatistics}
\alias{NormalyzerStatistics}
\title{Class representing a dataset for statistical processing in NormalyzerDE}
\usage{
NormalyzerStatistics(experimentObj, logTrans = FALSE)

NormalyzerStatistics(experimentObj, logTrans = FALSE)
}
\arguments{
\item{experimentObj}{Instance of SummarizedExperiment containing matrix
and design information as column data}

\item{logTrans}{Whether the input data should be log transformed}
}
\value{
nds Generated NormalyzerStatistics instance
}
\description{
Is initialized with an annotation matrix, a data matrix and a design
data frame. This object can subsequently be processed to generate statistical
values and in turn used to write a full matrix with additional statistical
information as well as a graphical report of the comparisons.
}
\section{Slots}{

\describe{
\item{\code{annotMat}}{Matrix containing annotation information}

\item{\code{dataMat}}{Matrix containing (normalized) expression data}

\item{\code{filteredDataMat}}{Filtered matrix with low-count rows removed}

\item{\code{designDf}}{Data frame containing design conditions}

\item{\code{filteringContrast}}{Vector showing which entries are filtered 
(due to low count)}

\item{\code{pairwiseCompsP}}{List with P-values for pairwise comparisons}

\item{\code{pairwiseCompsFdr}}{List with FDR-values for pairwise comparisons}

\item{\code{pairwiseCompsAve}}{List with average expression values}

\item{\code{pairwiseCompsFold}}{List with log2 fold-change values for pairwise 
comparisons}

\item{\code{contrasts}}{Spot for saving vector of last used contrasts}

\item{\code{condCol}}{Column containing last used conditions}

\item{\code{batchCol}}{Column containing last used batch conditions}
}}

\examples{
data(example_stat_summarized_experiment)
nst <- NormalyzerStatistics(example_stat_summarized_experiment)
}
