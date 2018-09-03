#' Small example dataset used to demonstrate code consistency in testing and
#' as dummy data in the vignette.
#' 
#' @format A data frame containing annotation and expression data
#' @keywords internal
"example_data"

#' Same data as in "example_data", but omitting the annotation meaning that it
#' only contains the expression data.
#'  
#' @format A data frame containing expression data
#' @keywords internal
"example_data_only_values"

#' Design matrix corresponding to the small example datasets.
#' 
#' @format A design matrix corresponding to the dataset "example_data"
#' @keywords internal
"example_design"

#' Same data as in "example_data", but normalized and ready for statistical
#' processing.
#' 
#' @format A normalized data frame ready for statistical processing
#' @keywords internal
"example_stat_data"

#' NormalyzerDE results object used internally for testing
#' 
#' @format An instance of the class NormalyzerResults
#' @keywords internal
"regression_test_nr"

#' SummarizedExperiment object prepared with design-matrix, data-matrix
#' and annotation columns loaded
#' 
#' @format An instance of the class SummarizedExperiment
#' @keywords internal
"example_summarized_experiment"
