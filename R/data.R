#' Small example dataset used to demonstrate code consistency in testing and
#' as dummy data in the vignette.
#' 
#' @format A data frame containing annotation and expression data
"example_data"

#' Same data as in "example_data", but omitting the annotation meaning that it
#' only contains the expression data.
#'  
#' @format A data frame containing expression data
"example_data_only_values"

#' Design matrix corresponding to the small example datasets.
#' 
#' @format A design matrix corresponding to the dataset "example_data"
"example_design"

#' Same data as in "example_data", but normalized and ready for statistical
#' processing.
#' 
#' @format A normalized data frame ready for statistical processing
"example_stat_data"

#' NormalyzerDE results object used internally for testing
#' 
#' @format An instance of the class NormalyzerResults
"regression_test_nr"
