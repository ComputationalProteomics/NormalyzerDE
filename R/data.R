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

#' SummarizedExperiment object prepared with design-matrix, data-matrix
#' and annotation columns loaded for raw data
#' 
#' @format An instance of the class SummarizedExperiment
#' @keywords internal
"example_summarized_experiment"

#' SummarizedExperiment object prepared with design-matrix, data-matrix
#' and annotation columns for normalized data
#' 
#' @format An instance of the class SummarizedExperiment with stats data
#' @keywords internal
"example_stat_summarized_experiment"

#' Full raw NormalyzerDE matrix used for internal testing
#' 
#' @format A data table ready for analysis in NormalyzerDE
#' @keywords internal
"example_wide_data"

#' Design matrix belonging together with example_wide_data. Used for
#' internal testing.
#' 
#' @format A design table ready for analysis in NormalyzerDE
#' @keywords internal
"example_wide_design"

