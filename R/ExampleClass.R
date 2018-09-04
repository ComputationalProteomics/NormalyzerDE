#' @rdname ExampleClass
#' @export
ExampleClass <- setClass("ExampleClass", representation(value_sum = "numeric", value_product = "numeric"))

#' Constructor for ExampleClass
#' 
#' @param first_number First numeric value
#' @param second_number Second numeric value
#' 
#' @return instance ExampleClass instance
#' 
#' @aliases ExampleClass
#'
#' @docType class
#'
#' @examples
#' instance <- ExampleClass(2,3)
#' @rdname ExampleClass
#' @importFrom utils packageVersion
#' @export
setMethod("ExampleClass", definition = function(first_number, second_number) { 
    
              value_sum <- first_number + second_number
              value_prod <- first_number * second_number
              
              instance <- new("ExampleClass",
                value_sum=value_sum,
                value_product=value_prod
              ) 

              instance
          })

#setGeneric("ExampleClass", function(first_number, second_number) { standardGeneric("ExampleClass") })
