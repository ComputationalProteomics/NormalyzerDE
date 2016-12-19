

#' S4 class to represent normalization evaluations
#' 
#' @slot avgcvmem TODO: Look into
#' @slot avgmadmem TODO: Look into
#' @slot avgvarmem TODO: Look into
#' @slot avgcvmempdiff TODO: Look into
#' @slot avgmadmempdiff TODO: Look into
#' @slot avgvarmempdiff TODO: Look into
#' @slot nonsiganfdrlist TODO: Look into
#' @slot nonsiganfdrlistcvpdiff TODO: Look into
#' @slot anfdr TODO: Look into
#' @slot kwfdr TODO: Look into
#' @slot avgpercorsum TODO: Look into
#' @slot avgspecorsum TODO: Look into
#' @export
NormalizationEvaluationResults <- setClass("NormalizationEvaluationResults",
                                            slots=c(
                                                avgcvmem = "matrix",
                                                avgmadmem = "matrix",
                                                avgvarmem = "matrix",
                                                avgcvmempdiff = "numeric",
                                                avgmadmempdiff = "numeric",
                                                avgvarmempdiff = "numeric",
                                                nonsiganfdrlist = "numeric",
                                                nonsiganfdrlistcvpdiff = "numeric",
                                                anfdr = "matrix",
                                                kwfdr = "matrix",
                                                
                                                avgpercorsum = "list",
                                                avgspecorsum = "list"
                                            ),
                                            prototype=prototype())
                                            
