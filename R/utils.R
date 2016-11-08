myprint <- function(..., debug=T) {
    if (debug) {
        print(paste("DEBUG: ", ...))
    }
    else {
        print(paste(...))
    }
}

returnOne <- function() {
    1
}

varprint <- function(variable, debug=T) {
    varName <- deparse(substitute(variable))
    
    if (debug) {
        cat(varName, variable, sprintf("[%s]", typeof(variable)), sprintf("[%s elements]", length(variable)), "\n")
    }
    else {
        cat("DEBUG: ", varName, variable, sprintf("[%s]", typeof(variable)), sprintf("[%s elements]", length(variable)), "\n")
    }
}
