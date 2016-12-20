#' Debugging utility used to output "DEBUG" together with output
#' @param ... Input to print
#' @param debug Append text DEBUG in front of statement, default is TRUE
myprint <- function(..., debug=T) {
    if (debug) {
        print(paste("DEBUG: ", ...))
    }
    else {
        print(paste(...))
    }
}

#' Debugging utility used to output information about variable together with
#'  its content
#' @param variable Target variable to print information about
#' @param debug Append text DEBUG in front of statement, default is TRUE
varprint <- function(variable, debug=T) {
    varName <- deparse(substitute(variable))
    
    if (debug) {
        cat(varName, variable, sprintf("[%s]", typeof(variable)), sprintf("[%s elements]", length(variable)), "\n")
    }
    else {
        cat("DEBUG: ", varName, variable, sprintf("[%s]", typeof(variable)), sprintf("[%s elements]", length(variable)), "\n")
    }
}


#' Create empty directory for run if not already present
#' 
#' @param jobName Name of the run.
#' @param outputDir Path to directory where to create the output directory.
#' @return Path to newly created directory.
#' @export
setupJobDir <- function(jobName, outputDir) {
    
    print("DEBUG: Setup job dir called")
    
    varprint(jobName)
    varprint(outputDir)
    
    if (is.null(outputDir)) {
        jobDir <- paste(getwd(), "/", jobName[1], sep="")
    }
    else {
        jobDir <- paste(outputDir, "/", jobName[1], sep="")
    }
    createDirectory(jobDir)
    
    varprint(jobDir)
    
    jobDir
}

#' Create directory, or return error if already present
#' 
#' @param targetPath Path where to attempt to create direcotry
createDirectory <- function(targetPath) {
    
    if (file.exists(targetPath)) {
        error_handler <- "Directory already exists"
        class(error_handler) <- "try-error"
        if (inherits(error_handler, "try-error")) {
            return(error_handler)
        }
        stop("Directory already exists")
    } 
    else {
        dir.create(targetPath)
    }
}