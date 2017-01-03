#' Debugging utility used to output "DEBUG" together with output
#' @param ... Input to print
#' @param debug Append text DEBUG in front of statement, default is TRUE
#' @return None
myprint <- function(..., debug=TRUE) {
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
#' @return None
varprint <- function(variable, debug=TRUE) {
    varName <- deparse(substitute(variable))
    
    if (debug) {
        cat(varName, variable, sprintf("[%s]", typeof(variable)), 
            sprintf("[%s elements]", length(variable)), "\n")
    }
    else {
        cat("DEBUG: ", varName, variable, sprintf("[%s]", typeof(variable)), 
            sprintf("[%s elements]", length(variable)), "\n")
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
#' @param targetPath Path where to attempt to create directory
#' @return None
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

#' Get number of seconds between two Sys.time() objects
#' 
#' @param start Start-time object
#' @param end End-time object
#' @return None
elapsedSecondsBetweenSystimes <- function(start, end) {
    
    startSecond <- strtoi(format(start, "%s"))
    endSecond <- strtoi(format(end, "%s"))
    elapsed <- end - start
    elapsed
}

#' Returns samples present only once in Normalyzer header
#' 
#' @param normalyzerDf Normalyzer data frame
#' @return Vector with names of non-replicated samples
getNonReplicatedFromDf <- function(normalyzerDf) {
    
    header <- normalyzerDf[1,]
    sampleValues <- header[which(header != "0")]
    headerCounts <- table(sampleValues)
    
    nonReplicatedSamples <- names(headerCounts[which(headerCounts == 1)])
    nonReplicatedSamples
}