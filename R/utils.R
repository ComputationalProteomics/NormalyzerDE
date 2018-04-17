#' Create empty directory for run if not already present
#' 
#' @param jobName Name of the run.
#' @param outputDir Path to directory where to create the output directory.
#' @return Path to newly created directory.
#' @export
setupJobDir <- function(jobName, outputDir) {
    
    if (is.null(outputDir)) {
        jobDir <- paste(getwd(), "/", jobName[1], sep="")
    }
    else {
        jobDir <- paste(outputDir, "/", jobName[1], sep="")
    }
    createDirectory(jobDir)
    
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
        dir.create(targetPath, recursive=TRUE)
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
    sampleValues <- header[which(as.numeric(header) > 0)]
    headerCounts <- table(sampleValues)
    
    nonReplicatedSamples <- names(headerCounts[which(headerCounts == 1)])
    nonReplicatedSamples
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL, single_legend=FALSE) {
    
    # library(grid)
    
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots == 1) {
        print(plots[[1]])
        
    } 
    else {
        # Set up the page
        grid::grid.newpage()
        grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow(layout), ncol(layout))))
        
        # Make each plot, in the correct location
        for (i in 1:numPlots) {
            # Get the i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
            
            print(plots[[i]], vp = grid::viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
        }
    }
}


grid_arrange_shared_legend <- function(plots, ncol = length(list(...)), nrow = 1, position = c("bottom", "right"), plot_info="[no extra info]") {
    
    # library(grid)
    # library(gridExtra)
    
    # plots <- list(...)
    position <- match.arg(position)
    g <- ggplot2::ggplotGrob(plots[[1]] + ggplot2::theme(legend.position = position))$grobs
    legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    lwidth <- sum(legend$width)
    gl <- lapply(plots, function(x) x + ggplot2::theme(legend.position="none"))
    gl <- c(gl, ncol = ncol, nrow = nrow)
    
    combined <- switch(position,
                       "bottom" = gridExtra::arrangeGrob(do.call(gridExtra::arrangeGrob, gl),
                                              legend,
                                              ncol = 1,
                                              heights = grid::unit.c(grid::unit(1, "npc") - lheight, lheight)),
                       "right" = gridExtra::arrangeGrob(do.call(gridExtra::arrangeGrob, gl),
                                             legend,
                                             ncol = 2,
                                             widths = grid::unit.c(grid::unit(1, "npc") - lwidth, lwidth)),
                       top = grid::textGrob("test", vjust = 1, gp = grid::gpar(fontface = "bold", cex = 1.5)))

    grid::grid.draw(combined)
    grid::grid.text(plot_info, x=0.1, y=0.9, just="left")
}

get_datestamp_string <- function() {
    format(Sys.time(), "%Y%m%d_%H%M%S")
}


getIndexList <- function(targetVector) {
    
    indexList <- list()
    uniq_vals <- unique(targetVector)
    for (val in uniq_vals) {
        indexList[[toString(val)]] <- which(targetVector == val)
    }
    indexList
}


#' Removes rows from matrix where replicates aren't represented by at least
#' one value
#' 
#' @param dataMatrix Matrix with expression values for entities in replicate 
#'  samples.
#' @param replicateHeader Header showing how samples in matrix are replicated.
#' @return Reduced matrix where rows without any number are excluded.
filterLinesWithEmptySamples <- function(dataMatrix, replicateHeader) {
    
    replicatesHaveDataCtr <- getAllSamplesHaveValuesContrast(dataMatrix, replicateHeader)
    dataMatrix[replicatesHaveDataCtr, ]
}


#' Get contrast vector (TRUE/FALSE-values) indicating whether all samples
#' specified in replicateHeader have at least one non-NA value
#' 
#' @param dataMatrix Matrix with expression values for entities in replicate 
#'  samples.
#' @param replicateHeader Header showing how samples in matrix are replicated.
#' @param minCount Minimum number of non-NA required for sample.
#' @return Contrast vector
getAllSamplesHaveValuesContrast <- function(dataMatrix, replicateHeader, minCount=1) {
    
    replicatesHaveData <- rep(TRUE, nrow(dataMatrix))
    indexList <- getIndexList(replicateHeader)
    
    for (sampleIndex in 1:length(names(indexList))) {
        
        repVal <- names(indexList)[sampleIndex]
        cols <- indexList[[repVal]]
        
        nbrNAperReplicate <- rowSums(is.na(dataMatrix[, cols, drop=FALSE]))
        nbrReplicates <- length(cols)
        nbrNonNA <- nbrReplicates - nbrNAperReplicate
        replicatesHaveData <- (nbrNonNA >= minCount & replicatesHaveData)
    }
    
    replicatesHaveData
}


#' Get contrast vector (TRUE/FALSE-values) indicating whether at least half
#' of values in row have non-NA values
#' 
#' @param dataMatrix Matrix with expression values for entities in replicate 
#'  samples.
#' @param replicateHeader Header showing how samples in matrix are replicated.
#' @return Contrast vector
getLowNALinesContrast <- function(dataMatrix, replicateHeader) {
    
    nbsNAperLine <- rowSums(is.na(dataMatrix))
    enoughNAContrast <- (nbsNAperLine < ncol(dataMatrix) / 2)
    
    enoughNAContrast
}


#' Get contrast vector (TRUE/FALSE-values) indicating whether both at least
#' half values are present, and each sample has at least one non-NA value
#' 
#' @param dataMatrix Matrix with expression values for entities in replicate 
#'  samples.
#' @param replicateHeader Header showing how samples in matrix are replicated.
#' @param var_filter_frac Variance filtering fraction.
#' @param minCount Minimum number of required values present in samples.
#' @return Contrast vector
getRowNAFilterContrast <- function(dataMatrix, 
                                   replicateHeader, 
                                   var_filter_frac=NULL,
                                   minCount=1) {
    
    lowNALinesContrast <- getLowNALinesContrast(dataMatrix, replicateHeader)
    samplesHaveValuesContrast <- getAllSamplesHaveValuesContrast(dataMatrix, replicateHeader, minCount=minCount)
    
    if (!is.null(var_filter_frac)) {
        samplesPassVarThres <- getVarFilteredContrast(dataMatrix, var_filter_frac)
    }
    else {
        samplesPassVarThres <- rep(TRUE, nrow(dataMatrix))
    }
    
    return (lowNALinesContrast & samplesHaveValuesContrast & samplesPassVarThres)
}


getVarFilteredContrast <- function(logMatrix, var_filter_frac) {
  
    feature_vars <- apply(logMatrix, 1, function(x) {stats::var(x, na.rm=TRUE)})
    frac_thres <- sort(feature_vars)[(length(feature_vars) - 1)  * var_filter_frac + 1]
    feature_vars > frac_thres
}


