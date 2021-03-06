#' Load raw data into dataframe
#' 
#' General function which allows specifying different types of input data
#' including "proteios", "maxquantpep" (peptide output from MaxQuant) and 
#' "maxquantprot" (protein output from MaxQuant) formats.
#' 
#' @param dataPath File path to design matrix.
#' @param inputFormat If input is given in standard NormalyzerDE format, 
#' Proteios format or in MaxQuant protein or peptide format
#' @return rawData Raw data loaded into data frame
#' @examples \dontrun{
#' df <- loadData("data.tsv")
#' }
loadData <- function(dataPath, inputFormat="default") {
    
    if (inputFormat == "default") {
        rawData <- loadRawDataFromFile(dataPath)
    } else if (inputFormat == "proteios") {
        rawData <- proteiosToNormalyzer(dataPath)
    } else if (inputFormat == "maxquantpep") {
        rawData <- maxQuantToNormalyzer(dataPath, protLevel=FALSE)
    } else if (inputFormat == "maxquantprot") {
        rawData <- maxQuantToNormalyzer(dataPath, protLevel=TRUE)
    } else {
        valids <- c("default", "proteios", "maxquantpep", "maxquantprot")
        stop("Unknown inputFormat: ", 
             inputFormat, 
             " valids are: ", 
             paste(valids, collapse=", ")
        )
    }
    
    rawData
}

#' Load raw design into dataframe
#' 
#' Takes a design path, loads the matrix and ensures that the sample column
#' is in character format and that the group column is in factor format.
#' 
#' @param designPath File path to design matrix.
#' @param sampleCol Column name for column containing sample names.
#' @param groupCol Column name for column containing condition levels.
#' @return designMatrix Design data loaded into data frame
#' @examples \dontrun{
#' df <- loadDesign("design.tsv")
#' }
loadDesign <- function(designPath, sampleCol="sample", groupCol="group") {
    
    designMatrix <- utils::read.table(
        designPath, 
        sep="\t", 
        stringsAsFactors=FALSE, 
        header=TRUE, 
        comment.char=""
    )
    
    if (!(sampleCol %in% colnames(designMatrix)) || !(groupCol %in% colnames(designMatrix))) {
        stop("Both sampleCol value and groupCol value must be present in the design matrix header.",
             "\nsampleCol: ", sampleCol, "\ngroupCol: ", groupCol,
             "\nDesign matrix header: ", paste(colnames(designMatrix), collapse=", "))
    }
    
    designMatrix[, sampleCol] <- as.character(designMatrix[, sampleCol])
    designMatrix[, groupCol] <- as.factor(as.character(designMatrix[, groupCol]))
    designMatrix
}

#' Prepare SummarizedExperiment object for raw data to be normalized containing 
#' data, design and annotation information
#' 
#' @param dataPath File path to data matrix.
#' @param designPath File path to design matrix.
#' @param inputFormat Type of matrix for data, can be either 'default',
#'   'proteios', 'maxquantprot' or 'maxquantpep'
#' @param zeroToNA If TRUE zeroes in the data is automatically converted to
#'   NA values
#' @param sampleColName Column name for column containing sample names
#' @param groupColName Column name for column containing condition levels
#' @return experimentObj SummarizedExperiment object loaded with the data
#' @export
#' @examples 
#' data_path <- system.file(package="NormalyzerDE", "extdata", "tiny_data.tsv")
#' design_path <- system.file(package="NormalyzerDE", "extdata", "tiny_design.tsv")
#' df <- setupRawDataObject(data_path, design_path)
setupRawDataObject <- function(dataPath, designPath, inputFormat="default", zeroToNA=FALSE, 
                               sampleColName="sample", groupColName="group") {
    
    rawDesign <- loadDesign(designPath, sampleCol=sampleColName, groupCol=groupColName)
    rawData <- loadData(dataPath, inputFormat=inputFormat)
    rdf <- rawData[2:nrow(rawData), ]
    colnames(rdf) <- rawData[1, ]
    
    verifyDesignMatrix(rdf, rawDesign, sampleColName)
        
    rawDesign[[sampleColName]] <- as.character(rawDesign[[sampleColName]])
    
    sdf <- rdf[, as.character(rawDesign[[sampleColName]])]
    adf <- rdf[, !(colnames(rdf) %in% as.character(rawDesign[[sampleColName]])), drop=FALSE]
    
    experimentObj <- SummarizedExperiment::SummarizedExperiment(
        assays=list(raw=as.matrix(sdf)),
        rowData=adf,
        colData=rawDesign,
        metadata=list(sample=sampleColName, group=groupColName)
    )
    experimentObj
}

#' Prepare SummarizedExperiment object for statistics data
#' 
#' @param dataPath Path to raw data matrix
#' @param designPath Path to design matrix
#' @param sampleColName Name for column in design matrix containing sample names
#' @return experimentObj Prepared instance of SummarizedExperiment
#' @export
#' @examples 
#' data_path <- system.file(package="NormalyzerDE", "extdata", "tiny_data.tsv")
#' design_path <- system.file(package="NormalyzerDE", "extdata", "tiny_design.tsv")
#' sumExpObj <- setupRawContrastObject(data_path, design_path, "sample")
setupRawContrastObject <- function(dataPath, designPath, sampleColName) {
    
    fullDf <- utils::read.csv(
        dataPath, 
        sep="\t", 
        stringsAsFactors=FALSE, 
        quote="", 
        comment.char="",
        check.names=FALSE
    )
    
    designDf <- utils::read.csv(
        designPath, 
        sep="\t",
        stringsAsFactors=FALSE,
        quote="",
        comment.char="",
        check.names=FALSE
    )
    
    verifyDesignMatrix(fullDf, designDf, sampleColName)
    
    sdf <- fullDf[, designDf[[sampleColName]]]
    adf <- fullDf[, !(colnames(fullDf) %in% as.character(designDf[[sampleColName]])), drop=FALSE]
    
    experimentObj <- SummarizedExperiment::SummarizedExperiment(
        assays=list(raw=as.matrix(sdf)),
        colData=designDf,
        rowData=adf,
        metadata=list(sample=sampleColName)
    )
    experimentObj
}

#' Verify that input data is in correct format, and if so, return a generated
#'  NormalyzerDE data object from that input data
#' 
#' This function performs a number of checks on the input data and provides
#' informative error messages if the data isn't fulfilling the required format.
#' Checks include verifying that the design matrix matches to the data matrix,
#' that the data matrix contains valid numbers and that samples have enough
#' values for analysis
#' 
#' @param jobName Name of ongoing run.
#' @param summarizedExp Summarized experiment input object
#' @param threshold Minimum number of features.
#' @param omitSamples Automatically omit invalid samples from analysis.
#' @param requireReplicates Require there to be at least to samples per
#'        condition
#' @param quiet Don't print output messages during processing 
#' @param noLogTransform Don't log-transform the provided data
#' @param tinyRunThres If less features in run, a limited run is performed
#' 
#' @return Normalyzer data object representing verified input data.
#' @export
#' @examples
#' data(example_summarized_experiment)
#' normObj <- getVerifiedNormalyzerObject("job_name", example_summarized_experiment)
getVerifiedNormalyzerObject <- function(
        jobName, 
        summarizedExp,
        threshold=15, 
        omitSamples=FALSE,
        requireReplicates=TRUE,
        quiet=FALSE,
        noLogTransform=FALSE,
        tinyRunThres=50
    ) {

    SummarizedExperiment::assay(summarizedExp) <- preprocessData(SummarizedExperiment::assay(summarizedExp), quiet=quiet)
    summarizedExp <- filterOnlyNARows(summarizedExp)
    
    # TODO: The getter seems to not be available, check later if temporary issue
    groupCol <- summarizedExp@metadata$group
    sampleCol <- summarizedExp@metadata$sample
    designMatrix <- as.data.frame(SummarizedExperiment::colData(summarizedExp))
    
    if (!groupCol %in% colnames(designMatrix)) {
        stop("Given groupCol: '", groupCol, "' was not present among design matrix columns")
    }

    groups <- SummarizedExperiment::colData(summarizedExp)[[groupCol]]
    samples <- SummarizedExperiment::colData(summarizedExp)[[sampleCol]]
    dataMatrix <- SummarizedExperiment::assay(summarizedExp)
    annotationMatrix <- as.matrix(data.frame(lapply(
        SummarizedExperiment::rowData(summarizedExp),
        as.character
    ), stringsAsFactors=FALSE, check.names=FALSE))

    verifyDesignMatrix(dataMatrix, designMatrix, sampleCol)
    verifyValidNumbers(dataMatrix, groups, noLogTransform=noLogTransform, quiet=quiet)
    
    processedRawData <- getReplicateSortedData(dataMatrix, groups)
    # repSortedRawData <- getReplicateSortedData(dataMatrix, groups)
    # processedRawData <- preprocessData(repSortedRawData, quiet=quiet)
    
    lowCountSampleFiltered <- getLowCountSampleFiltered(
        processedRawData, 
        groups,
        threshold=threshold, 
        stopIfTooFew=!omitSamples
    )
    
    designMatrix <- designMatrix[designMatrix[[sampleCol]] %in% colnames(lowCountSampleFiltered), ]
    
    # If no samples left after omitting, stop
    verifyMultipleSamplesPresent(
        lowCountSampleFiltered, 
        groups, 
        requireReplicates=requireReplicates, 
        quiet=quiet
    )
    
    validateSampleReplication(
        lowCountSampleFiltered, 
        groups, 
        requireReplicates=requireReplicates, 
        quiet=quiet
    )
    
    nds <- NormalyzerDataset(
        jobName=jobName,
        designMatrix=designMatrix,
        rawData=processedRawData,
        annotationData=annotationMatrix,
        sampleNameCol=sampleCol,
        groupNameCol=groupCol,
        tinyRunThres=tinyRunThres,
        quiet=quiet
    )
    
    nds
}


filterOnlyNARows <- function(summarizedExp) {

    dataMatrix <- SummarizedExperiment::assay(summarizedExp)
    
    nonFullNAContr <- rowSums(is.na(SummarizedExperiment::assay(summarizedExp))) != ncol(summarizedExp)
    if (length(which(!nonFullNAContr)) > 0) {
        message(length(which(!nonFullNAContr)), " entries with only NA values omitted")
        summarizedExp <- summarizedExp[nonFullNAContr, ]
    }
    
    summarizedExp
}

#' Try reading raw Normalyzer matrix from provided filepath
#' 
#' @param inputPath Path to Normalyzer data.
#' @return Table containing raw data from input file.
#' @keywords internal
loadRawDataFromFile <- function(inputPath) {
    
    tryCatch(
        rawData <- as.matrix(utils::read.table(inputPath, 
                                               header=FALSE, 
                                               sep="\t", 
                                               stringsAsFactors=FALSE,
                                               quote="",
                                               comment.char="")),
        error=function(e) {
            message("Error encountered for input file:", inputPath, ", error: ", e)
            stop("Please provide a valid input file.\n")
        },
        
        warning=function(w) {
            message("An issue was encountered when attempting to load:", 
                    inputPath,
                    "Warning:",
                    w)
            stop("Please investigate the warning and provide a valid input file.\n")
        }
    )
    
    rawData
}


#' Verify that input fields conform to the expected formats
#' 
#' @param rawDataOnly Dataframe with input data.
#' @param groups Condition levels for comparisons.
#' @return None
#' @keywords internal
verifyValidNumbers <- function(rawDataOnly, groups, noLogTransform=FALSE, quiet=FALSE) {
    
    # Fields expected to contain numbers in decimal or scientific notation, or containing NA or null
    validPatterns <- c("\\d+(\\.\\d+)?", "NA", "\"NA\"", "null", "\\d+(\\.\\d+)?[eE]([\\+\\-])?\\d+$", "")
    
    regexPattern <- sprintf("^(%s)$", paste(validPatterns, collapse="|"))
    nonMatchIndices <- grep(regexPattern, rawDataOnly, perl=TRUE, ignore.case=TRUE, invert = TRUE)
    naIndices <- which(is.na(rawDataOnly))
    invalidNonNAIndices <- nonMatchIndices[!nonMatchIndices %in% naIndices]
    rowsWithIssues <- unique((invalidNonNAIndices-1) %% nrow(rawDataOnly) + 1)
    
    if (length(invalidNonNAIndices) > 0) {
        errorString <- paste(
            "Invalid values encountered in input data.",
            "Only valid data is numeric (dot-decimal, not comma) and NA- or na-fields",
            "Invalid fields: ",
            paste(rawDataOnly[unique(invalidNonNAIndices)], collapse=" "),
            "These were encountered for row numbers (showing max 10 first):",
            paste(utils::head(rowsWithIssues, 10), collapse=", "),
            "Content of the first row with issues: ",
            paste(rawDataOnly[rowsWithIssues[1], ], collapse=", "),
            "Aborting...",
            sep="\n"
        )
        
        stop(errorString)
    }
    
    if (!noLogTransform) {
        belowOnePattern <- c("^0.\\d+$")
        belowOneMatches <- grep(belowOnePattern, rawDataOnly, perl=TRUE)
        rowsWithIssues <- unique((belowOneMatches-1) %% nrow(rawDataOnly) + 1)
        firstIssueRow <- rawDataOnly[rowsWithIssues[1], ]

        if (length(belowOneMatches) > 0) {
            errorString <- paste(
                "Encountered below-one values in raw data. As the data is log-transformed ",
                "during processing this will lead to negative values which in turn will ",
                "crash processing. Consider using the 'noLogTransform' option if your data ",
                "already is normally distributed or scaling all values if appropriate.", 
                "First ten row numbers where this issue was encountered, excluding header row: ",
                paste(utils::head(rowsWithIssues, 10), collapse=", "),
                paste0("Content of first row (row ", rowsWithIssues[1], ") with issues: "),
                paste(firstIssueRow, collapse=", "),
                sep="\n"
            )
            stop(errorString)
        }
    }

    if (!quiet) message("Input data checked. All fields are valid.")
}


#' Verify that design matrix setup matches the data matrix
#' 
#' @param fullMatrix Dataframe with input data.
#' @param designMatrix Dataframe with design setup.
#' @param sampleCol Column in design matrix containing sample IDs.
#' 
#' @return None
#' @keywords internal
verifySummarizedExperiment <- function(summarizedExp, sampleCol) {
    
    fullMatrix <- cbind(
        data.frame(SummarizedExperiment::rowData(summarizedExp)),
        SummarizedExperiment::assay(summarizedExp)
    )
    
    designMatrix <- data.frame(SummarizedExperiment::colData(summarizedExp))
    
    verifyDesignMatrix(
        fullMatrix,
        designMatrix,
        sampleCol
    )
}


#' Verify that design matrix setup matches the data matrix
#' 
#' @param fullMatrix Dataframe with input data.
#' @param designMatrix Dataframe with design setup.
#' @param sampleCol Column in design matrix containing sample IDs.
#' 
#' @return None
#' @keywords internal
verifyDesignMatrix <- function(fullMatrix, designMatrix, sampleCol) {

    if (!(sampleCol %in% colnames(designMatrix))) {
        stop("Design matrix header must contain sampleCol name. \n", 
             "Provided sampleCol was: ", sampleCol, " \n", 
             "Following header was found in the design matrix: ", paste(colnames(designMatrix), collapse=", "))
    }
    
    designColnames <- designMatrix[, sampleCol]
    
    if (!all(designColnames %in% colnames(fullMatrix))) {
        errorString <- paste(
            "Not all columns present in design matrix are present in data matrix. \n",
            " The following elements were not found in the data matrix header: \n\n",
            paste(base::setdiff(designColnames, colnames(fullMatrix)), collapse=", "), " \n",
            "\n Please carefully check that the column names in your data matrix", 
            "matches the sample column in the design matrix"
        )
        stop(errorString)
    }

    dataMatrix <- fullMatrix[, designMatrix[, sampleCol]]
    dataColumns <- dataMatrix[, designColnames]
    
    if (length(designColnames) != ncol(dataColumns)) {
        errorString <- paste(
            "Number of samples does not match number of selected columns",
            "Found number of columns:",
            ncol(dataColumns),
            "Expected number of columns:",
            length(designColnames),
            "Are all columns in the design matrix present in the data matrix?"
        )
        stop(errorString)
    }
    
    if (length(unique(designColnames)) != length(designColnames)) {
        errorString <- paste(
            "Sample labels must be unique, found: ",
            length(unique(designColnames)),
            "unique, expected:",
            length(designColnames)
        )
        stop(errorString)
    }
}




#' Replace empty values (0 or empty field) with NA in input data
#' 
#' @param dataMatrix Matrix with raw data.
#' @param quiet Don't show diagnostic messages
#' @return Parsed rawdata where 0 values are replaced with NA
#' @keywords internal
preprocessData <- function(dataMatrix, quiet=FALSE) {

    zeroFields <- length(dataMatrix[!is.na(dataMatrix) & dataMatrix == 0])
    emptyFields <- length(dataMatrix[!is.na(dataMatrix) & dataMatrix == ""])
    nullFields <- length(dataMatrix[!is.na(dataMatrix) & dataMatrix == "null"])
    
    if (zeroFields != 0) {
        if (!quiet) {
            message(zeroFields, " fields with '0' were replaced by 'NA'")
        }
        dataMatrix[dataMatrix == 0] <- NA
    }
    
    if (emptyFields != 0) {
        if (!quiet) {
            message(emptyFields, " empty fields were replaced by 'NA'")
        }
        dataMatrix[dataMatrix == ""] <- NA
    }
    
    if (nullFields != 0) {
        if (!quiet) {
            message(nullFields, " 'null' fields were replaced by 'NA'")
        }
        dataMatrix[dataMatrix == "null"] <- NA
    }

    dataMatrix
}

#' Verify that samples contain at least a lowest number of values
#' 
#' @param dataMatrix Dataframe with processed input data.
#' @param groups Vector containing condition levels.
#' @param threshold Lowest number of allowed values in a column.
#' @param stopIfTooFew Abort run if lower than threshold number of values in
#'        column
#' @return None
#' @keywords internal
getLowCountSampleFiltered <- function(dataMatrix, groups, threshold=15, stopIfTooFew=TRUE) {
    
    sampleIndices <- seq_along(groups)
    numberOfValues <- colSums(!is.na(dataMatrix))
    notPassingThreshold <- which(numberOfValues < threshold)
    
    if (length(notPassingThreshold) == length(numberOfValues)) {
        errorString <- paste(
            "None of the samples had enough valid non-NA values",
            "Found number of non-NA values:",
            paste(numberOfValues[notPassingThreshold], collapse=" "),
            "Current threshold: ",
            threshold,
            "You could try lowering the threshold by adjusting \"sampleAbundThres\"",
            "Be aware that this will likely lead to downstream crashes.",
            sep="\n"
        )
        stop(errorString)
    }
    else if (length(notPassingThreshold) > 0) {
        errorString <- paste(
            "Following samples does not contain enough non-NA values:",
            paste(sampleIndices[notPassingThreshold], collapse=" "),
            "Found number of values:",
            paste(numberOfValues[notPassingThreshold], collapse=" "),
            "Current threshold:",
            threshold,
            "You can force processing without this sample by specifying the",
            "option \"omitLowAbundSamples\"",
            sep="\n")
        
        if (stopIfTooFew) {
            stop(errorString)
        }
        else {
            warning(errorString)
        }
    }
    
    if (length(notPassingThreshold > 0)) {
        naSamplesOmittedDf <- dataMatrix[, -sampleIndices[notPassingThreshold]]
    }
    else {
        dataMatrix
    }
}


#' Check whether all samples have replicates
#' 
#' @param dataMatrix Prepared matrix containing expression data.
#' @param groups Vector containing condition levels
#' @param requireReplicates By default stops processing if not all samples
#'  have replicates
#' @return None
#' @keywords internal
validateSampleReplication <- function(dataMatrix, groups, requireReplicates=TRUE, quiet=FALSE) {
    
    headerCounts <- table(groups)
    nonReplicatedSamples <- names(headerCounts[headerCounts == 1])

    if (length(nonReplicatedSamples) > 0) {
        
        errorString <- paste(
            "Following group conditions does not have replicates:",
            paste(nonReplicatedSamples, collapse=" "),
            "By default this is not allowed.",
            "You can force limited processing of non-replicated data by",
            "setting the \"requireReplicates\" option to FALSE",
            sep="\n")
        
        if (requireReplicates) {
            stop(errorString)
        }
        else if (!quiet) {
            warning(errorString)
        }
    }
    else {
        if (!quiet) message("Sample replication check: All samples have replicates")
    }
}


#' Check whether more than one sample is present
#' 
#' @param dataMatrix Prepared dataframe.
#' @param groups Vector containing condition levels
#' @param requireReplicates By default stops processing if not all samples
#'  have replicates
#' @return None
#' @keywords internal
verifyMultipleSamplesPresent <- function(dataMatrix, groups, requireReplicates=TRUE, quiet=FALSE) {
    
    samples <- groups[as.numeric(as.factor(groups)) > 0]
    distinctSamples <- unique(samples)

    if (length(samples) < 2) {
        errorString <- paste(
            "At least two samples are required to run Normalyzer",
            "Here, we found:",
            paste(samples, collapse=" "),
            sep="\n")
        
        stop(errorString)
    }
            
    if (length(distinctSamples) == 1) {
        
        errorString <- paste(
            "Found less than two distinct sample groups. Following was found:",
            paste(distinctSamples, collapse=" "),
            "For full processing two or more sample groups are required",
            "You can force limited processing for one sample group by setting the", 
            "\"requireReplicates\" option to FALSE\n",
            sep="\n")
        
        if (requireReplicates) {
            stop(errorString)
        }
        else if (!quiet) {
            warning(errorString)
        }
    }
    else if (length(distinctSamples) == 0) {
        stop("No replicate groups found. Double check your input file,
             and check that your data haven't been filtered out in preceeding
             input validation steps.")
    }
    else {
        if (!quiet) message("Sample check: More than one sample group found")
    }
}



