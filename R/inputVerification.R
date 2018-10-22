#' Load raw data into dataframe
#' 
#' General function which allows specifying different types of input data
#' including "proteios", "maxquantpep" (peptide output from MaxQuant) and 
#' "maxquantprot" (protein output from MaxQuant) formats.
#' 
#' @param dataPath File path to design matrix.
#' @param inputFormat If input is given in standard NormalyzerDE format, 
#' Proteios format or in MaxQuant protein or peptide format
#' @param zeroToNA Automatically convert zeroes to NA values
#' @return rawData Raw data loaded into data frame
#' @examples \dontrun{
#' df <- loadData("data.tsv")
#' }
loadData <- function(dataPath, inputFormat="default", zeroToNA=FALSE) {
    
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
    
    if (zeroToNA) {
        rawData[rawData == "0"] <- NA
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
    rawData <- loadData(dataPath, inputFormat=inputFormat, zeroToNA=zeroToNA)
    rdf <- rawData[2:nrow(rawData), ]
    colnames(rdf) <- rawData[1, ]
    
    verifyDesignMatrix(rdf, rawDesign, sampleColName)
        
    rawDesign[[sampleColName]] <- as.character(rawDesign[[sampleColName]])
    
    sdf <- rdf[, as.character(rawDesign[[sampleColName]])]
    adf <- rdf[, !(colnames(rdf) %in% as.character(rawDesign[[sampleColName]]))]
    
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
    adf <- fullDf[, !(colnames(fullDf) %in% as.character(designDf[[sampleColName]]))]
    
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
        quiet=FALSE
    ) {

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
    ), stringsAsFactors=FALSE))

    verifyDesignMatrix(dataMatrix, designMatrix, sampleCol)
    verifyValidNumbers(dataMatrix, groups, quiet=quiet)
    
    repSortedRawData <- getReplicateSortedData(dataMatrix, groups)
    processedRawData <- preprocessData(repSortedRawData)
    
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
        rawData=dataMatrix,
        annotationData=annotationMatrix,
        sampleNameCol=sampleCol,
        groupNameCol=groupCol,
        quiet=quiet
    )
    
    nds
}


filterOnlyNARows <- function(summarizedExp) {

    dataMatrix <- SummarizedExperiment::assay(summarizedExp)
    
    nonFullNAContr <- rowSums(is.na(SummarizedExperiment::assay(summarizedExp))) != ncol(summarizedExp)
    if (length(which(!nonFullNAContr)) > 0) {
        warning(
            length(which(!nonFullNAContr)), 
            " entries with only NA values found and were omitted")
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
            stop("Please provide a valid input file.")
        },
        
        warning=function(w) {
            message("An issue was encountered when attempting to load:", 
                    inputPath,
                    "Warning:",
                    w)
            stop("Please investigate the warning and provide a valid input file.")
        }
    )
    
    rawData
}


#' Verify that input fields conform to the expected formats
#' 
#' @param rawDataOnly Dataframe with input data.
#' @param groups Condition levels for comparisons.
#' @return Parsed rawdata where 0 values are replaced with NA
#' @keywords internal
verifyValidNumbers <- function(rawDataOnly, groups, quiet=FALSE) {
    
    # Fields expected to contain numbers in decimal or scientific notation, or containing NA or null
    validPatterns <- c("\\d+(\\.\\d+)?", "NA", "\"NA\"", "null", "\\d+(\\.\\d+)?[eE]([\\+\\-])?\\d+$")
    
    regexPattern <- sprintf("^(%s)$", paste(validPatterns, collapse="|"))
    matches <- grep(regexPattern, rawDataOnly, perl=TRUE, ignore.case=TRUE)
    nonMatches <- stats::na.omit(rawDataOnly[-matches])
    
    if (length(nonMatches) > 0) {
        errorString <- paste(
            "Invalid values encountered in input data.",
            "Only valid data is numeric and NA- or na-fields",
            "Invalid fields: ",
            paste(unique(nonMatches), collapse=" "),
            "Aborting...",
            sep="\n"
        )
        
        stop(errorString)
    }
    
    zeroRegexPattern <- c("^0$")
    zeroMatches <- grep(zeroRegexPattern, rawDataOnly, perl=TRUE)
    if (length(zeroMatches) > 0) {
        errorString <- paste(
            "Encountered zeroes in data. Must be replaced with NA before processing.",
            "This can be done automatically setting the zeroToNA-flag option to TRUE.",
            sep="\n"
        )
        stop(errorString)
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
            paste(setdiff(designColnames, colnames(fullMatrix)), collapse=", "), " \n",
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




#' Replace 0 values with NA in input data
#' 
#' @param dataMatrix Matrix with raw data.
#' @return Parsed rawdata where 0 values are replaced with NA
#' @keywords internal
preprocessData <- function(dataMatrix) {

    dataMatrix[dataMatrix == 0] <- NA
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
            "Following samples does not have replicates:",
            paste(nonReplicatedSamples, collapse=" "),
            "By default this is not allowed.",
            "You can force limited processing of non-replicated data by",
            "specifying the \"requireReplicates\" option",
            sep="\n")
        
        if (requireReplicates) {
            stop(errorString)
        }
        else {
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
            "Found less than two distinct samples. Following was found:",
            paste(distinctSamples, collapse=" "),
            "For full processing two samples are required",
            "You can force limited processing by turning of the", 
            "\"requireReplicates\" option",
            sep="\n")
        
        if (requireReplicates) {
            stop(errorString)
        }
        else {
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


#' #' Setup Normalyzer dataset from given raw data
#' #' 
#' #' @param jobName Name of ongoing run.
#' #' @param designMatrix Dataframe containing condition matrix.
#' #' @param fullRawMatrix Dataframe with unparsed input data.
#' #' @param sampleNameCol Name of column in design matrix containing sample names.
#' #' @param groupNameCol Name of column in design matrix contaning conditions.
#' #' @return Data object representing loaded data.
#' #' @keywords internal
#' generateNormalyzerDataset <- function(jobName, designMatrix, fullRawMatrix, sampleNameCol, groupNameCol, quiet=FALSE) {
#'     
#'     nds <- NormalyzerDataset(
#'         jobName=jobName, 
#'         rawData=fullRawMatrix, 
#'         designMatrix=designMatrix,
#'         sampleNameCol=sampleNameCol, 
#'         groupNameCol=groupNameCol)
#'     nds <- setupValues(nds, quiet=quiet)
#'     nds
#' }

