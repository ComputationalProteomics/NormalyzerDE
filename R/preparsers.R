
proteiosToNormalyzer <- function(proteiosFp, sep="\t") {
    
    valuesDf <- as.matrix(
        utils::read.table(proteiosFp, 
        header=FALSE, 
        sep="\t", 
        stringsAsFactors=FALSE,
        quote=""))
    fullMat <- as.matrix(valuesDf)
    fullMat
}

maxQuantToNormalyzer <- function(maxQuantFp, protLevel, sep="\t") {
    
    pepIntensityPattern <- "Intensity\\."
    
    if (!protLevel) {
        annotCols <- c("Sequence", "Mass", "Proteins", "Leading.razor.protein", "PEP", "Charges")
        matrixType <- "peptide.txt"
    }
    else {
        annotCols <- c("Protein.IDs", "Majority.protein.IDs", "Fasta.headers")
        matrixType <- "proteinGroups.txt"
    }
    
    fullDf <- utils::read.csv(maxQuantFp, sep=sep, stringsAsFactors=FALSE, comment.char="", quote="", header=TRUE)
    cnames <- colnames(fullDf)
    intensityCols <- cnames[grepl(pepIntensityPattern, cnames)]
    
    headerNames <- c(annotCols, intensityCols)
    headerNamesTrimmed <- gsub("Intensity.", "", headerNames)

    if (!(all(annotCols %in% colnames(fullDf)))) {
        stop(
            "Performing MaxQuant processing for the ", matrixType, " matrix\n", 
            "Didn't find all of the following expected columns (with spaces instead of dots): \n",
            paste(annotCols, collapse=", "),
            "\nColumns in the input data:\n",
            paste(colnames(fullDf), collapse=", ")
        )
    }
    
    valuesDf <- fullDf[, c(annotCols, intensityCols)]
    rawDf <- rbind(headerNamesTrimmed, valuesDf)
    rawMat <- as.matrix(rawDf)
    
    return(rawMat)
}

