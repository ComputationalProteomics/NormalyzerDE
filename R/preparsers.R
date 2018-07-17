
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
    }
    else {
        annotCols <- c("Protein.IDs", "Majority.protein.IDs", "Fasta.headers")
    }
    
    fullDf <- utils::read.csv(maxQuantFp, sep=sep, stringsAsFactors=FALSE, comment.char="", quote="", header=TRUE)
    cnames <- colnames(fullDf)
    intensityCols <- cnames[which(grepl(pepIntensityPattern, cnames))]
    
    headerNames <- c(annotCols, intensityCols)
    headerNamesTrimmed <- gsub("Intensity.", "", headerNames)

    valuesDf <- fullDf[, c(annotCols, intensityCols)]
    rawDf <- rbind(headerNamesTrimmed, valuesDf)
    rawMat <- as.matrix(rawDf)
    
    return(rawMat)
}

