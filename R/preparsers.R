


legacyNormalyzerToDesign <- function(legacyMatrixFp, sep="\t") {

    fullRawData <- as.matrix(
           utils::read.table(legacyMatrixFp, header=FALSE, sep="\t",
                             stringsAsFactors=FALSE, quote="", comment.char=""))

    allAnnotGroups <- as.numeric(fullRawData[1,])
    parsedDf <- fullRawData[-1:-2,]
    colnames(parsedDf) <- fullRawData[2,]
    
    sampleColumns <- fullRawData[2, allAnnotGroups > 0]
    
    designDf <- data.frame(sample=sampleColumns, group=allAnnotGroups[allAnnotGroups > 0], stringsAsFactors=F)
    
    return(designDf)
}


proteoisToNormalyzer <- function(proteiosFp, sep="\t") {
    
    values_df <- read.csv(proteiosFp, skip=9, sep=sep, header=F, na.strings="null")
    
    head_lines <- readLines(file(proteiosFp, open="r"), 9)
    
    sample_line_nbr <- 6
    annot_line_nbr <- 9
    
    sample_fields <- strsplit(trimws(gsub("Sample name", "", head_lines[[sample_line_nbr]])), "\t")[[1]]
    annot_fields <- strsplit(trimws(head_lines[[annot_line_nbr]]), "\t")[[1]]
    
    colnames(values_df) <- c(annot_fields, sample_fields)

    return(values_df)    
}

maxQuantToNormalyzer <- function(maxQuantFp, sep="\t") {
    
    pep_intensity_pattern <- "Intensity\\."
    annot_cols <- c("Sequence", "Mass", "Proteins", "Leading.razor.protein", "PEP", "Charges")
    
    full_df <- read.csv(maxQuantFp, sep=sep)
    intensity_cols <- cnames[which(grepl(pep_intensity_pattern, cnames))]
    
    all_cols <- c(annot_cols, intensity_cols)
    all_cols_trimmed <- gsub("Intensity.", "", all_cols)

    target_df <- full_df[, c(annot_cols, intensity_cols)]
    colnames(target_df) <- all_cols_trimmed
    
    return(target_df)
}

