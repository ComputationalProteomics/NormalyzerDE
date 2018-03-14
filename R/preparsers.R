


legacyNormalyzerToDesign <- function(legacyMatrixFp, sep="\t") {

    fullRawData <- as.matrix(
           utils::read.table(legacyMatrixFp, header=FALSE, sep="\t",
                             stringsAsFactors=FALSE, quote="", comment.char=""))

    allAnnotGroups <- as.numeric(fullRawData[1,])
    parsedDf <- fullRawData[-1:-2,]
    colnames(parsedDf) <- fullRawData[2,]
    
    sampleColumns <- fullRawData[2, allAnnotGroups > 0]
    
    designDf <- data.frame(sample=sampleColumns, group=allAnnotGroups[allAnnotGroups > 0], stringsAsFactors=F)
    designDf
}


proteoisToNormalyzerLegacy <- function(proteiosFp, sep="\t") {
    
    values_df <- read.csv(proteiosFp, skip=9, sep=sep, header=F, na.strings="null", stringsAsFactors=FALSE, comment.char="", quote="")
    
    con <- file(proteiosFp, open="r")
    head_lines <- readLines(con, 9)
    close(con)
    
    sample_line_nbr <- 6
    annot_line_nbr <- 9
    
    sample_fields <- strsplit(trimws(gsub("Sample name", "", head_lines[[sample_line_nbr]])), "\t")[[1]]
    annot_fields <- strsplit(trimws(head_lines[[annot_line_nbr]]), "\t")[[1]]

    header_names <- c(annot_fields, sample_fields)
    values_df <- rbind(header_names, values_df)
    full_m <- as.matrix(values_df)
    full_m
}

proteiosToNormalyzer <- function(proteiosFp, sep="\t") {
    values_df <- as.matrix(
        utils::read.table(inputPath, 
        header=FALSE, 
        sep="\t", 
        stringsAsFactors=FALSE,
        quote=""))
    full_m <- as.matrix(values_df)
    full_m
}

maxQuantToNormalyzer <- function(maxQuantFp, protLevel, sep="\t") {
    
    pep_intensity_pattern <- "Intensity\\."
    
    if (!protLevel) {
        annot_cols <- c("Sequence", "Mass", "Proteins", "Leading.razor.protein", "PEP", "Charges")
    }
    else {
        annot_cols <- c("Protein.IDs", "Majority.protein.IDs", "Fasta.headers")
    }
    
    full_df <- read.csv(maxQuantFp, sep=sep, stringsAsFactors=FALSE, comment.char="", quote="", header=T)
    cnames <- colnames(full_df)
    intensity_cols <- cnames[which(grepl(pep_intensity_pattern, cnames))]
    
    header_names <- c(annot_cols, intensity_cols)
    header_names_trimmed <- gsub("Intensity.", "", header_names)

    values_df <- full_df[, c(annot_cols, intensity_cols)]
    # colnames(target_df) <- all_cols_trimmed
    raw_df <- rbind(header_names_trimmed, values_df)
    raw_m <- as.matrix(raw_df)
    
    return(raw_m)
}

