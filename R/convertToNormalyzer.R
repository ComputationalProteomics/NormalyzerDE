

proteoisToNormalyzer <- function(proteios_fp) {
    
    values_df <- read.csv(proteios_fp, skip=9, sep="\t", header=F, na.strings="null")
    
    head_lines <- readLines(file(proteios_fp, open="r"), 9)
    
    sample_line_nbr <- 6
    annot_line_nbr <- 9
    
    sample_fields <- strsplit(trimws(gsub("Sample name", "", head_lines[[sample_line_nbr]])), "\t")[[1]]
    annot_fields <- strsplit(trimws(head_lines[[annot_line_nbr]]), "\t")[[1]]
    
    colnames(values_df) <- c(annot_fields, sample_fields)

    return(values_df)    
}

maxQuantToNormalyzer <- function(maxquant_fp) {
    
    pep_intensity_pattern <- "Intensity\\."
    annot_cols <- c("Sequence", "Mass", "Proteins", "Leading.razor.protein", "PEP", "Charges")
    
    full_df <- read.csv(maxquant_fp, sep="\t")
    intensity_cols <- cnames[which(grepl(pep_intensity_pattern, cnames))]
    
    all_cols <- c(annot_cols, intensity_cols)
    all_cols_trimmed <- gsub("Intensity.", "", all_cols)

    target_df <- full_df[, c(annot_cols, intensity_cols)]
    colnames(target_df) <- all_cols_trimmed
    
    return(target_df)
}

