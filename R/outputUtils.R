#' Write normalization matrices to file
#' 
#' @param nr Normalyzer results 
#' @return None
writeNormalizedDatasets <- function(nr, jobdir, include_pvals=FALSE, include_pairwise_comparisons=FALSE, include_cv_col=FALSE, include_anova_p=FALSE) {
    
    nds <- nr@nds
    ner <- nr@ner
    
    methodnames <- getUsedMethodNames(nr)
    methodlist <- getNormalizationMatrices(nr)
    annotationColumns <- nds@annotationValues
    
    for (sampleIndex in 1:length(methodnames)) {
        
        currentMethod <- methodnames[sampleIndex]
        filePath <- paste(jobdir, "/", currentMethod, "-normalized.txt", sep="")
        outputTable <- cbind(annotationColumns, methodlist[[sampleIndex]])

        if (include_pvals) {
            
            anova_col <- ner@anovaFDRWithNA[,sampleIndex]
            kw_col <- ner@krusWalFDRWithNA[,sampleIndex]
            
            if (nrow(outputTable) != length(anova_col) || nrow(outputTable) != length(kw_col)) {
                stop(paste("Table row count:", nrow(outputTable), "must match p-value vector lengths for anova:", 
                           length(anova_col), "and and kruskal wallis:", length(kw_col)))
            }
            outputTable <- cbind(outputTable, anova=anova_col, kruskal_wallis=kw_col)
        }
        
        if (include_anova_p) {
            anova_p <- ner@anova_p[,sampleIndex]
            
            if (nrow(outputTable) != length(anova_p)) {
                stop(paste("Table row count:", nrow(outputTable), "must match p-value vector length for anova: ", length(anova_p)))
            }
            
            outputTable <- cbind(outputTable, anova_p=anova_p)
        }
        
        if (include_pairwise_comparisons) {
            
            for (comp in names(nr@ner@pairwise_comps)) {
                
                comp_col_p <- ner@pairwise_comps[[comp]][,sampleIndex]
                comp_col_fdr <- ner@pairwise_comps_fdr[[comp]][,sampleIndex]
                
                new_colnames <- c(colnames(outputTable), paste("comp", comp, "p", sep="_"), paste("comp", comp, "fdr", sep="_"))
                outputTable <- cbind(outputTable, comp_col_p, comp_col_fdr)
                colnames(outputTable) <- new_colnames
            }
        }
        
        if (include_cv_col) {
            cv_col <- ner@featureCVPerMethod[,sampleIndex]
            outputTable <- cbind(outputTable, CV=cv_col)
        }

        utils::write.table(outputTable, file=filePath, sep="\t", row.names=FALSE, quote=FALSE)
    }
    
    if (!all(is.na(nr@houseKeepingVars))) {
        hkVarsName <- "housekeeping-variables.tsv"
        hkFilePath <- paste(jobdir, "/", hkVarsName, sep="")
        
        utils::write.table(file=hkFilePath, nr@houseKeepingVars, sep="\t", 
                           row.names=FALSE, quote=FALSE)
    }
    
    rawdata_name <- "submitted_rawdata.tsv"
    rawFilePath <- paste(jobdir, "/", rawdata_name, sep="")
    rawOutputTable <- cbind(annotationColumns, nds@filterrawdata)
    
    utils::write.table(rawOutputTable, file=rawFilePath, sep="\t", row.names=FALSE, quote=FALSE)
}