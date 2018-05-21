#' Write normalization matrices to file
#' 
#' @param nr Results object.
#' @param jobdir Path to output directory.
#' @param includePairwiseComparisons Include p-values for pairwise comparisons.
#' @param includeCvCol Include CV column in output.
#' @param includeAnovaP Include ANOVA p-value in output.
#' @param normSuffix String used to name output together with normalization names.
#' @param rawdataName Name of output raw data file.
#' @param hkVarsName Name out output housekeeping variables file.
#' @return None
#' @export
#' @examples
#' normObj <- getVerifiedNormalyzerObject("data.tsv", "jobName", "design.tsv")
#' normResults <- normMethods(normObj)
#' normResultsWithEval <- analyzeNormalizations(normObj)
#' writeNormalizedDatasets(normResultsWithEval, "path/to/output")
writeNormalizedDatasets <- function(nr, jobdir, includePairwiseComparisons=FALSE, 
                                    includeCvCol=FALSE, includeAnovaP=FALSE,
                                    normSuffix="-normalized.txt",
                                    rawdataName="submitted_rawdata.txt",
                                    hkVarsName="housekeeping-variables.tsv") {
    
    nds <- nr@nds
    ner <- nr@ner

    methodnames <- getUsedMethodNames(nr)
    methodlist <- getNormalizationMatrices(nr)
    annotationColumns <- nds@annotationValues
    
    for (sampleIndex in 1:length(methodnames)) {
        
        currentMethod <- methodnames[sampleIndex]
        filePath <- paste(jobdir, "/", currentMethod, normSuffix, sep="")
        outputTable <- cbind(annotationColumns, methodlist[[sampleIndex]])

        if (includeAnovaP) {
            anovaP <- ner@anovaP[,sampleIndex]
            
            if (nrow(outputTable) != length(anovaP)) {
                stop(paste("Table row count:", nrow(outputTable), 
                           "must match p-value vector length for anova: ", 
                           length(anovaP)))
            }
            
            outputTable <- cbind(outputTable, anovaP=anovaP)
        }
        
        if (includePairwiseComparisons) {
            
            for (comp in names(nr@ner@pairwiseComps)) {
                
                compColP <- ner@pairwiseComps[[comp]][, sampleIndex]
                compColFdr <- ner@pairwiseCompsFdr[[comp]][, sampleIndex]
                
                newColnames <- c(colnames(outputTable), paste("comp", comp, "p", sep="_"), paste("comp", comp, "fdr", sep="_"))
                outputTable <- cbind(outputTable, compColP, compColFdr)
                colnames(outputTable) <- newColnames
            }
        }
        
        if (includeCvCol) {
            cvCol <- ner@featureCVPerMethod[, sampleIndex]
            outputTable <- cbind(outputTable, CV=cvCol)
        }

        utils::write.table(outputTable, file=filePath, sep="\t", row.names=FALSE, quote=FALSE)
    }
    
    rawFilePath <- paste(jobdir, "/", rawdataName, sep="")
    rawOutputTable <- cbind(annotationColumns, nds@filterrawdata)
    
    utils::write.table(rawOutputTable, file=rawFilePath, sep="\t", row.names=FALSE, quote=FALSE)
}