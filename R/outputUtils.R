#' Write normalization matrices to file
#' 
#' Outputs each of the normalized datasets to the specified directory.
#' 
#' @param nr Results object.
#' @param jobdir Path to output directory.
#' @param includePairwiseComparisons Include p-values for pairwise comparisons.
#' @param includeCvCol Include CV column in output.
#' @param includeAnovaP Include ANOVA p-value in output.
#' @param normSuffix String used to name output together with normalization names.
#' @param rawdataName Name of output raw data file.
#' @return None
#' @export
#' @examples
#' data(example_data)
#' data(example_design)
#' normObj <- getVerifiedNormalyzerObject("job_name", example_design, example_data)
#' normResults <- normMethods(normObj)
#' normResultsWithEval <- analyzeNormalizations(normResults)
#' outputDir <- tempdir()
#' writeNormalizedDatasets(normResultsWithEval, outputDir)
writeNormalizedDatasets <- function(nr, jobdir, includePairwiseComparisons=FALSE, 
                                    includeCvCol=FALSE, includeAnovaP=FALSE,
                                    normSuffix="-normalized.txt",
                                    rawdataName="submitted_rawdata.txt") {
    
    nds <- nds(nr)
    ner <- ner(nr)

    methodnames <- getUsedMethodNames(nr)
    methodlist <- getNormalizationMatrices(nr)
    annotationColumns <- annotationValues(nds)
    
    for (sampleIndex in seq_along(methodnames)) {
        
        currentMethod <- methodnames[sampleIndex]
        filePath <- paste(jobdir, "/", currentMethod, normSuffix, sep="")
        outputTable <- cbind(annotationColumns, methodlist[[sampleIndex]])

        if (includeAnovaP) {
            anovaP <- anovaP(ner)[,sampleIndex]
            
            if (nrow(outputTable) != length(anovaP)) {
                stop("Table row count: ", nrow(outputTable), 
                     " must match p-value vector length for anova: ", 
                     length(anovaP))
            }
            
            outputTable <- cbind(outputTable, anovaP=anovaP)
        }
        
        if (includePairwiseComparisons) {
            
            ner <- ner(nr)
            for (comp in names(pairwiseComps(ner))) {
                
                compColP <- pairwiseComps(ner)[[comp]][, sampleIndex]
                compColFdr <- pairwiseCompsFdr(ner)[[comp]][, sampleIndex]
                
                newColnames <- c(
                    colnames(outputTable), 
                    paste("comp", comp, "p", sep="_"), 
                    paste("comp", comp, "fdr", sep="_")
                )
                outputTable <- cbind(outputTable, compColP, compColFdr)
                colnames(outputTable) <- newColnames
            }
        }
        
        if (includeCvCol) {
            cvCol <- featureCVPerMethod(ner)[, sampleIndex]
            outputTable <- cbind(outputTable, CV=cvCol)
        }

        utils::write.table(
            outputTable, file=filePath, sep="\t", row.names=FALSE, quote=FALSE)
    }
    
    rawFilePath <- paste(jobdir, "/", rawdataName, sep="")
    rawOutputTable <- cbind(annotationColumns, filterrawdata(nds))
    
    utils::write.table(
        rawOutputTable, 
        file=rawFilePath, 
        sep="\t", 
        row.names=FALSE, 
        quote=FALSE
    )
}