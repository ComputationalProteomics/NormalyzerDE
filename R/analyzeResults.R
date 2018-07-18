#' Calculate measures for normalization results
#' 
#' @param nr Normalyzer results object with calculated results.
#' @param comparisons Target sample contrasts to run.
#' @param categoricalAnova ANOVA can be categorical or numerical.
#' @return Normalyzer results with attached evaluation results object.
#' @export
#' @examples
#' data(example_data)
#' data(example_design)
#' normObj <- getVerifiedNormalyzerObject("job_name", example_design, example_data)
#' normResults <- normMethods(normObj)
#' normResultsWithEval <- analyzeNormalizations(normResults)
analyzeNormalizations <- function(nr, comparisons=NULL, categoricalAnova=FALSE) {
    
    nds <- nr@nds
    nr@ner <- setupNormalizationEvaluationObject(
        nr, 
        comparisons=comparisons,
        categoricalAnova=categoricalAnova)

    nr
}

#' Pearson and Spearman correlation calculations for methods and samples
#' Calculates internal correlation per condition
#' 
#' @param nr Normalyzer results object with calculated results.
#' @param ner Normalyzer evaluation object.
#' @return ner Normalyzer evaluation object with attached evaluation results.
#' @keywords internal
calculateCorrelations <- function(nr, ner) {
    
    methodlist <- getNormalizationMatrices(nr)
    allReplicateGroups <- nr@nds@sampleReplicateGroups
    sampleGroupsWithReplicates <- nr@nds@samplesGroupsWithReplicates

    ner@avgpercorsum <- calculateSummarizedCorrelationVector(
        methodlist, allReplicateGroups, sampleGroupsWithReplicates, "pearson")
    ner@avgspecorsum <- calculateSummarizedCorrelationVector(
        methodlist, allReplicateGroups, sampleGroupsWithReplicates, "spearman")
    ner
}

calculateSummarizedCorrelationVector <- function(
    methodlist, allReplicateGroups, sampleGroupsWithReplicates, corrType) {
    
    validCorrTypes <- c("pearson", "spearman")
    if (!corrType %in% validCorrTypes) {
        stop(paste(
            "Unknown correlation type:", 
            corrType, 
            "valid are:", 
            paste(validCorrTypes, collapse=", ")))
    }
    
    avgCorSum <- list()
    
    for (i in seq_len(length(methodlist))) {
        
        corSum <- vector()
        methodData <- as.matrix(methodlist[[i]])
        
        for (groupNbr in seq_len(length(sampleGroupsWithReplicates))) {
            
            specificReplicateVals <- as.matrix(
                methodData[, which(allReplicateGroups == sampleGroupsWithReplicates[groupNbr])])
            class(specificReplicateVals) <- "numeric"
            corVals <- stats::cor(
                specificReplicateVals , 
                use="pairwise.complete.obs", 
                method=corrType)
            # spearCor <- stats::cor(
            #     specificReplicateVals , 
            #     use="pairwise.complete.obs", 
            #     method="spearman")
            
            for (index in seq_len(ncol(specificReplicateVals) - 1)) {
                corSum <- c(corSum, corVals[index, -(seq_len(index))])
            }
        }
        
        avgCorSum[[i]] <- corSum
    }
    
    avgCorSum
}


#' Setup normalization evaluation object 
#' 
#' @param nr Normalyzer results object to be evaluated
#' @param comparisons Group comparisons to make
#' @param categoricalAnova If comparison should be made numerically or groupwise
#' @return Normalization evaluation object
#' @keywords internal
setupNormalizationEvaluationObject <- function(nr, 
                                               comparisons=NULL, 
                                               categoricalAnova=FALSE) {
    
    ner <- NormalyzerEvaluationResults()
    ner <- calculateCV(ner, nr)
    
    singleRepRun <- nr@nds@singleReplicateRun
    
    if (!singleRepRun) {
        ner <- calculateMAD(ner, nr)
        ner <- calculateAvgVar(ner, nr)
        ner <- calculateSignificanceMeasures(
            ner, 
            nr, 
            categoricalAnova=categoricalAnova)
        
    }
    
    ner <- calculateCorrelations(nr, ner)
    ner
}

