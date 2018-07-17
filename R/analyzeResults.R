#' Calculate measures for normalization results
#' 
#' @param nr Normalyzer results object with calculated results.
#' @param comparisons Target sample contrasts to run.
#' @param categoricalAnova ANOVA can be categorical or numerical.
#' @param varFilterFrac Perform variance filtering before tests.
#'
#' @return Normalyzer results with attached evaluation results object.
#' @export
#' @examples
#' data(example_data)
#' data(example_design)
#' normObj <- getVerifiedNormalyzerObject("job_name", example_design, example_data)
#' normResults <- normMethods(normObj)
#' normResultsWithEval <- analyzeNormalizations(normResults)
analyzeNormalizations <- function(nr, 
                                  comparisons=NULL, 
                                  categoricalAnova=FALSE,
                                  varFilterFrac=NULL) {
    
    nds <- nr@nds
    nr@ner <- setupNormalizationEvaluationObject(nr, 
                                                 comparisons=comparisons,
                                                 categoricalAnova=categoricalAnova,
                                                 varFilterFrac=varFilterFrac)

    nr
}

#' Pearson and Spearman correlation calculations for methods and samples
#' Calculates internal correlation per condition
#' 
#' @param nr Normalyzer results object with calculated results.
#' @param ner Normalyzer evaluation object.
#' @return Normalyzer evaluation object with attached evaluation results.
#' @keywords internal
calculateCorrelations <- function(nr, ner) {
    
    methodlist <- getNormalizationMatrices(nr)
    allReplicateGroups <- nr@nds@sampleReplicateGroups

    avgpercorsum <- list()
    avgspecorsum <- list()

    for (i in seq_len(length(methodlist))) {
        
        pearCorSum <- vector()
        spearCorSum <- vector()
        methodData <- as.matrix(methodlist[[i]])
        sampleGroupsWithReplicates <- nr@nds@samplesGroupsWithReplicates

        for (groupNbr in seq_len(length(sampleGroupsWithReplicates))) {
            
            specificReplicateVals <- as.matrix(methodData[, which(allReplicateGroups == sampleGroupsWithReplicates[groupNbr])])
            class(specificReplicateVals) <- "numeric"
            pearCor <- stats::cor(specificReplicateVals , use="pairwise.complete.obs", method="pearson")
            spearCor <- stats::cor(specificReplicateVals , use="pairwise.complete.obs", method="spearman")
            
            for (index in seq_len(ncol(specificReplicateVals) - 1)) {
                pearCorSum <- c(pearCorSum, pearCor[index, -(seq_len(index))])
                spearCorSum <- c(spearCorSum, spearCor[index, -(seq_len(index))])
            }
        }
        
        avgpercorsum[[i]] <- pearCorSum
        avgspecorsum[[i]] <- spearCorSum
    }
    
    ner@avgpercorsum <- avgpercorsum
    ner@avgspecorsum <- avgspecorsum
    ner
}


#' Setup normalization evaluation object 
#' 
#' @param nr Normalyzer results object to be evaluated
#' @param comparisons Group comparisons to make
#' @param categoricalAnova If comparison should be made numerically or groupwise
#' @param varFilterFrac Filtering high-variance samples
#' @return Normalization evaluation object
#' @keywords internal
setupNormalizationEvaluationObject <- function(nr, 
                                               comparisons=NULL, 
                                               categoricalAnova=FALSE,
                                               varFilterFrac=NULL) {
    
    ner <- NormalyzerEvaluationResults()
    ner <- calculateCV(ner, nr)
    
    singleRepRun <- nr@nds@singleReplicateRun
    
    if (!singleRepRun) {
        ner <- calculateMAD(ner, nr)
        ner <- calculateAvgVar(ner, nr)
        ner <- calculateSignificanceMeasures(ner, 
                                             nr, 
                                             categoricalAnova=categoricalAnova,
                                             varFilterFrac=varFilterFrac)
        
    }
    
    ner <- calculateCorrelations(nr, ner)
    ner
}

