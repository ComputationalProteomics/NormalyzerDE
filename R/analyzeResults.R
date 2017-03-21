#' Calculate measures for normalization results
#' 
#' @param nr Normalyzer results object with calculated results.
#' @param name Name of the ongoing processing run.
#'
#' @return Normalyzer results with attached evaluation results object.
#' @export
analyzeNormalizations <- function(nr, name, comparisons=NULL, categorical_anova=FALSE) {
    
    nds <- nr@nds
    nr@ner <- setupNormalizationEvaluationObject(nr, 
                                                 comparisons=comparisons,
                                                 categorical_anova=categorical_anova)
    print("Analysis finished. Next, preparing plots and report.")
    nr
}

#' Pearson and Spearman correlation calculations for methods and samples
#' ! Needs some further investigations (why? speed or something else?)
#' 
#' @param nr Normalyzer results object with calculated results.
#' @param ner Normalyzer evaluation object.
#' @return Normalyzer evaluation object with attached evaluation results.
calculateCorrelations <- function(nr, ner) {
    
    methodlist <- getNormalizationMatrices(nr)
    filterED <- nr@nds@sampleReplicateGroups
    
    graphics::par(mfrow=c(1,1))
    avgpercorsum <- list()
    avgspecorsum <- list()
    corsum <- vector()
    
    for (i in 1:length(methodlist)) {
        percorsum <- vector()
        specorsum <- vector()
        
        flag1 <- 1
        datastore <- as.matrix(methodlist[[i]])
        un <- unique(filterED)
        
        for (uq in 1:length(un)) {
            dt <- as.matrix(datastore[, which(filterED == un[uq])])
            class(dt) <- "numeric"
            percor <- stats::cor(dt, use="pairwise.complete.obs", method="pearson")
            spercor <- stats::cor(dt, use="pairwise.complete.obs", method="spearman")
            
            for (rn in 1:(ncol(dt) - 1)) {
                percorsum <- c(percorsum, percor[rn, -(1:rn)])
                specorsum <- c(specorsum, spercor[rn, -(1:rn)])
            }
        }
        
        avgpercorsum[[i]] <- percorsum
        avgspecorsum[[i]] <- specorsum
    }
    
    ner@avgpercorsum <- avgpercorsum
    ner@avgspecorsum <- avgspecorsum
    ner
}


#' Setup normalization evaluation object 
#' 
#' @param nr Normalyzer results object to be evaluated
#' @return Normalization evaluation object
setupNormalizationEvaluationObject <- function(nr, comparisons=NULL, categorical_anova=FALSE) {
    
    ner <- NormalizationEvaluationResults()
    ner <- calculateCV(ner, nr)
    
    singleRepRun <- nr@nds@singleReplicateRun 
    
    if (!singleRepRun) {
        ner <- calculateMAD(ner, nr)
        ner <- calculateAvgVar(ner, nr)
        ner <- calculateSignificanceMeasures(ner, nr, categorical_anova=categorical_anova)
        
        if (!is.null(comparisons)) {
            ner <- calculatePairwiseComparisons(ner, nr, comparisons)
        }
    }
    
    ner <- calculateCorrelations(nr, ner)
    ner
}

