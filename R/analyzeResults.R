#' Calculate measures for normalization results
#' 
#' @param nr Normalyzer results object with calculated results.
#' @param name Name of the ongoing processing run.
#'
#' @return Normalyzer results with attached evaluation results object.
#' @export
analyzeNormalizations <- function(nr, name) {
    
    print("DEBUG: analyzeAndPlot entered")

    nds <- nr@nds
    nr@ner <- setupNormalizationEvaluationObject(nr)

    print("Analysis finished. Next, preparing plots and report.")
    
    nr
}

#' Pearson and Spearman correlation calculations for methods and samples
#' ! Needs some further investigations 
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


#' Removes rows from matrix where replicates aren't represented by at least
#' one value
#' 
#' @param dataMatrix Matrix with expression values for entities in replicate 
#'  samples.
#' @param replicateHeader Header showing how samples in matrix are replicated.
#' @return Reduced matrix where rows without any number are excluded.
filterLinesWithEmptySamples <- function(dataMatrix, replicateHeader) {
    
    firstIndices <- getFirstIndicesInVector(replicateHeader)
    lastIndices <- getFirstIndicesInVector(replicateHeader, reverse=TRUE)
    
    replicatesHaveData <- rep(TRUE, nrow(dataMatrix))
    
    for (i in 1:length(firstIndices)) {
        
        firstIndex <- firstIndices[i]
        lastIndex <- lastIndices[i]
        
        nbrNAperReplicate <- rowSums(is.na(dataMatrix[, firstIndex:lastIndex]))
        nbrReplicates <- lastIndex - firstIndex + 1
        replicatesHaveData <- (nbrNAperReplicate < nbrReplicates & replicatesHaveData)
    }
    
    dataMatrix[replicatesHaveData, ]
}


#' Setup normalization evaluation object 
#' 
#' @param nr Normalyzer results object to be evaluated
#' @return Normalization evaluation object
setupNormalizationEvaluationObject <- function(nr) {
    
    ner <- NormalizationEvaluationResults()
    ner <- calculateCV(ner, nr)
    
    singleRepRun <- nr@nds@singleReplicateRun 
    
    if (!singleRepRun) {
        ner <- calculateMAD(ner, nr)
        ner <- calculateAvgVar(ner, nr)
        ner <- calculateSignificanceMeasures(ner, nr)
    }
    
    ner <- calculateCorrelations(nr, ner)
    ner
}

