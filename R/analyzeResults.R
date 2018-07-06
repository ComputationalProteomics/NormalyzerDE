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
#' normObj <- getVerifiedNormalyzerObject("data.tsv", "job_name", "design.tsv")
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
#' ! Needs some further investigations (why? speed or something else?)
#' 
#' @param nr Normalyzer results object with calculated results.
#' @param ner Normalyzer evaluation object.
#' @return Normalyzer evaluation object with attached evaluation results.
calculateCorrelations <- function(nr, ner) {
    
    methodlist <- getNormalizationMatrices(nr)
    allReplicateGroups <- nr@nds@sampleReplicateGroups

    graphics::par(mfrow=c(1,1))
    avgpercorsum <- list()
    avgspecorsum <- list()
    corsum <- vector()
    
    for (i in 1:length(methodlist)) {
        pearCorSum <- vector()
        spearCorSum <- vector()
        
        # flag1 <- 1
        datastore <- as.matrix(methodlist[[i]])
        sampleGroupsWithReplicates <- nr@nds@samplesGroupsWithReplicates
        # allUniqueRepGroups <- unique(allReplicateGroups)
        
        # browser()
        
        for (groupNbr in 1:length(sampleGroupsWithReplicates)) {
            
            # print(paste("Calculating for group:", groupNbr))
            # 
            # browser()
            
            dt <- as.matrix(datastore[, which(allReplicateGroups == sampleGroupsWithReplicates[groupNbr])])
            class(dt) <- "numeric"
            pearCor <- stats::cor(dt, use="pairwise.complete.obs", method="pearson")
            spearCor <- stats::cor(dt, use="pairwise.complete.obs", method="spearman")
            
            for (rn in 1:(ncol(dt) - 1)) {
                pearCorSum <- c(pearCorSum, pearCor[rn, -(1:rn)])
                spearCorSum <- c(spearCorSum, spearCor[rn, -(1:rn)])
            }
        }
        
        # browser()

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

