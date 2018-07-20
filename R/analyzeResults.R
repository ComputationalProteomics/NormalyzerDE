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
    
    ner <- calculateCorrelations(ner, nr)
    nr@ner <- ner
    
    nr
}


calculateReplicateCV <- function(methodList, sampleReplicateGroups) {
    
    methodCount <- length(methodList)
    numberFeatures <- nrow(methodList[[1]])
    conditionLevels <- length(levels(as.factor(unlist(sampleReplicateGroups))))
    avgCVPerNormAndReplicates <- matrix(nrow=conditionLevels, ncol=methodCount, byrow=TRUE)
    
    for (methodIndex in seq_len(methodCount)) {
        
        processedDataMatrix <- methodList[[methodIndex]]
        featureCondCVs <- matrix(
            nrow=nrow(processedDataMatrix), 
            ncol=length(levels(as.factor(unlist(sampleReplicateGroups)))), 
            byrow=TRUE)
        
        for (i in seq_len(nrow(processedDataMatrix))) {
            
            featureCVs <- RcmdrMisc::numSummary(
                processedDataMatrix[i, ], 
                statistics=c("cv"), 
                groups=unlist(sampleReplicateGroups))
            featureCondCVs[i, ] <- featureCVs$table
        }
        
        # Sum up CV for all features for each condition
        summedCondCVs <- apply(featureCondCVs, 2, mean, na.rm=TRUE)
        avgCVPerNormAndReplicates[, methodIndex] <- summedCondCVs * 100
    }
    
    avgCVPerNormAndReplicates
}

calculateFeatureCV <- function(methodList) {
    
    methodCount <- length(methodList)
    numberFeatures <- nrow(methodList[[1]])
    methodFeatureCVMatrix <- matrix(nrow=numberFeatures, ncol=methodCount)

    for (methodIndex in seq_len(methodCount)) {
        
        processedDataMatrix <- methodList[[methodIndex]]
        # Calculate sample-wise CV
        cv <- function(row) {
            stDev <- stats::sd(row, na.rm=TRUE)
            meanVal <- mean(unlist(row), na.rm=TRUE)
            stDev / mean(meanVal) * 100
        }
        featureCVs <- apply(processedDataMatrix, 1, cv)
        methodFeatureCVMatrix[, methodIndex] <- featureCVs
    }
    
    methodFeatureCVMatrix
}







#' Calculates correlation values between replicates for each condition matrix.
#' Finally returns a list containing the results for all matrices.
#' 
#' @param methodlist List containing normalized matrices for each normalization
#'   method
#' @param allReplicateGroups Vector with condition groups matching the columns
#'   found in the normalization methods
#' @param sampleGroupsWithReplicates Unique vector with condition groups
#'   present in two or more samples
#' @param corrType Type of correlation (Pearson or Spearman)
#' @return avgCorSum
#' @export
#' @examples
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
        
        methodData <- as.matrix(methodlist[[i]])
        corSum <- calculateCorrSum(
            methodData, allReplicateGroups, 
            sampleGroupsWithReplicates, corrType)
        avgCorSum[[i]] <- corSum
    }
    
    avgCorSum
}

#' Calculates internal correlations for each condition having at least two
#' samples and returns a vector with correlation values corresponding to each
#' condition
#' 
#' @param methodData
#' @param allReplicateGroups Full condition header corresponding to data tables
#'   columns
#' @param sampleGroupsWithReplicates Unique conditions where number of
#'   replicates exceeds one
#' @param corrType Type of correlation (Pearson or Spearman)
#' @return corSums
#' @export
#' @examples
calculateCorrSum <- function(methodData, allReplicateGroups, 
                             sampleGroupsWithReplicates, corrType) {
    
    corSums <- vector()
    for (groupNbr in seq_len(length(sampleGroupsWithReplicates))) {
        
        specificReplicateVals <- as.matrix(
            methodData[, which(allReplicateGroups == sampleGroupsWithReplicates[groupNbr])])
        class(specificReplicateVals) <- "numeric"
        corVals <- stats::cor(
            specificReplicateVals , 
            use="pairwise.complete.obs", 
            method=corrType)
        
        for (index in seq_len(ncol(specificReplicateVals) - 1)) {
            corSums <- c(corSums, corVals[index, -(seq_len(index)), drop="FALSE"])
        }
    }
    
    corSums
}

















