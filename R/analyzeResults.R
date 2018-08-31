#' Calculate measures for normalization results
#' 
#' This function prepares an NormalyzerEvaluationResults object containing
#' the evaluation measures CV (coefficient of variance), MAD (median absolute
#' deviation), average variance, significance measures (ANOVA between
#' condition groups) and correlation between replicates.
#' 
#' @param nr Normalyzer results object with calculated results.
#' @param categoricalAnova Whether categorical or numerical (ordered) ANOVA
#' should be calculated.
#' @return Normalyzer results with attached evaluation results object.
#' @export
#' @examples
#' data(example_data)
#' data(example_design)
#' normObj <- getVerifiedNormalyzerObject("job_name", example_design, example_data)
#' normResults <- normMethods(normObj)
#' normResultsWithEval <- analyzeNormalizations(normResults)
analyzeNormalizations <- function(nr, categoricalAnova=FALSE) {
    
    # browser()
    
    nds <- nds(nr)
    ner <- NormalyzerEvaluationResults()
    ner <- calculateCV(ner, nr)
    singleRepRun <- singleReplicateRun(nds)
    
    if (!singleRepRun) {
        ner <- calculateMAD(ner, nr)
        ner <- calculateAvgVar(ner, nr)
        ner <- calculateSignificanceMeasures(
            ner, 
            nr, 
            categoricalAnova=categoricalAnova)
    }
    
    ner <- calculateCorrelations(ner, nr)
    ner(nr) <- ner
    
    nr
}

#' Calculate CV per replicate group and normalization technique
#' 
#' Iterates through each normalization method and calculate average CV values
#' per replicate group.
#' 
#' @param methodList List containing normalized matrices.
#' @param sampleReplicateGroups Condition header.
#' @return avgCVPerNormAndReplicates Matrix with group CVs as rows and
#'   normalization technique as columns
#' @keywords internal
calculateReplicateCV <- function(methodList, sampleReplicateGroups) {
    
    methodCount <- length(methodList)
    numberFeatures <- nrow(methodList[[1]])
    conditionLevels <- length(levels(as.factor(unlist(sampleReplicateGroups))))
    avgCVPerNormAndReplicates <- matrix(
        nrow=conditionLevels, 
        ncol=methodCount, 
        byrow=TRUE
    )
    
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

#' Calculate CV values for each feature. Iterates through each normalization 
#' method and calculates a matrix of CV values where each column correspond to 
#' a method and each row corresponds to a feature.
#' 
#' @param methodList List containing normalized matrices.
#' @param sampleReplicateGroups Condition header.
#' @return methodFeatureCVMatrix Matrix with feature as rows and normalization
#'   method as columns
#' @keywords internal
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

#' Calculate average MAD (Median Absolute Deviation) for each feature in
#' each condition and then calculates the average for each replicate group
#' 
#' @param methodList List containing normalized matrices.
#' @param sampleReplicateGroups Condition header.
#' @return condAvgMadMat Matrix with average MAD for each biological condition.
#' @keywords internal
calculateAvgMadMem <- function(methodList, sampleReplicateGroups) {
    
    methodCount <- length(methodList)
    conditionLevels <- length(levels(as.factor(unlist(sampleReplicateGroups))))
    condAvgMadMat <- matrix(nrow=conditionLevels, ncol=methodCount, byrow=TRUE)
    indexList <- getIndexList(sampleReplicateGroups)
    
    for (methodIndex in seq_len(methodCount)) {
        
        processedDataMatrix <- methodList[[methodIndex]]
        
        featureMedianAbsDevMat <- matrix(
            nrow=nrow(processedDataMatrix), 
            ncol=length(levels(as.factor(unlist(sampleReplicateGroups)))), 
            byrow=TRUE)
        
        for (sampleIndex in seq_len(length(names(indexList)))) {
            repVal <- names(indexList)[sampleIndex]
            cols <- indexList[[repVal]]
            featureMedianAbsDevMat[, sampleIndex] <- apply(
                processedDataMatrix[, cols], 
                1, 
                function(x) { stats::mad(x, na.rm=TRUE) })
        }
        
        condAvgMadMat[, methodIndex] <- apply(
            featureMedianAbsDevMat, 
            2, 
            mean, 
            na.rm=TRUE)
    }
    condAvgMadMat
}

#' Calculate average variance for each feature in each condition and then
#' calculate the average for each replicate group
#' 
#' @param methodList List containing normalized matrices.
#' @param sampleReplicateGroups Condition header.
#' @return avgVarianceMat Matrix with average variance for each biological
#' condition
#' @keywords internal
calculateAvgReplicateVariation <- function(methodList, sampleReplicateGroups) {

    methodCount <- length(methodList)
    indexList <- getIndexList(sampleReplicateGroups)
    conditionLevelCount <- length(unique(sampleReplicateGroups))
    avgVarianceMat <- matrix(
        nrow=conditionLevelCount, 
        ncol=methodCount, 
        byrow=TRUE)

    for (methodIndex in seq_len(methodCount)) {
        
        replicateGroupVariance <- vector()
        rowVariances <- vector()
        processedDataMatrix <- methodList[[methodIndex]]
        
        for (sampleIndex in seq_len(length(names(indexList)))) {
            
            repVal <- names(indexList)[sampleIndex]
            cols <- indexList[[repVal]]
            
            rowNonNACount <- apply(
                processedDataMatrix[, cols], 
                1, 
                function(x) { sum(!is.na(x)) }) - 1
            
            rowVariances <- rowNonNACount * apply(
                processedDataMatrix[, cols], 
                1, 
                function(x) { stats::var(x, na.rm=TRUE) })
            
            replicateGroupVariance <- c(
                replicateGroupVariance, 
                sum(rowVariances, na.rm=TRUE) / sum(rowNonNACount, na.rm=TRUE))
        }
        
        avgVarianceMat[, methodIndex] <- replicateGroupVariance
    }
    
    avgVarianceMat
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
#' @keywords internal
calculateSummarizedCorrelationVector <- function(
    methodlist, allReplicateGroups, sampleGroupsWithReplicates, corrType) {
    
    validCorrTypes <- c("pearson", "spearman")
    if (!corrType %in% validCorrTypes) {
        stop("Unknown correlation type: ", 
             corrType, 
             " valid are: ", 
             paste(validCorrTypes, collapse=", "))
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
#' @param methodData Expression data matrix
#' @param allReplicateGroups Full condition header corresponding to data tables
#'   columns
#' @param sampleGroupsWithReplicates Unique conditions where number of
#'   replicates exceeds one
#' @param corrType Type of correlation (Pearson or Spearman)
#' @return corSums
#' @keywords internal
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
            corSums <- c(
                corSums, 
                corVals[index, -(seq_len(index)), drop="FALSE"]
            )
        }
    }
    
    corSums
}

#' Calculates ANOVA p-values comparing the different condition groups
#' returning a vector of resulting p-values with NA-values where too few
#' values were present in at least one of the groups.
#' 
#' @param methodList List containing normalized matrices
#' @param sampleReplicateGroups Condition header
#' @param categoricalANOVA Whether the ANOVA should be calculated using ordered
#' or categorical groups
#' @return avgVarianceMat Matrix with average variance for each biological
#' condition
#' @keywords internal
calculateANOVAPValues <- function(methodList, 
                                  sampleReplicateGroups, 
                                  categoricalANOVA) {
    
    anovaPVals <- list()
    methodCount <- length(methodList)
    anovaPValsWithNA <- matrix(NA, ncol=methodCount, nrow=nrow(methodList[[1]]))
  
    for (methodIndex in seq_len(methodCount)) {
        
        processedDataMatrix <- methodList[[methodIndex]]
        naFilterContrast <- getRowNAFilterContrast(
            processedDataMatrix,
            sampleReplicateGroups,
            2)
        
        naFilteredData <- processedDataMatrix[naFilterContrast,]
        
        if (categoricalANOVA) {
            testLevels <- factor(sampleReplicateGroups)
        }
        else {
            testLevels <- sampleReplicateGroups
        }
        
        anovaPValCol <- apply(
            naFilteredData, 
            1, 
            function(sampleIndex) {
                summary(stats::aov(unlist(sampleIndex)~testLevels))[[1]][[5]][1]
            } 
        )
        
        anovaPValsWithNA[naFilterContrast, methodIndex] <- anovaPValCol
    }
    
    anovaPValsWithNA
}

#' Uses a list of FDR-values to extract features with low variance in the
#' log2-transformed dataset. This is then used to calculate the average CV
#' for these 'lowly variable' features in each normalization approach
#' 
#' @param referenceFDR List of FDR values used as non-normalized reference
#' @param methodList List containing normalized matrices
#' @return lowVarFeaturesAverageCVs Average CV values for lowly variable
#' features in each normalization approach
#' @keywords internal
findLowlyVariableFeatures <- function(referenceFDR, methodList) {
    
    fivePercCount <- 5 * length(referenceFDR) / 100
    pThresLowVal <- min(utils::head(rev(sort(referenceFDR)), n=fivePercCount))
    nbrAbovePThres <- sum(referenceFDR >= pThresLowVal)
    
    if (is.infinite(pThresLowVal) || nbrAbovePThres == 0) {
        warning("Too few successful ANOVA calculations to generate lowly ",
                "variable features, skipping")
        return(NULL)
    }
        
    lowlyVariableFeatures <- which(referenceFDR >= pThresLowVal)
    lowVarFeaturesAverageCVs <- vector()
    methodCount <- length(methodList)
    
    for (mlist in seq_len(methodCount)) {
        lowVarFeatures <- methodList[[mlist]][lowlyVariableFeatures, ]
        lowVarFeaturesAverageCVs[mlist] <- mean(
            apply(
                lowVarFeatures, 
                1, 
                function(sampleIndex) raster::cv(sampleIndex, na.rm=TRUE)), 
            na.rm=TRUE)
    }
    
    lowVarFeaturesAverageCVs
}
