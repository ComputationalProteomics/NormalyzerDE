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

#' Calculate CV per replicate group and normalization technique
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

#' Calculate CV values for each feature
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

#' General function for calculating percentage difference of average column
#' means in matrix
#' 
#' @param targetMat Matrix for which column means should be compared
#' @return percDiffVector Vector with percentage difference, where first element
#'   always will be 100
#' @keywords internal
calculatePercentageAvgDiffInMat <- function(targetMat) {

    calculatePercDiff <- function (sampleIndex, mat) {
        mean(mat[, sampleIndex]) * 100 / mean(mat[, 1])
    }
    
    percDiffVector <- vapply(
        seq_len(ncol(targetMat)), 
        calculatePercDiff,
        0,
        mat=targetMat)
    
    percDiffVector
}


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
#' @export
#' @keywords internal
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
#' @param methodData Expression data matrix
#' @param allReplicateGroups Full condition header corresponding to data tables
#'   columns
#' @param sampleGroupsWithReplicates Unique conditions where number of
#'   replicates exceeds one
#' @param corrType Type of correlation (Pearson or Spearman)
#' @return corSums
#' @export
#' @keywords internal
calculateCorrSum <- function(
    methodData, allReplicateGroups, sampleGroupsWithReplicates, corrType) {
    
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

calculateANOVAPValues <- function(methodList, sampleReplicateGroups, categoricalANOVA) {
    
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
            function(sampleIndex) 
                summary(stats::aov(unlist(sampleIndex)~testLevels))[[1]][[5]][1])
        
        anovaPValsWithNA[naFilterContrast, methodIndex] <- anovaPValCol
    }
    
    anovaPValsWithNA
}

findLowlyVariableFeatures <- function(log2AnovaFDR, methodList) {
    
    fivePercCount <- 5 * length(log2AnovaFDR) / 100
    pThresLowVal <- min(utils::head(rev(sort(log2AnovaFDR)), n=fivePercCount))
    nbrAbovePThres <- sum(log2AnovaFDR >= pThresLowVal)
    
    if (is.infinite(pThresLowVal) || nbrAbovePThres == 0) {
        warning(paste("Too few successful ANOVA calculations to generate lowly",
                      "variable features, skipping"))
        return(NULL)
    }
        
    lowlyVariableFeatures <- which(log2AnovaFDR >= pThresLowVal)
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
