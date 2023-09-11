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
#' data(example_summarized_experiment)
#' normObj <- getVerifiedNormalyzerObject("job_name", example_summarized_experiment)
#' normResults <- normMethods(normObj)
#' normResultsWithEval <- analyzeNormalizations(normResults)
analyzeNormalizations <- function(nr, categoricalAnova=FALSE) {

    ner <- NormalyzerEvaluationResults(nr)
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
    
    calculateFeatureCVs <- function(feature, groups) {

        featureCVs <- vapply(
            unique(groups), 
            function(group) {
                targetFeatures = feature[groups == group]
                stats::sd(targetFeatures, na.rm=TRUE) / mean(targetFeatures, na.rm=TRUE) 
            }, 
            c(1)
        )
        featureCVs
    }
    
    calculateMethodReplicateCVs <- function(methodData, groups) {
        
        cvPerFeatureAndGroup <- apply(
            methodData,
            1,
            calculateFeatureCVs,
            groups=sampleReplicateGroups
        )
        
        if (length(unique(groups)) > 1) {
            # Transpose to get groups as column heads
            featureCondCVs <- t(cvPerFeatureAndGroup)
        }
        else if (length(unique(groups)) == 1) {
            # For one feature the CVs are dropped to a vector
            featureCondCVs <- data.frame(cvPerFeatureAndGroup)
            colnames(featureCondCVs) <- groups[1]
        }
        else {
            stop("Unknown state encountered for groups:", paste(groups, collapse=", "))
        }
        
        # Calculate mean CV for all features for each condition
        summedCondCVs <- colMeans(featureCondCVs, na.rm=TRUE)
        summedCondCVs * 100
    }
    
    avgCVPerNormAndReplicates <- vapply(
        methodList,
        calculateMethodReplicateCVs,
        rep(0, length(unique(sampleReplicateGroups))),
        groups=sampleReplicateGroups
    )
    
    # If only one group, this is needed to get data in correct shape and format
    if (length(unique(sampleReplicateGroups)) == 1) {
        avgCVPerNormAndReplicates <- t(as.matrix(avgCVPerNormAndReplicates))
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

    calculateFeatureCVVector <- function(processedDataMatrix) {
        cv <- function(row) {
            stDev <- stats::sd(row, na.rm=TRUE)
            meanVal <- mean(unlist(row), na.rm=TRUE)
            stDev / mean(meanVal) * 100
        }
        featureCVs <- apply(processedDataMatrix, 1, cv)
        featureCVs
    }
    
    methodFeatureCVMatrix <- vapply(
        methodList,
        calculateFeatureCVVector,
        rep(0, numberFeatures)
    )

    colnames(methodFeatureCVMatrix) <- names(methodList)
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
    
    groupIndexList <- getIndexList(sampleReplicateGroups)

    calculateAvgFeatureMadForGroup <- function(groupIndices, methodData) {
        groupData <- methodData[, groupIndices]
        featureMAD <- matrixStats::rowMads(groupData, na.rm=TRUE)
        featureMAD
    }
        
    calculateGroupMadForMethod <- function(methodData, groups, indexList) {

        # Extracts groups of replicates and calculate MAD for each feature
        featureMADMat <- vapply(
            indexList,
            calculateAvgFeatureMadForGroup,
            rep(0, nrow(methodData)),
            methodData=methodData
        )
        
        methodRepGroupMADMean <- colMeans(featureMADMat, na.rm=TRUE)
        methodRepGroupMADMean
    }
    
    condAvgMadMat <- vapply(
        methodList,
        calculateGroupMadForMethod,
        rep(0, length(unique(sampleReplicateGroups))),
        groups=sampleReplicateGroups,
        indexList=groupIndexList
    ) 

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

    groupIndexList <- getIndexList(sampleReplicateGroups)
    
    calculateReplicateGroupVariance <- function(groupIndices, methodData) {
        
        groupData <- methodData[, groupIndices]
        rowNonNACount <- rowSums(!is.na(groupData)) - 1
        rowVariances <- rowNonNACount * matrixStats::rowVars(groupData, na.rm=TRUE)
        replicateGroupVariance <- sum(rowVariances, na.rm=TRUE) / sum(rowNonNACount, na.rm=TRUE)
        replicateGroupVariance
    }
    
    avgVarianceMat <- vapply(
        methodList,
        function(methodData) {
            replicateGroupVariance <- vapply(
                groupIndexList,
                calculateReplicateGroupVariance,
                0,
                methodData=methodData
            )
            replicateGroupVariance
        },
        rep(0, length(unique(sampleReplicateGroups)))
    )
    
    avgVarianceMat
}

#' Calculates correlation values between replicates for each condition matrix.
#' Finally returns a matrix containing the results for all dataset
#' 
#' @param methodlist List containing normalized matrices for each normalization
#'   method
#' @param allReplicateGroups Vector with condition groups matching the columns
#'   found in the normalization methods
#' @param sampleGroupsWithReplicates Unique vector with condition groups
#'   present in two or more samples
#' @param corrType Type of correlation (Pearson or Spearman)
#' @return avgCorSum Matrix with column corresponding to normalization
#' approaches and rows corresponding to replicate group
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
    
    corr_combination_count <- function(allReplicateGroups) {
        replicate_counts <- table(allReplicateGroups)
        sum(vapply(
            replicate_counts, 
            function(count) { (count * (count-1)) / 2 }, 
            0))
    }
    
    avgCorSum <- vapply(
        methodlist,
        calculateCorrSum,
        rep(0, corr_combination_count(allReplicateGroups)),
        allReplicateGroups=allReplicateGroups,
        sampleGroupsWithReplicates=sampleGroupsWithReplicates,
        corrType=corrType
    )
    
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
    for (groupNbr in seq_along(sampleGroupsWithReplicates)) {

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
    
    calculateANOVAsForMethod <- function(methodData, sampleReplicateGroups) {

        naFilterContrast <- getRowNAFilterContrast(
            methodData,
            sampleReplicateGroups,
            2)
        
        naFilteredData <- methodData[naFilterContrast,, drop=FALSE]
        
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

        anovaPValsWithNACol <- rep(NA, nrow(methodData))        
        anovaPValsWithNACol[naFilterContrast] <- anovaPValCol
        anovaPValsWithNACol
    }
    
    anovaPValsWithNA <- vapply(
        methodList,
        calculateANOVAsForMethod,
        rep(0, nrow(methodList[[1]])),
        sampleReplicateGroups=sampleReplicateGroups
    )
    
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
findLowlyVariableFeaturesCVs <- function(referenceFDR, methodList) {
    
    referenceFDRWoNA <- stats::na.omit(referenceFDR)
    
    fivePercCount <- 5 * length(referenceFDRWoNA) / 100
    pThresLowVal <- min(utils::head(rev(sort(referenceFDRWoNA)), n=fivePercCount))
    nbrAbovePThres <- sum(referenceFDRWoNA >= pThresLowVal)
    
    if (is.infinite(pThresLowVal) || nbrAbovePThres == 0) {
        warning("Too few successful ANOVA calculations to generate lowly ",
                "variable features, skipping")
        return(NULL)
    }
    
    lowlyVariableFeatures <- referenceFDR >= pThresLowVal
    
    if (length(lowlyVariableFeatures) != nrow(methodList[[1]])) {
        stop("Lowly variable features contrast needs to be same length as ",
             "number of rows in matrix")
    }
    
    lowVarFeaturesAverageCVs <- vector()
    methodCount <- length(methodList)

    calculateAverageCVs <- function(methodData, lowVarContrast) {
        
        lowVarFeatures <- methodData[lowVarContrast, ]
        averageCVs <- mean(
            matrixStats::rowSds(lowVarFeatures, na.rm=TRUE) / rowMeans(lowVarFeatures, na.rm=TRUE) * 100,
            na.rm=TRUE)
        averageCVs
    }
    
    lowVarFeaturesAverageCVs <- vapply(
        methodList,
        calculateAverageCVs,
        0,
        lowVarContrast=lowlyVariableFeatures
    )
    
    lowVarFeaturesAverageCVs
}






