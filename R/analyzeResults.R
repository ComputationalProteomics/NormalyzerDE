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
    
    methodList <- getNormalizationMatrices(nr)
    methodNames <- getUsedMethodNames(nr)
    rawData <- nds@rawData
    filterRawData <- nds@filterrawdata
    sampleReplicateGroups <- nds@sampleReplicateGroups
    getEDdata <- nds@inputHeaderValues
    HKflag <- !is.null(nr@houseKeepingVars)
    
    pooledVarMem <- vector()
    medianOfSamples <- vector()
    dataMadSamplesMem <- vector()
    dataLabels <- vector()
    
    methodCount <- length(methodNames)
    rowCount <- length(levels(as.factor(unlist(sampleReplicateGroups))))
    
    avgcvmem <- matrix(nrow=rowCount, ncol=methodCount, byrow=TRUE)
    avgmadmem <- matrix(nrow=rowCount, ncol=methodCount, byrow=TRUE)
    avgvarmem <- matrix(nrow=rowCount, ncol=methodCount, byrow=TRUE)
    datamadpeptmem <- matrix(nrow=nrow(filterRawData), 
                             ncol=methodCount, byrow=TRUE)
    
    anovaPVal <- vector()
    anovaFDR <- vector()
    krusWalPVal <- vector()
    krusValFDR <- vector()
    
    firstIndices <- getFirstIndicesInVector(sampleReplicateGroups, reverse=FALSE)
    lastIndices <- getFirstIndicesInVector(sampleReplicateGroups, reverse=TRUE)
    
    print(firstIndices)
    print(lastIndices)

    # Here, we iterate through the different normalizations
    for (methodIndex in 1:methodCount) {
        
        # print(paste("Processing: ", methodlist[[meti]]))
        print(paste("Processing methodlist index: ", methodIndex))
        
        # datastore <- methods::slot(nr, slotNameList[[meti]])
        datastore <- methodList[[methodIndex]]
        templabel <- methodNames[methodIndex]
        
        datamadsamples <- NA
        datamadpeptides <- NA
        pooledvar <- NA
        dataskew <- NA
        datacv <- NA
        tempcv <- NA
        tempskew <- NA
        datakurt <- NA
        pvobj1 <- 0
        pvobj2 <- 0
        
        dataLabels <- c(dataLabels, templabel)
        
        # MAD rows
        # sampleIndex <- 1
        # startColIndex <- 1
        # endColIndex <- 1
        # flag <- 1
        # flag1 <- 1
        # count <- 0
        medianAbsDevMem <- matrix(nrow=nrow(datastore), 
                                  ncol=length(levels(as.factor(unlist(sampleReplicateGroups)))), 
                                  byrow=TRUE)
        tempcv <- vector()
        replicateGroupVariance <- vector()
        rowVariances <- vector()
        rowNonNACount <- vector()
        
        print("Before internal loop")
        
        for (sampleIndex in 1:length(firstIndices)) {
            
            startColIndex <- firstIndices[sampleIndex]
            endColIndex <- lastIndices[sampleIndex]
            
            medianAbsDevMem[, sampleIndex] <- apply(datastore[, startColIndex:endColIndex], 1, function(x) { stats::mad(x, na.rm=TRUE) })
            rowNonNACount <- apply(datastore[, startColIndex:endColIndex], 1, function(x) { sum(!is.na(x)) }) - 1
            
            rowVariances <- rowNonNACount * apply(datastore[, startColIndex:endColIndex], 1, function(x) { stats::var(x, na.rm=TRUE) })

            replicateGroupVariance <- c(replicateGroupVariance, sum(rowVariances, na.rm=TRUE) / sum(rowNonNACount, na.rm=TRUE))
        }
        
        avgvarmem[, methodIndex] <- replicateGroupVariance
        temmadmatsum <- apply(medianAbsDevMem, 2, mean, na.rm=TRUE)
        avgmadmem[, methodIndex] <- temmadmatsum
        
        tempcvmat <- matrix(nrow=nrow(datastore), ncol=length(levels(as.factor(unlist(sampleReplicateGroups)))), byrow=TRUE)
        
        for (i in 1:nrow(datastore)) {
            
            tempcv <- RcmdrMisc::numSummary(datastore[i, ], statistics=c("cv"), groups=unlist(sampleReplicateGroups))
            tempcvmat[i, ] <- tempcv$table
        }
        
        temcvmatsum <- apply(tempcvmat, 2, mean, na.rm=TRUE)
        avgcvmem[, methodIndex] <- temcvmatsum * 100
        
        # ANOVA
        nbsNAperLine <- rowSums(is.na(datastore))
        
        # Retrieving lines where at least half is non-NAs
        datastoretmp <- datastore[nbsNAperLine < (ncol(datastore) / 2), ]
        dataStoreReplicateNAFiltered <- filterLinesWithEmptySamples(datastoretmp, sampleReplicateGroups)
        
        anovaPVal <- cbind(anovaPVal, apply(dataStoreReplicateNAFiltered, 1, function(sampleIndex) summary(stats::aov(unlist(sampleIndex)~sampleReplicateGroups))[[1]][[5]][1]))
        anovaFDR <- cbind(anovaFDR, stats::p.adjust(anovaPVal[, methodIndex], method="BH"))
        
        # Kruskal Wallis
        krusWalPVal <- cbind(krusWalPVal, apply(dataStoreReplicateNAFiltered, 1, function(sampleIndex) stats::kruskal.test(unlist(sampleIndex)~sampleReplicateGroups, na.action="na.exclude")[[3]][1]))
        krusValFDR <- cbind(krusValFDR, stats::p.adjust(krusWalPVal[, methodIndex], method="BH")) 
    }
    
    # make sure first method is data2log2
    avgcvmempdiff <- sapply(1:ncol(avgcvmem), function (sampleIndex) (mean(avgcvmem[, sampleIndex]) * 100) / mean(avgcvmem[, 1]))
    avgmadmempdiff <- sapply(1:ncol(avgmadmem), function (sampleIndex) (mean(avgmadmem[, sampleIndex]) * 100) / mean(avgmadmem[, 1]))
    avgvarmempdiff <- sapply(1:ncol(avgvarmem), function (sampleIndex) (mean(avgvarmem[, sampleIndex]) * 100) / mean(avgvarmem[, 1]))
    
    # finds top 5% of least DE variables in log2 data based on ANOVA
    # generates error if it doesnt find leastDE peptides
    if (sum(anovaFDR[, 1] >= min(utils::head(rev(sort(anovaFDR[, 1])), n=(5 * nrow(anovaFDR) / 100)))) > 0) {
        nonsiganfdrlist <- which(anovaFDR[, 1] >= min(utils::head(rev(sort(anovaFDR[, 1])), n=(5 * nrow(anovaFDR) / 100))))
    }
    
    nonsiganfdrlistcv <- vector()
    for (mlist in 1:methodCount) {
        tmpdata <- methodList[[mlist]][nonsiganfdrlist, ]
        nonsiganfdrlistcv[mlist] <- mean(apply(tmpdata, 1, function(sampleIndex) raster::cv(sampleIndex, na.rm=TRUE)), na.rm=TRUE)
    }
    
    nonsiganfdrlistcvpdiff <- sapply(1:length(nonsiganfdrlistcv), function(sampleIndex) (nonsiganfdrlistcv[sampleIndex] * 100) / nonsiganfdrlistcv[1])
    print("Finished analysis, preparing plots and report...")
    
    ner <- NormalizationEvaluationResults()
    
    ner@avgcvmem <- avgcvmem
    ner@avgmadmem <- avgmadmem
    ner@avgvarmem <- avgvarmem
    ner@avgcvmempdiff <- avgcvmempdiff
    ner@avgmadmempdiff <- avgmadmempdiff
    ner@avgvarmempdiff <- avgvarmempdiff
    ner@nonsiganfdrlist <- nonsiganfdrlist
    ner@nonsiganfdrlistcvpdiff <- nonsiganfdrlistcvpdiff
    ner@anfdr <- anovaFDR
    ner@kwfdr <- krusValFDR
    
    ner <- calculateCorrelations(nr, ner)
    
    nr@ner <- ner
    
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