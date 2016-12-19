#' Calculate measures for normalization results
#' 
#' @param nr Normalyzer results object with calculated results.
#' @param name Name of the ongoing processing run.
#' @return Normalyzer results with attached evaluation results object.
#' @export
#' @examples
#' normMethods(normalyzer_results_object, "my_run")
analyzeNormalizations <- function(nr, name) {
    
    print("DEBUG: analyzeAndPlot entered")
    
    # Analysis takes too long, rewriting needed
    currentjob <- name
    
    # What is this?
    if (is.na(currentjob[2])) {
        currentjob[2] <- currentjob[1]
    }
    
    nds <- nr@nds
    
    methodlist <- getNormalizationMatrices(nr)
    methodnames <- getMethodNames(nr)
    getrawdata <- nds@rawData
    filterrawdata <- nds@filterrawdata
    filterED <- nds@sampleReplicateGroups
    getEDdata <- nds@inputHeaderValues
    HKflag <- !is.null(nr@houseKeepingVars)
    
    pooledvarmem <- vector()
    Medianofsamples <- vector()
    datamadsamplesmem <- vector()
    datalabels <- vector()
    
    methodCount <- length(methodnames)
    rowCount <- length(levels(as.factor(unlist(filterED))))
    
    avgcvmem <- matrix(nrow=rowCount, ncol=methodCount, byrow=T)
    avgmadmem <- matrix(nrow=rowCount, ncol=methodCount, byrow=T)
    avgvarmem <- matrix(nrow=rowCount, ncol=methodCount, byrow=T)
    datamadpeptmem <- matrix(nrow=nrow(filterrawdata), ncol=methodCount, byrow=T)
    
    anpvalue <- vector()
    anfdr <- vector()
    kwpvalue <- vector()
    kwfdr <- vector()
    
    for (meti in 1:methodCount) {
        
        # print(paste("Processing: ", methodlist[[meti]]))
        print(paste("Processing methodlist index: ", meti))
        
        # datastore <- slot(nr, slotNameList[[meti]])
        datastore <- methodlist[[meti]]
        templabel <- methodnames[meti]
        
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
        
        datalabels <- c(datalabels,templabel)
        #MAD rows
        x <- 1
        z <- 1
        y <- 1
        flag <- 1
        flag1 <- 1
        count <- 0
        madmem <- matrix(nrow=nrow(datastore), ncol=length(levels(as.factor(unlist(filterED)))), byrow=T)
        tempcv <- vector()
        varmem <- vector()
        tempvar <- vector()
        nonmissingmat <- vector()
        
        for (i in 1:length(filterED)) {
            if (x != filterED[i] || i == length(filterED)) {
                
                y <- i - 1
                if (i == length(filterED)) {
                    y <- i
                }
                
                if (flag == 1) {
                    count <- count + 1
                    madmem[, count] <- apply(datastore[, z:y], 1, function(x) { mad(x, na.rm=T) })
                    nonmissingmat <- (apply(datastore[, z:y], 1, function(x) { ((sum(!is.na(x)))) })) - 1
                    tempvar <- nonmissingmat * apply(datastore[, z:y], 1, function(x) { var(x, na.rm=TRUE) })
                }
                
                if (flag == 2) {  
                    count <- count + 1
                    madmem[, count] <- apply(datastore[, z:y], 1, function(x) { mad(x, na.rm=T) })  
                    nonmissingmat <- (apply(datastore[, z:y], 1, function(x) { ((sum(!is.na(x)))) })) - 1
                    tempvar <- nonmissingmat * apply(datastore[, z:y], 1, function(x) { var(x, na.rm=TRUE) })
                }
                
                varmem <- c(varmem, ((sum(tempvar, na.rm=T)) / (sum(nonmissingmat, na.rm=T))))
                z <- i
                x <- filterED[i]
                flag <- 2
            }
        }
        
        avgvarmem[, meti] <- varmem
        temmadmatsum <- apply(madmem, 2, mean, na.rm=T)
        avgmadmem[, meti] <- temmadmatsum
        
        tempcvmat<-matrix(nrow=nrow(datastore), ncol=length(levels(as.factor(unlist(filterED)))), byrow=T)
        
        for (i in 1:nrow(datastore)) {
            
            tempcv <- numSummary(datastore[i, ], statistics=c("cv"), groups=unlist(filterED))
            tempcvmat[i, ] <- tempcv$table
        }
        
        temcvmatsum <- apply(tempcvmat, 2, mean, na.rm=T)
        avgcvmem[, meti] <- ((temcvmatsum * 100))
        
        # ANOVA
        nbsNAperLine <- rowSums(is.na(datastore))
        
        # Retrieving lines where at least half is non-NAs
        datastoretmp <- datastore[nbsNAperLine < (ncol(datastore) / 2), ]
        dataStoreReplicateNAFiltered <- filterLinesWithEmptySamples(datastoretmp, filterED)
        
        anpvalue <- cbind(anpvalue, apply(dataStoreReplicateNAFiltered, 1, function(x) summary(aov(unlist(x)~filterED))[[1]][[5]][1]))
        anfdr <- cbind(anfdr, p.adjust(anpvalue[, meti], method="BH"))
        
        # Kruskal Wallis
        kwpvalue <- cbind(kwpvalue, apply(dataStoreReplicateNAFiltered, 1, function(x) kruskal.test(unlist(x)~filterED, na.action="na.exclude")[[3]][1]))
        kwfdr <- cbind(kwfdr, p.adjust(kwpvalue[, meti], method="BH")) 
    }
    
    # make sure first method is data2log2
    avgcvmempdiff <- sapply(1:ncol(avgcvmem), function (x) (mean(avgcvmem[, x]) * 100) / mean(avgcvmem[, 1]))
    avgmadmempdiff <- sapply(1:ncol(avgmadmem), function (x) (mean(avgmadmem[, x]) * 100) / mean(avgmadmem[, 1]))
    avgvarmempdiff <- sapply(1:ncol(avgvarmem), function (x) (mean(avgvarmem[, x]) * 100) / mean(avgvarmem[, 1]))
    
    # finds top 5% of least DE variables in log2 data based on ANOVA
    # generates error if it doesnt find leastDE peptides
    if (sum(anfdr[, 1] >= min(head(rev(sort(anfdr[, 1])), n=(5 * nrow(anfdr) / 100)))) > 0) {
        nonsiganfdrlist <- which(anfdr[, 1] >= min(head(rev(sort(anfdr[, 1])), n=(5 * nrow(anfdr) / 100))))
    }
    
    nonsiganfdrlistcv <- vector()
    for (mlist in 1:methodCount) {
        tmpdata <- methodlist[[mlist]][nonsiganfdrlist, ]
        nonsiganfdrlistcv[mlist] <- mean(apply(tmpdata, 1, function(x) cv(x, na.rm=T)), na.rm=T)
    }
    
    nonsiganfdrlistcvpdiff <- sapply(1:length(nonsiganfdrlistcv), function(x) (nonsiganfdrlistcv[x] * 100) / nonsiganfdrlistcv[1])
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
    ner@anfdr <- anfdr
    ner@kwfdr <- kwfdr
    
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
    
    par(mfrow=c(1,1))
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
            percor <- cor(dt, use="pairwise.complete.obs", method="pearson")
            spercor <- cor(dt, use="pairwise.complete.obs", method="spearman")
            
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
    lastIndices <- getFirstIndicesInVector(replicateHeader, reverse=T)
    
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