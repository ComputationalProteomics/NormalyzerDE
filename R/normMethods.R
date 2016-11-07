
retrieveRawData <- function(datafile) {
    if (class(datafile) == "character") {
        getrawdata <- as.matrix((read.table(datafile, header=F, sep="\t", stringsAsFactors=F, quote="")))    
    } 
    else if (class(datafile) == "data.frame") {
        getrawdata <- as.matrix(datafile)
    } 
    else if (class(datafile) == "matrix") {
        getrawdata <- datafile  
    }
    getrawdata
}

setupJobDir <- function(jobdir) {

    if (file.exists(jobdir)) {
        abc <- "Directory already exists"
        class(abc) <- "try-error"
        if (inherits(abc, "try-error")) {
            return(abc)
        }
        stop("Directory already exists")
    } 
    else {
        dir.create(jobdir)
    }
}

getReplicateSortedData <- function(getrawdata) {
    b <- NULL
    b <- as.factor(getrawdata[1,])
    l <- levels(b)
    b <- NULL
    for (i in 1:length(l)) {
        b <- cbind(b, getrawdata[, which(getrawdata[1,] == l[as.numeric(i)])])
    }
    getrawdata <- b
    getrawdata
}

## Evaluate whether the input format is in correct format
## Abort processing if this isn't the case
parseDataForErrors <- function(getrawdata) {
    
    checkrep <- getrawdata[1,]
    repunique <- unique(checkrep)
    
    for (i in 1:length(repunique)) {
        if (repunique[i] != 0) {
            if (length(grep(repunique[i],checkrep)) < 2) {
                abc <- paste("Number of replicates are less than 2 for the group ", repunique[i], sep="")
                class(abc) <- "try-error"
                if (inherits(abc,"try-error")) {
                    return(abc)
                }
                stop(paste("Number of replicates are less than 2 for the group ", repunique[i], sep=""))
            }
        }
    }
    checkrep
}

## Replaces zeroes with NAs
preprocessData <- function(getrawdata) {
    
    rep0 <- getrawdata[-1,]
    rep0[which(rep0==0)] <- NA
    getrawdata <- rbind(getrawdata[1,], rep0)
    
    getrawdata
}

setupFilterRawData <- function(getrawdata, getEDdata, filterED) {
    
    filterrawdata <- getrawdata[, -(1:(length(getEDdata) - length(filterED)))]
    colnames(filterrawdata) <- getrawdata[2, -(1:(length(getEDdata) - length(filterED)))]
    filterrawdata <- (as.matrix((filterrawdata[-(1:2), ])))
    class(filterrawdata) <- "numeric"
    
    filterrawdata
}

initializeNormObject <- function(normObj, filterrawdata) {
    
    normObj@data2log2 <- log2(filterrawdata)
    normObj@data2GI <- matrix(nrow=nrow(filterrawdata), ncol=ncol(filterrawdata), byrow=T)
    normObj@data2ctr <- matrix(nrow=nrow(filterrawdata), ncol=ncol(filterrawdata), byrow=T)
    normObj@data2med <- matrix(nrow=nrow(filterrawdata), ncol=ncol(filterrawdata), byrow=T) 
    normObj@data2mean <- matrix(nrow=nrow(filterrawdata), ncol=ncol(filterrawdata), byrow=T)
    
    normObj
}

calculateHKdataForNormObj <- function(normObj, getEDdata, filterED, Hkvar, filterrawdata) {
    
    Hkvartemp <- as.matrix(Hkvar[, -(1:(length(getEDdata) - length(filterED)))])
    class(Hkvartemp) <- "numeric"
    colmedianctr <- apply(Hkvartemp, 2, FUN="mean")
    
    for(i in 1:nrow(filterrawdata)) {
        normObj@data2ctr[i,] <- unlist(sapply(1:ncol(filterrawdata), 
                                              function(zd) {(filterrawdata[i, zd] / colmedianctr[zd]) * (mean(colmedianctr))}))
    }
    normObj@data2ctrlog <- log2(normObj@data2ctr)
    colnames(normObj@data2ctrlog) <- colnames(normObj@data2log2)
    
    normObj
}

## TODO: Understand this better. Really.
furtherNormalizationsForNormObject <- function(normObj, filterrawdata, colsum, medofdata, meanofdata) {
    
    avgcolsum <- median(colsum)
    for(i in 1:nrow(filterrawdata)) {
        normObj@data2GI[i,] <- unlist(sapply(1:ncol(filterrawdata), function(zd) {(filterrawdata[i,zd] / colsum[zd]) * (avgcolsum)}))
        normObj@data2med[i,] <- unlist(sapply(1:ncol(filterrawdata), function(zd) {(filterrawdata[i,zd] / medofdata[zd]) * mean(medofdata)}))
        normObj@data2mean[i,] <- unlist(sapply(1:ncol(filterrawdata), function(zd) {(filterrawdata[i,zd] / meanofdata[zd]) * mean(meanofdata)}))
    } 
    
    normObj@data2GI <- log2(normObj@data2GI)
    normObj@data2med <- log2(normObj@data2med)
    normObj@data2mean <- log2(normObj@data2mean)
    colnames(normObj@data2GI) <- colnames(normObj@data2log2)
    colnames(normObj@data2med) <- colnames(normObj@data2log2)
    colnames(normObj@data2mean) <- colnames(normObj@data2log2)
    
    normObj
}

getMethodList <- function(HKflag, normObj) {
    
    if (HKflag) {
        methodlist <- list(normObj@data2log2, normObj@data2GI, normObj@data2med, normObj@data2mean, normObj@data2ctrlog)
    }
    else {
        methodlist <- list(normObj@data2log2, normObj@data2GI, normObj@data2med, normObj@data2mean)
    }
    methodlist
}

getMethodNames <- function(HKflag) {
    
    if (HKflag) {
        methodnames <- c("Log2", "TI-G", "MedI-G", "AI-G", "NF-G")
    }
    else {
        methodnames <- c("Log2", "TI-G", "MedI-G", "AI-G")
    }
    methodnames
}

performVSNNormalization <- function(normObj, filterrawdata) {
    
    # sink(sinkfile, type="output")
    cat("\n###Warning messages for VSN-G###\n")
    # sink(sinkfile, type="message")
    normObj@data2vsn <- justvsn(filterrawdata)
    
    normObj
}

performQuantileNormalization <- function(normObj) {
    normObj@data2quantile <- normalize.quantiles((normObj@data2log2), copy=T)
    normObj
}

performSMADNormalization <- function(normObj) {
    mediandata <- apply(normObj@data2log2, 2, "median", na.rm=T)
    maddata <- apply(normObj@data2log2, 2, function(x) mad(x, na.rm=T))
    normObj@data2mad <- t(apply(normObj@data2log2, 1, function(x) ((x - mediandata) / maddata)))
    normObj@data2mad <- normObj@data2mad + mean(mediandata)
    normObj
}

performCyclicLoessNormalization <- function(normObj) {
    normObj@data2loess <- normalizeCyclicLoess(normObj@data2log2, method="fast")
    normObj
}

performGlobalRLRNormalization <- function(normObj) {
    
    mediandata <- apply(normObj@data2log2, 1, "median", na.rm=T)
    isFirstSample <- TRUE

    for (j in 1:ncol(normObj@data2log2)) {

        LRfit <- rlm(as.matrix(normObj@data2log2[, j])~mediandata, na.action=na.exclude)
        Coeffs <- LRfit$coefficients
        coeffs2 <- Coeffs[2]
        coeffs1 <- Coeffs[1]

        if (isFirstSample) {
            value_to_assign <- (normObj@data2log2[, j] - coeffs1) / coeffs2
            globalfittedRLR <- (normObj@data2log2[, j] - coeffs1) / coeffs2
            isFirstSample <- FALSE
        }
        else {
            globalfittedRLR <- cbind(globalfittedRLR, (normObj@data2log2[, j] - coeffs1) / coeffs2)
        }
    }

    normObj@globalfittedRLR <- globalfittedRLR
    colnames(normObj@globalfittedRLR) <- colnames(normObj@data2log2)

    normObj
}


## Retrieve index vector with first or last occurences of each element
getFirstIndicesInVector <- function(targetVector, reverse=F) {
    
    encounteredNumbers <- c()
    firstIndices <- c()
    
    if (!reverse) {
        startIndex <- 1
        endIndex <- length(targetVector)
    }
    else {
        startIndex <- length(targetVector)
        endIndex <- 1
    }
    
    for (i in startIndex:endIndex) {
        targetValue <- targetVector[i]
        
        if (!is.element(targetValue, encounteredNumbers)) {
            
            encounteredNumbers <- append(encounteredNumbers, targetValue)
            
            if (!reverse) {
                firstIndices <- append(firstIndices, i)
            }
            else {
                firstIndices <- append(i, firstIndices)
            }
        }
    }
    
    firstIndices
}

normMethods <- function(datafile, currentjob) {

    IGNORE_NORMFINDER = F  # Leads to downstreams crash for now
    
    print("DEBUG: normMethods entered !!")

    getrawdata <- retrieveRawData(datafile)
    
    print(getwd())
    print(currentjob)
    
    jobdir <- paste(getwd(), "/", currentjob[1], sep="")

    ## Setup steps
    print("DEBUG: normMethods setup")
    setupJobDir(jobdir)
    getrawdata <- getReplicateSortedData(getrawdata)
    checkrep <- parseDataForErrors(getrawdata)
    getrawdata <- preprocessData(getrawdata)
    
    ## TODO: Collect much of this in S4 class
    print("DEBUG: normMethods house keeping")
    getEDdata <- ((getrawdata[1,]))
    
    # TODO: Get rid of filterrawdata1, meaning compared to filterrawdata?
    filterrawdata1 <- getrawdata
    countna <- rowSums(!is.na(filterrawdata1[, which(checkrep > 0)]))
    filterrawdata1 <- filterrawdata1[countna >= (1 * ncol(filterrawdata1[, which(checkrep>0)])), ]
    sampleReplicateGroups <- as.numeric(getEDdata[-which(getEDdata < 1)])

    if (nrow(filterrawdata1) < 1000 && !IGNORE_NORMFINDER) {
        # won't work without replicates

        Hkvar <- normfinder(filterrawdata1, getEDdata)
        HKflag <- TRUE
    } 
    else {
        HKflag <- FALSE
    }

    print("DEBUG: normMethods after house keeping")
        
    sampleReplicateGroups <- as.numeric(getEDdata[-which(getEDdata < 1)])
    filterrawdata <- setupFilterRawData(getrawdata, getEDdata, sampleReplicateGroups)

    normObj <- NormalyzerObject()

    # CONVERT TO LOG2 - Total intensity normalization
    normObj <- initializeNormObject(normObj, filterrawdata)
    
    colsum <- colSums(filterrawdata, na.rm=T)
    medofdata <- apply(filterrawdata, 2, FUN="median", na.rm=T)  
    meanofdata <- apply(filterrawdata, 2, FUN="mean", na.rm=T) 

    if (HKflag) {
        normObj <- calculateHKdataForNormObj(normObj, getEDdata, sampleReplicateGroups, Hkvar, filterrawdata)
    }

    normObj <- furtherNormalizationsForNormObject(normObj, filterrawdata, colsum, medofdata, meanofdata)

    methodlist <- getMethodList(HKflag, normObj)
    methodnames <- getMethodNames(HKflag)

    print("DEBUG: normMethods before sinkfile")
    
    # sinkfile <- file(paste(jobdir, "/warnings-generated.txt", sep=""), open="wt")
    # sink(sinkfile, type="message", append=F)
    
    print("DEBUG: normMethods after pre-existing sinkfile")

    
    # Perform other norm. if the dataset is not small  
    if (nrow(filterrawdata) > 50) {

        print("DEBUG: vsn")        
        normObj <- performVSNNormalization(normObj, filterrawdata)
        print("DEBUG: Quantile")
        normObj <- performQuantileNormalization(normObj)
        print("DEBUG: SMAD")
        normObj <- performSMADNormalization(normObj)
        print("DEBUG: Cyclic Loess")
        normObj <- performCyclicLoessNormalization(normObj)
        print("DEBUG: Global RLR")
        normObj <- performGlobalRLRNormalization(normObj)

        # NORMALIZATION within REPLICATES RLR, VSN and Loess
        
        # Jakob translation: We are normalyzing within samples belonging to same group delimitation,
        # set by the header replicate indices
        
        # sink(sinkfile, type="output")
        cat("\n### Warning messages for VSN-R ###\n")
        # sink(sinkfile, type="warning")
        warn <- NULL
        # currentReplicateGroup <- 1
        startColIndex <- 1  # Used to slice columns, start column
        endColIndex <- 1    # Used to slice columns, end column
        isFirstFittedLR <- TRUE
        data2limloess1 <- NULL
        
        print("DEBUG: Before loop")
        print("DEBUG: sampleReplicateGroups")
        print(sampleReplicateGroups)
        
        firstIndices <- getFirstIndicesInVector(sampleReplicateGroups, reverse=F)
        lastIndices <- getFirstIndicesInVector(sampleReplicateGroups, reverse=T)
        stopifnot(length(firstIndices) == length(lastIndices))

        normObj@fittedLR <- matrix(, nrow=nrow(filterrawdata), ncol=0)
        normObj@data2limloess <- matrix(, nrow=nrow(filterrawdata), ncol=0)
        normObj@data2vsnrep <- matrix(, nrow=nrow(filterrawdata), ncol=0)
        
        for (sampleIndex in 1:length(firstIndices)) {
            
            print(paste("DEBUG: Current sampleIndex:", sampleIndex))
            print(paste("DEBUG: Inner IF statement executed for index: ", sampleIndex))

            startColIndex <- firstIndices[sampleIndex]
            endColIndex <- lastIndices[sampleIndex]

            # Median based LR normalization
            mediandata <- apply(normObj@data2log2[, startColIndex:endColIndex], 1, "median", na.rm=T)

            for (currentCol in startColIndex:endColIndex) {
                
                LRfit <- rlm(as.matrix(normObj@data2log2[, currentCol])~mediandata, na.action=na.exclude)
                Coeffs <- LRfit$coefficients
                LRcoeff2 <- Coeffs[2]
                LRcoeff1 <- Coeffs[1]
                normObj@fittedLR <- cbind(normObj@fittedLR, (normObj@data2log2[, currentCol] - LRcoeff1) / LRcoeff2)
            }
            
            normObj@data2limloess <- cbind(normObj@data2limloess, normalizeCyclicLoess(normObj@data2log2[, startColIndex:endColIndex], method="fast"))
            normObj@data2vsnrep <- cbind(normObj@data2vsnrep, justvsn(as.matrix(filterrawdata[, startColIndex:endColIndex])))
        }

        print("DEBUG: After loop")
        
        print(sprintf("DEBUG: filterrawdata lines: %i", nrow(filterrawdata)))
        print(sprintf("DEBUG: normObj@data2log2 lines: %i", nrow(normObj@data2log2)))
        
        print(sprintf("DEBUG: normObj@data2limloess lines: %i", nrow(normObj@data2limloess)))
        print(sprintf("DEBUG: normObj@fittedLR lines: %i", nrow(normObj@fittedLR)))
        print(sprintf("DEBUG: normObj@data2vsnrep lines: %i", nrow(normObj@data2vsnrep)))

        print("DEBUG: ------- AFTER SPRINTF STATEMENTS")
                
        #sink(sinkfile,type="output")
        closeAllConnections()
        colnames(normObj@fittedLR) <- colnames(normObj@data2log2)
        colnames(normObj@data2quantile) <- colnames(normObj@data2log2)
        
        if (HKflag) {
            methodlist <- list(normObj@data2log2, normObj@data2limloess, normObj@fittedLR, normObj@data2vsnrep, normObj@data2loess, 
                               normObj@globalfittedRLR, normObj@data2vsn, normObj@data2GI, normObj@data2med, normObj@data2mean, 
                               normObj@data2ctrlog, normObj@data2quantile)
            methodnames <- c("Log2", "Loess-R", "RLR-R", "VSN-R", "Loess-G", "RLR-G", "VSN-G", "TI-G", "MedI-G", "AI-G", "NF-G", "Quantile")
        } 
        else {
            methodlist <- list(normObj@data2log2, normObj@data2limloess, normObj@fittedLR, normObj@data2vsnrep, normObj@data2loess, 
                               normObj@globalfittedRLR, normObj@data2vsn, normObj@data2GI, normObj@data2med, normObj@data2mean, normObj@data2quantile)
            methodnames <- c("Log2", "Loess-R", "RLR-R", "VSN-R", "Loess-G", "RLR-G", "VSN-G", "TI-G", "MedI-G", "AI-G", "Quantile")
        }
        
        print("METHODLISTLENGTH: ")
        print(length(methodlist))
    }
    
    #Create tmp dir for the current job
    for (sampleIndex in 1:length(methodlist)) {
        write.table(file=paste(jobdir, "/", methodnames[sampleIndex], "-normalized.txt", sep=""), 
                    cbind(getrawdata[-(1:2), (1:(length(getEDdata) - length(sampleReplicateGroups)))], methodlist[[sampleIndex]]), sep="\t", row.names=F, col.names=getrawdata[2,], quote=F)
    }

    if (HKflag) {
        write.table(file=paste(jobdir, "/housekeeping-variables.txt", sep=""), Hkvar, sep="\t", row.names=F, col.names=getrawdata[2,], quote=F)
    }
    
    write.table(file=paste(jobdir, "/submitted_rawdata.txt", sep=""), 
                cbind(getrawdata[-(1:2), (1:(length(getEDdata) - length(sampleReplicateGroups)))], filterrawdata), sep="\t", row.names=F,
                col.names=getrawdata[2,], quote=F)
    methodlist <- list(methodlist, methodnames, getrawdata, filterrawdata, sampleReplicateGroups, HKflag)
    
    # sink(NULL)
    
    return(methodlist)
}
