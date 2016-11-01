
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
                abc <- paste("Number of replicates are less than 2 for the group ", repunique[i],sep="")
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

performVSNNormalization <- function(normObj) {
    
    sink(sinkfile, type="output")
    cat("\n###Warning messages for VSN-G###\n")
    sink(sinkfile, type="message")
    normObj@data2vsn <- justvsn((filterrawdata))
    
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

normMethods <- function(datafile, currentjob) {

    print("DEBUG: normMethods entered !!")

    getrawdata <- retrieveRawData(datafile)
    
    print(getwd())
    print(currentjob)
    
    jobdir <- paste(getwd(), "/", currentjob[1], sep="")

    ## Setup steps
    setupJobDir(jobdir)
    getrawdata <- getReplicateSortedData(getrawdata)
    checkrep <- parseDataForErrors(getrawdata)
    getrawdata <- preprocessData(getrawdata)
    
    ## TODO: Collect much of this in S4 class
    HKflag <- T
    getEDdata <- ((getrawdata[1,]))
    filterrawdata1 <- getrawdata
    countna <- rowSums(!is.na(filterrawdata1[, which(checkrep > 0)]))
    filterrawdata1 <- filterrawdata1[countna >= (1 * ncol(filterrawdata1[, which(checkrep>0)])), ]
    filterED <- as.numeric(getEDdata[-which(getEDdata < 1)])

    if (nrow(filterrawdata1) < 1000) {
        # won't work without replicates
        Hkvar <- normfinder(filterrawdata1, getEDdata)
    } 
    else {
        HKflag <- F 
    }
    
    filterED <- as.numeric(getEDdata[-which(getEDdata < 1)])
    filterrawdata <- setupFilterRawData(getrawdata, getEDdata, filterED)
    
    normObj <- NormalyzerObject()

    # CONVERT TO LOG2 - Total intensity normalization
    normObj <- initializeNormObject(normObj, filterrawdata)
    
    
    colsum <- colSums(filterrawdata, na.rm=T)
    medofdata <- apply(filterrawdata, 2, FUN="median", na.rm=T)  
    meanofdata <- apply(filterrawdata, 2, FUN="mean", na.rm=T) 

    if (HKflag) {
        normObj <- calculateHKdataForNormObj(normObj, getEDdata, filterED, Hkvar, filterrawdata)
    }

    normObj <- furtherNormalizationsForNormObject(normObj, filterrawdata, colsum, medofdata, meanofdata)

    methodlist <- getMethodList(HKflag, normObj)
    methodnames <- getMethodNames(HKflag)

    sinkfile <- file(paste(jobdir, "/warnings-generated.txt", sep=""), open="wt")
    sink(sinkfile, type="message", append=F)
    
    print(paste("DEBUG: Number of rows of filtered raw data"), nrow(filterrawdata))
    
    # Perform other norm. if the dataset is not small  
    if (nrow(filterrawdata) > 100) {
        
        normObj <- performVSNNormalization(normObj, filterrawdata)
        normObj <- performQuantileNormalization(normObj)
        normObj <- performSMADNormalization(normObj)

        # SMAD normalization
        # mediandata <- apply(normObj@data2log2, 2, "median", na.rm=T)
        # maddata <- apply(normObj@data2log2, 2, function(x) mad(x, na.rm=T))
        # normObj@data2mad <- t(apply(normObj@data2log2, 1, function(x) ((x - mediandata) / maddata)))
        # normObj@data2mad <- normObj@data2mad + mean(mediandata)
        
        # Global loess
        normObj@data2loess <- normalizeCyclicLoess(normObj@data2log2, method="fast")
        
        # Global RLR
        mediandata <- apply(normObj@data2log2, 1, "median", na.rm=T)
        flag1 <- 1
        
        for (j in 1:ncol(normObj@data2log2)) {
            LRfit <- rlm(as.matrix(normObj@data2log2[, j])~mediandata, na.action=na.exclude)
            Coeffs <- LRfit$coefficients
            a<-Coeffs[2]
            b<-Coeffs[1]
            if (flag1 == 1) {
                globalfittedRLR <- (normObj@data2log2[,j]-b) / a
                flag1 <- 2
            }
            else {
                globalfittedRLR <- cbind(globalfittedRLR, (normObj@data2log2[, j] - b) / a)
            }
        }
        colnames(globalfittedRLR) <- colnames(normObj@data2log2)
        
        # NORMALIZATION within REPLICATES RLR, VSN and Loess
        {
            sink(sinkfile, type="output")
            cat("\n###Warning messages for VSN-R###\n")
            sink(sinkfile, type="message")
            warn <- NULL
            x <- 1
            z <- 1
            y <- 1
            flag <- 1
            flag1 <- 1
            data2limloess1 <- NULL
            for (i in 1:length(filterED)) {
                if (x!=filterED[i] || i==length(filterED)) {
                    
                    y <- i-1
                    if (i == length(filterED)) {
                        y <- i
                    }
                    
                    if (flag == 1) {
                        
                        # median based LR normalization
                        mediandata <- apply(normObj@data2log2[, z:y], 1, "median", na.rm=T)
                        for (j in z:y) {
                            LRfit <- rlm(as.matrix(normObj@data2log2[,j]) ~mediandata,na.action=na.exclude)
                            Coeffs <- LRfit$coefficients
                            a <- Coeffs[2]
                            b <- Coeffs[1]
                            if (flag1 == 1) {
                                fittedLR <- (normObj@data2log2[, j] - b) / a
                                flag1 <- 2
                            }
                            else {
                                fittedLR <- cbind(fittedLR, (normObj@data2log2[, j] - b) / a)
                            }
                        }
                        normObj@data2vsnrep <- justvsn(as.matrix(filterrawdata[, z:y]))
                        
                        normObj@data2limloess <- normalizeCyclicLoess(normObj@data2log2[, z:y], method="fast")
                    }
                    if (flag == 2) {  
                        # median based LR normalization
                        mediandata <- apply(normObj@data2log2[, z:y], 1, "median", na.rm=T)
                        for (j in z:y) {
                            LRfit <- rlm(as.matrix(normObj@data2log2[, j]) ~mediandata, na.action=na.exclude)
                            Coeffs <- LRfit$coefficients
                            a <- Coeffs[2]
                            b <- Coeffs[1]
                            fittedLR <- cbind(fittedLR, (normObj@data2log2[, j] - b) / a)
                        }
                        normObj@data2limloess <- cbind(normObj@data2limloess, normalizeCyclicLoess(normObj@data2log2[, z:y], method="fast"))
                        normObj@data2vsnrep <- cbind(normObj@data2vsnrep, justvsn(as.matrix(filterrawdata[, z:y])))
                    }
                    z <- i
                    x <- filterED[i]
                    flag <- 2
                }
            }
        }
        
        #sink(sinkfile,type="output")
        closeAllConnections()
        colnames(fittedLR)<-colnames(normObj@data2log2)
        colnames(normObj@data2quantile)<-colnames(normObj@data2log2)
        
        if (HKflag) {
            methodlist <- list(normObj@data2log2, normObj@data2limloess, fittedLR, normObj@data2vsnrep, normObj@data2loess, globalfittedRLR, normObj@data2vsn,
                               normObj@data2GI, normObj@data2med, normObj@data2mean, normObj@data2ctrlog, normObj@data2quantile)
            methodnames <- c("Log2", "Loess-R", "RLR-R", "VSN-R", "Loess-G", "RLR-G", "VSN-G", "TI-G", "MedI-G", "AI-G", "NF-G", "Quantile")
        } 
        else {
            methodlist <- list(normObj@data2log2, normObj@data2limloess, fittedLR, normObj@data2vsnrep, normObj@data2loess, globalfittedRLR, normObj@data2vsn, normObj@data2GI, normObj@data2med, normObj@data2mean, normObj@data2quantile)
            methodnames <- c("Log2", "Loess-R", "RLR-R", "VSN-R", "Loess-G", "RLR-G", "VSN-G", "TI-G", "MedI-G", "AI-G", "Quantile")
        }
    }
    
    #Create tmp dir for the current job
    for (i in 1:length(methodlist)) {
        write.table(file=paste(jobdir, "/", methodnames[i], "-normalized.txt", sep=""), 
                    cbind(getrawdata[-(1:2), (1:(length(getEDdata) - length(filterED)))], methodlist[[i]]), sep="\t", row.names=F,col.names=getrawdata[2,],quote=F)
    }

    if (HKflag) {
        write.table(file=paste(jobdir, "/housekeeping-variables.txt", sep=""), Hkvar, sep="\t", row.names=F, col.names=getrawdata[2,], quote=F)
    }
    
    write.table(file=paste(jobdir, "/submitted_rawdata.txt", sep=""), cbind(getrawdata[-(1:2), (1:(length(getEDdata) - length(filterED)))], filterrawdata), sep="\t", row.names=F,
                col.names=getrawdata[2,], quote=F)
    methodlist <- list(methodlist, methodnames, getrawdata, filterrawdata, filterED, HKflag)
    return(methodlist)
}
