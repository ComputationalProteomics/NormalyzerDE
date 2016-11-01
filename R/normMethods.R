
normMethods<-function(datafile,currentjob) {
    
    print("DEBUG: normMethods entered")
    
    if (class(datafile) == "character") {
        getrawdata <- as.matrix((read.table(datafile,header=F,sep="\t",stringsAsFactors=F,quote="")))    
    } 
    else if (class(datafile) == "data.frame") {
        getrawdata <- as.matrix(datafile)
    } 
    else if (class(datafile) == "matrix") {
        getrawdata <- datafile  
    }
    
    jobdir <- paste(getwd(), "/", currentjob[1], sep="")
    
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
    
    #Sort the uploaded data based on replicates
    b <- NULL
    b <- as.factor(getrawdata[1,])
    l <- levels(b)
    b <- NULL
    
    for (i in 1:length(l)) {
        b<-cbind(b,getrawdata[, which(getrawdata[1,] == l[as.numeric(i)])])
    }
    
    getrawdata <- b
    
    #Parse data for errors
    
    checkrep <- getrawdata[1,]
    repunique <- unique(checkrep)
    for (i in 1:length(repunique)) {
        if(repunique[i] != 0) {
            if(length(grep(repunique[i],checkrep)) < 2) {
                abc <- paste("Number of replicates are less than 2 for the group ", repunique[i],sep="")
                class(abc) <- "try-error"
                if (inherits(abc,"try-error")) {
                    return(abc)
                }
                stop(paste("Number of replicates are less than 2 for the group ", repunique[i], sep=""))
            }
        }
    }
    
    #replace 0 with NA
    rep0 <- getrawdata[-1,]
    rep0[which(rep0==0)] <- NA
    getrawdata <- rbind(getrawdata[1,], rep0)
    
    HKflag <- T
    getEDdata <- ((getrawdata[1,]))
    filterrawdata1 <- getrawdata
    countna <- rowSums(!is.na(filterrawdata1[,which(checkrep > 0)]))
    filterrawdata1 <- filterrawdata1[countna >= (1 * ncol(filterrawdata1[, which(checkrep>0)])), ]
    filterED <- as.numeric(getEDdata[-which(getEDdata < 1)])

    if(nrow(filterrawdata1) < 1000) {
        #wont work without replicates
        Hkvar<-normfinder(filterrawdata1,getEDdata)
    } 
    else {
        HKflag=F 
    }
    
    filterED <- as.numeric(getEDdata[-which(getEDdata < 1)])
    filterrawdata <- getrawdata[, -(1:(length(getEDdata) - length(filterED)))]
    colnames(filterrawdata) <- getrawdata[2, -(1:(length(getEDdata) - length(filterED)))]
    filterrawdata <- (as.matrix((filterrawdata[-(1:2), ])))
    class(filterrawdata) <- "numeric"
    
    #CONVERT TO LOG2
    data2log2 <- log2((filterrawdata))
    #Total intensity normalization
    data2GI <- matrix(nrow=nrow(filterrawdata), ncol=ncol(filterrawdata), byrow=T)
    data2ctr <- matrix(nrow=nrow(filterrawdata), ncol=ncol(filterrawdata), byrow=T)
    data2med <- matrix(nrow=nrow(filterrawdata), ncol=ncol(filterrawdata), byrow=T) 
    data2mean <- matrix(nrow=nrow(filterrawdata), ncol=ncol(filterrawdata), byrow=T)  
    
    colsum <- colSums(filterrawdata, na.rm=T)
    medofdata <- apply(filterrawdata, 2, FUN="median", na.rm=T)  
    meanofdata <- apply(filterrawdata, 2, FUN="mean", na.rm=T) 

    if (HKflag == T) {
        Hkvartemp <- as.matrix(Hkvar[, -(1:(length(getEDdata) - length(filterED)))])
        class(Hkvartemp) <- "numeric"
        colmedianctr <- apply(Hkvartemp, 2, FUN="mean")
        
        for(i in 1:nrow(filterrawdata)) {
            data2ctr[i,] <- unlist(sapply(1:ncol(filterrawdata), function(zd) {(filterrawdata[i, zd] / colmedianctr[zd]) * (mean(colmedianctr))}))
        }
        data2ctrlog <- log2(data2ctr)
        colnames(data2ctrlog) <- colnames(data2log2)
    }
    
    avgcolsum <- median(colsum)  
    for(i in 1:nrow(filterrawdata)) {
        data2GI[i,] <- unlist(sapply(1:ncol(filterrawdata), function(zd) {(filterrawdata[i,zd] / colsum[zd]) * (avgcolsum)}))
        data2med[i,] <- unlist(sapply(1:ncol(filterrawdata), function(zd) {(filterrawdata[i,zd] / medofdata[zd]) * mean(medofdata)}))
        data2mean[i,] <- unlist(sapply(1:ncol(filterrawdata), function(zd) {(filterrawdata[i,zd] / meanofdata[zd]) * mean(meanofdata)}))
    } 
    
    data2GI <- log2(data2GI)
    data2med <- log2(data2med)
    data2mean <- log2(data2mean)
    colnames(data2GI) <- colnames(data2log2)
    colnames(data2med) <- colnames(data2log2)
    colnames(data2mean) <- colnames(data2log2)
    
    if (HKflag == T) {
        methodlist <- list(data2log2, data2GI, data2med, data2mean, data2ctrlog)
        methodnames <- c("Log2", "TI-G", "MedI-G", "AI-G", "NF-G")
    }
    else {
        methodlist <- list(data2log2, data2GI, data2med, data2mean)
        methodnames <- c("Log2", "TI-G", "MedI-G", "AI-G")
    }
    
    sinkfile <- file(paste(jobdir, "/warnings-generated.txt", sep=""), open="wt")
    sink(sinkfile, type="message", append=F)
    
    #Perform other norm. if the dataset is not small  
    if(nrow(filterrawdata) > 100) {
        #VSN NORMALIZATION
        sink(sinkfile, type="output")
        cat("\n###Warning messages for VSN-G###\n")
        sink(sinkfile, type="message")
        data2vsn <- justvsn((filterrawdata))
        
        #QUANTILE NORMALIZATION
        data2quantile <- normalize.quantiles((data2log2), copy=T)
        
        #SMAD normalization
        mediandata <- apply(data2log2, 2, "median", na.rm=T)
        maddata <- apply(data2log2, 2, function(x) mad(x, na.rm=T))
        data2mad <- t(apply(data2log2, 1, function(x) ((x - mediandata) / maddata)))
        data2mad <- data2mad + mean(mediandata)
        
        #Global loess
        data2loess <- normalizeCyclicLoess(data2log2, method="fast")
        
        #Global RLR
        mediandata <- apply(data2log2, 1, "median", na.rm=T)
        flag1 <- 1
        
        for (j in 1:ncol(data2log2)) {
            LRfit <- rlm(as.matrix(data2log2[, j]) ~mediandata, na.action=na.exclude)
            Coeffs <- LRfit$coefficients
            a<-Coeffs[2]
            b<-Coeffs[1]
            if(flag1 == 1) {
                globalfittedRLR <- (data2log2[,j]-b) / a
                flag1 <- 2
            }
            else {
                globalfittedRLR <- cbind(globalfittedRLR, (data2log2[, j] - b) / a)
            }
        }
        colnames(globalfittedRLR) <- colnames(data2log2)
        
        #NORMALIZATION within REPLICATES RLR, VSN and  Loess
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
                    if(flag == 1) {
                        #median based LR normalization
                        
                        mediandata <- apply(data2log2[, z:y], 1, "median", na.rm=T)
                        for(j in z:y) {
                            LRfit <- rlm(as.matrix(data2log2[,j]) ~mediandata,na.action=na.exclude)
                            Coeffs <- LRfit$coefficients
                            a <- Coeffs[2]
                            b <- Coeffs[1]
                            if(flag1 == 1) {
                                fittedLR <- (data2log2[, j] - b) / a
                                flag1 <- 2
                            }
                            else {
                                fittedLR <- cbind(fittedLR, (data2log2[, j] - b) / a)
                            }
                        }
                        data2vsnrep <- justvsn(as.matrix(filterrawdata[, z:y]))
                        
                        data2limloess <- normalizeCyclicLoess(data2log2[, z:y], method="fast")
                    }
                    if (flag == 2) {  
                        #median based LR normalization
                        mediandata <- apply(data2log2[, z:y], 1, "median", na.rm=T)
                        for(j in z:y) {
                            LRfit <- rlm(as.matrix(data2log2[,j]) ~mediandata, na.action=na.exclude)
                            Coeffs <- LRfit$coefficients
                            a <- Coeffs[2]
                            b <- Coeffs[1]
                            fittedLR <- cbind(fittedLR, (data2log2[, j]-b) / a)
                        }
                        data2limloess <- cbind(data2limloess, normalizeCyclicLoess(data2log2[, z:y], method="fast"))
                        data2vsnrep <- cbind(data2vsnrep, justvsn(as.matrix(filterrawdata[, z:y])))
                    }
                    z <- i
                    x <- filterED[i]
                    flag <- 2
                }
            }
        }
        
        #sink(sinkfile,type="output")
        closeAllConnections()
        colnames(fittedLR)<-colnames(data2log2)
        colnames(data2quantile)<-colnames(data2log2)
        
        if (HKflag == T) {
            methodlist <- list(data2log2, data2limloess, fittedLR, data2vsnrep, data2loess, globalfittedRLR, data2vsn, 
                               data2GI, data2med, data2mean, data2ctrlog, data2quantile)
            methodnames <- c("Log2", "Loess-R", "RLR-R", "VSN-R", "Loess-G", "RLR-G", "VSN-G", "TI-G", "MedI-G", "AI-G", "NF-G", "Quantile")
        } 
        else {
            methodlist <- list(data2log2, data2limloess, fittedLR, data2vsnrep, data2loess, globalfittedRLR, data2vsn, data2GI, data2med, data2mean, data2quantile)
            methodnames <- c("Log2", "Loess-R", "RLR-R", "VSN-R", "Loess-G", "RLR-G", "VSN-G", "TI-G", "MedI-G", "AI-G", "Quantile")
        }
    }
    #Create tmp dir for the current job
    
    for (i in 1:length(methodlist)) {
        write.table(file=paste(jobdir, "/", methodnames[i], "-normalized.txt", sep=""), 
                    cbind(getrawdata[-(1:2), (1:(length(getEDdata) - length(filterED)))], methodlist[[i]]), sep="\t", row.names=F,col.names=getrawdata[2,],quote=F)
    }

    if(HKflag == T) {
        write.table(file=paste(jobdir, "/housekeeping-variables.txt", sep=""), Hkvar, sep="\t", row.names=F, col.names=getrawdata[2,], quote=F)
    }
    
    write.table(file=paste(jobdir, "/submitted_rawdata.txt", sep=""), cbind(getrawdata[-(1:2), (1:(length(getEDdata) - length(filterED)))], filterrawdata), sep="\t", row.names=F,
                col.names=getrawdata[2,], quote=F)
    methodlist <- list(methodlist, methodnames, getrawdata, filterrawdata, filterED, HKflag)
    return(methodlist)
}
