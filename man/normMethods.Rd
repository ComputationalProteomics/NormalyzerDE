\name{normMethods}
\alias{normMethods}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
normMethods
}
\description{
Normalizes data with different methods. This method is NOT built for direct execution.
}
\usage{
normMethods(file, currentjob)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{
%%     ~~Describe \code{file} here~~
}
  \item{currentjob}{
%%     ~~Describe \code{currentjob} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (file, currentjob) 
{
    currentjob = "name"
    file = "submitted_rawdata.txt"
    getrawdata <- as.matrix((read.table(file, header = F, sep = "\t", 
        stringsAsFactors = F)))
    b <- NULL
    b <- as.factor(getrawdata[1, ])
    l <- levels(b)
    b <- NULL
    for (i in 1:length(l)) {
        b <- cbind(b, getrawdata[, which(getrawdata[1, ] == l[as.numeric(i)])])
    }
    getrawdata <- b
    HKflag = T
    getEDdata <- ((getrawdata[1, ]))
    filterrawdata1 <- getrawdata
    countna <- rowSums(!is.na(filterrawdata1))
    filterrawdata1 <- filterrawdata1[countna >= (1 * ncol(filterrawdata1)), 
        ]
    filterED <- as.numeric(getEDdata[-which(getEDdata < 1)])
    if (nrow(filterrawdata1) < 1000) {
        setwd("~/normomics")
        source("normfinder-pipeline.R")
        Hkvar <- normfinder(filterrawdata1, getEDdata)
    }
    else {
        HKflag = F
    }
    filterED <- as.numeric(getEDdata[-which(getEDdata < 1)])
    filterrawdata <- getrawdata[, -(1:(length(getEDdata) - length(filterED)))]
    colnames(filterrawdata) <- getrawdata[2, -(1:(length(getEDdata) - 
        length(filterED)))]
    filterrawdata <- (as.matrix((filterrawdata[-(1:2), ])))
    class(filterrawdata) <- "numeric"
    data2log2 <- log2((filterrawdata))
    data2GI <- matrix(nrow = nrow(filterrawdata), ncol = ncol(filterrawdata), 
        byrow = T)
    data2ctr <- matrix(nrow = nrow(filterrawdata), ncol = ncol(filterrawdata), 
        byrow = T)
    data2med <- matrix(nrow = nrow(filterrawdata), ncol = ncol(filterrawdata), 
        byrow = T)
    data2mean <- matrix(nrow = nrow(filterrawdata), ncol = ncol(filterrawdata), 
        byrow = T)
    colsum <- colSums(filterrawdata, na.rm = T)
    medofdata <- apply(filterrawdata, 2, FUN = "median", na.rm = T)
    meanofdata <- apply(filterrawdata, 2, FUN = "mean", na.rm = T)
    if (HKflag == T) {
        Hkvartemp <- as.matrix(Hkvar[, -(1:(length(getEDdata) - 
            length(filterED)))])
        class(Hkvartemp) <- "numeric"
        colmedianctr <- apply(Hkvartemp, 2, FUN = "mean")
        for (i in 1:nrow(filterrawdata)) {
            data2ctr[i, ] <- unlist(sapply(1:ncol(filterrawdata), 
                function(zd) {
                  (filterrawdata[i, zd]/colmedianctr[zd]) * (mean(colmedianctr))
                }))
        }
        data2ctrlog <- log2(data2ctr)
        colnames(data2ctrlog) <- colnames(data2log2)
    }
    avgcolsum <- median(colsum)
    for (i in 1:nrow(filterrawdata)) {
        data2GI[i, ] <- unlist(sapply(1:ncol(filterrawdata), 
            function(zd) {
                (filterrawdata[i, zd]/colsum[zd]) * (avgcolsum)
            }))
        data2med[i, ] <- unlist(sapply(1:ncol(filterrawdata), 
            function(zd) {
                (filterrawdata[i, zd]/medofdata[zd]) * mean(medofdata)
            }))
        data2mean[i, ] <- unlist(sapply(1:ncol(filterrawdata), 
            function(zd) {
                (filterrawdata[i, zd]/meanofdata[zd]) * mean(meanofdata)
            }))
    }
    data2GI <- log2(data2GI)
    data2med <- log2(data2med)
    data2mean <- log2(data2mean)
    colnames(data2GI) <- colnames(data2log2)
    colnames(data2med) <- colnames(data2log2)
    colnames(data2mean) <- colnames(data2log2)
    if (HKflag == T) {
        methodlist <- list(data2log2, data2GI, data2med, data2mean, 
            data2ctrlog)
        methodnames <- c("Log2", "TI-G", "MedI-G", "AI-G", "NF-G")
    }
    else {
        methodlist <- list(data2log2, data2GI, data2med, data2mean)
        methodnames <- c("Log2", "TI-G", "MedI-G", "AI-G")
    }
    if (nrow(filterrawdata) > 100) {
        data2vsn <- justvsn((filterrawdata))
        data2quantile <- normalize.quantiles((data2log2), copy = T)
        mediandata <- apply(data2log2, 2, "median", na.rm = T)
        maddata <- apply(data2log2, 2, function(x) mad(x, na.rm = T))
        data2mad <- t(apply(data2log2, 1, function(x) ((x - mediandata)/maddata)))
        data2mad <- data2mad + mean(mediandata)
        data2loess <- normalizeCyclicLoess(data2log2, method = "fast")
        mediandata <- apply(data2log2, 1, "median", na.rm = T)
        flag1 = 1
        for (j in 1:ncol(data2log2)) {
            LRfit <- rlm(as.matrix(data2log2[, j]) ~ mediandata, 
                na.action = na.exclude)
            Coeffs <- LRfit$coefficients
            a <- Coeffs[2]
            b <- Coeffs[1]
            if (flag1 == 1) {
                globalfittedRLR <- (data2log2[, j] - b)/a
                flag1 = 2
            }
            else {
                globalfittedRLR <- cbind(globalfittedRLR, (data2log2[, 
                  j] - b)/a)
            }
        }
        colnames(globalfittedRLR) <- colnames(data2log2)
        {
            x <- 1
            z <- 1
            y <- 1
            flag <- 1
            flag1 <- 1
            data2limloess1 <- NULL
            for (i in 1:length(filterED)) {
                if (x != filterED[i] || i == length(filterED)) {
                  y <- i - 1
                  if (i == length(filterED)) {
                    y <- i
                  }
                  if (flag == 1) {
                    mediandata <- apply(data2log2[, z:y], 1, 
                      "median", na.rm = T)
                    for (j in z:y) {
                      LRfit <- rlm(as.matrix(data2log2[, j]) ~ 
                        mediandata, na.action = na.exclude)
                      Coeffs <- LRfit$coefficients
                      a <- Coeffs[2]
                      b <- Coeffs[1]
                      if (flag1 == 1) {
                        fittedLR <- (data2log2[, j] - b)/a
                        flag1 = 2
                      }
                      else {
                        fittedLR <- cbind(fittedLR, (data2log2[, 
                          j] - b)/a)
                      }
                    }
                    data2vsnrep <- justvsn(as.matrix(filterrawdata[, 
                      z:y]))
                    data2limloess <- normalizeCyclicLoess(data2log2[, 
                      z:y], method = "fast")
                  }
                  if (flag == 2) {
                    mediandata <- apply(data2log2[, z:y], 1, 
                      "median", na.rm = T)
                    for (j in z:y) {
                      LRfit <- rlm(as.matrix(data2log2[, j]) ~ 
                        mediandata, na.action = na.exclude)
                      Coeffs <- LRfit$coefficients
                      a <- Coeffs[2]
                      b <- Coeffs[1]
                      fittedLR <- cbind(fittedLR, (data2log2[, 
                        j] - b)/a)
                    }
                    data2limloess <- cbind(data2limloess, normalizeCyclicLoess(data2log2[, 
                      z:y], method = "fast"))
                    data2vsnrep <- cbind(data2vsnrep, justvsn(as.matrix(filterrawdata[, 
                      z:y])))
                  }
                  z <- i
                  x <- filterED[i]
                  flag = 2
                }
            }
        }
        colnames(fittedLR) <- colnames(data2log2)
        colnames(data2quantile) <- colnames(data2log2)
        if (HKflag == T) {
            methodlist <- list(data2log2, data2limloess, fittedLR, 
                data2vsnrep, data2loess, globalfittedRLR, data2vsn, 
                data2GI, data2med, data2mean, data2ctrlog, data2quantile)
            methodnames <- c("Log2", "Loess-R", "RLR-R", "VSN-R", 
                "Loess-G", "RLR-G", "VSN-G", "TI-G", "MedI-G", 
                "AI-G", "NF-G", "Quantile")
        }
        else {
            methodlist <- list(data2log2, data2limloess, fittedLR, 
                data2vsnrep, data2loess, globalfittedRLR, data2vsn, 
                data2GI, data2med, data2mean, data2quantile)
            methodnames <- c("Log2", "Loess-R", "RLR-R", "VSN-R", 
                "Loess-G", "RLR-G", "VSN-G", "TI-G", "MedI-G", 
                "AI-G", "Quantile")
        }
    }
    jobdir <- paste(getwd(), "/", currentjob[1], sep = "")
    dir.create(jobdir)
    for (i in 1:length(methodlist)) {
        write.table(file = paste(jobdir, "/", methodnames[i], 
            "-normalized.txt", sep = ""), cbind(getrawdata[-(1:2), 
            (1:(length(getEDdata) - length(filterED)))], methodlist[[i]]), 
            sep = "\t", row.names = F, col.names = getrawdata[2, 
                ], quote = F)
    }
    if (HKflag == T) {
        write.table(file = paste(jobdir, "/housekeeping-variables.txt", 
            sep = ""), Hkvar, sep = "\t", row.names = F, col.names = getrawdata[2, 
            ], quote = F)
    }
    write.table(file = paste(jobdir, "/submitted_rawdata.txt", 
        sep = ""), cbind(getrawdata[-(1:2), (1:(length(getEDdata) - 
        length(filterED)))], filterrawdata), sep = "\t", row.names = F, 
        col.names = getrawdata[2, ], quote = F)
    methodlist <- list(methodlist, methodnames, getrawdata, filterrawdata, 
        filterED, HKflag)
    return(methodlist)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
