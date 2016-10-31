\name{analyzeAndPlot}
\alias{analyzeAndPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
analyzeAndPlot
}
\description{
Analyzes normalized data and makes plots. This method is not built for direct execution. To evalute your data use the function "analyzeAndPlotForOneObject()". 
}
\usage{
analyzeAndPlot(normalizeddata, name). 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{normalizeddata}{
%%     ~~Describe \code{normalizeddata} here~~
}
  \item{name}{
%%     ~~Describe \code{name} here~~
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
function (normalizeddata, name) 
{
    {
        currentjob <- name
        if (is.na(currentjob[2])) {
            currentjob[2] = currentjob[1]
        }
        methodlist <- normalizeddata[[1]]
        methodnames <- normalizeddata[[2]]
        getrawdata <- normalizeddata[[3]]
        filterrawdata <- normalizeddata[[4]]
        filterED <- normalizeddata[[5]]
        getEDdata <- ((getrawdata[1, ]))
        HKflag <- normalizeddata[[6]]
        pooledvarmem <- vector()
        Medianofsamples <- vector()
        datamadsamplesmem <- vector()
        datalabels <- vector()
        avgcvmem <- matrix(nrow = length(levels(as.factor(unlist(filterED)))), 
            ncol = length(methodlist), byrow = T)
        avgmadmem <- matrix(nrow = length(levels(as.factor(unlist(filterED)))), 
            ncol = length(methodlist), byrow = T)
        avgvarmem <- matrix(nrow = length(levels(as.factor(unlist(filterED)))), 
            ncol = length(methodlist), byrow = T)
        datamadpeptmem <- matrix(nrow = nrow(filterrawdata), 
            ncol = length(methodlist), byrow = T)
        anpvalue <- anfdr <- kwpvalue <- kwfdr <- vector()
        for (meti in 1:length(methodlist)) {
            datastore <- (methodlist[[meti]])
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
            datalabels <- c(datalabels, templabel)
            x <- 1
            z <- 1
            y <- 1
            flag <- 1
            flag1 <- 1
            count <- 0
            madmem <- matrix(nrow = nrow(datastore), ncol = length(levels(as.factor(unlist(filterED)))), 
                byrow = T)
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
                    madmem[, count] <- apply(datastore[, z:y], 
                      1, function(x) {
                        mad(x, na.rm = T)
                      })
                    nonmissingmat <- (apply(datastore[, z:y], 
                      1, function(x) {
                        ((sum(!is.na(x))))
                      })) - 1
                    tempvar <- nonmissingmat * apply(datastore[, 
                      z:y], 1, function(x) {
                      var(x, na.rm = TRUE)
                    })
                  }
                  if (flag == 2) {
                    count <- count + 1
                    madmem[, count] <- apply(datastore[, z:y], 
                      1, function(x) {
                        mad(x, na.rm = T)
                      })
                    nonmissingmat <- (apply(datastore[, z:y], 
                      1, function(x) {
                        ((sum(!is.na(x))))
                      })) - 1
                    tempvar <- nonmissingmat * apply(datastore[, 
                      z:y], 1, function(x) {
                      var(x, na.rm = TRUE)
                    })
                  }
                  varmem <- c(varmem, ((sum(tempvar, na.rm = T))/(sum(nonmissingmat, 
                    na.rm = T))))
                  z <- i
                  x <- filterED[i]
                  flag = 2
                }
            }
            avgvarmem[, meti] <- varmem
            temmadmatsum <- apply(madmem, 2, mean, na.rm = T)
            avgmadmem[, meti] <- temmadmatsum
            tempcvmat <- matrix(nrow = nrow(datastore), ncol = length(levels(as.factor(unlist(filterED)))), 
                byrow = T)
            for (i in 1:nrow(datastore)) {
                tempcv <- numSummary(datastore[i, ], statistics = c("cv"), 
                  groups = unlist(filterED))
                tempcvmat[i, ] <- tempcv$table
            }
            temcvmatsum <- apply(tempcvmat, 2, mean, na.rm = T)
            avgcvmem[, meti] <- ((temcvmatsum * 100))
            an <- rowSums(is.na(datastore))
            datastoretmp <- datastore[an < (ncol(datastore)/2), 
                ]
            anpvalue <- cbind(anpvalue, apply(datastoretmp, 1, 
                function(x) summary(aov(unlist(x) ~ filterED))[[1]][[5]][1]))
            anfdr <- cbind(anfdr, p.adjust(anpvalue[, meti], 
                method = "BH"))
            kwpvalue <- cbind(kwpvalue, apply(datastoretmp, 1, 
                function(x) kruskal.test(unlist(x) ~ filterED, 
                  na.action = "na.exclude")[[3]][1]))
            kwfdr <- cbind(kwfdr, p.adjust(kwpvalue[, meti], 
                method = "BH"))
        }
        avgcvmempdiff <- sapply(1:ncol(avgcvmem), function(x) (mean(avgcvmem[, 
            x]) * 100)/mean(avgcvmem[, 1]))
        avgmadmempdiff <- sapply(1:ncol(avgmadmem), function(x) (mean(avgmadmem[, 
            x]) * 100)/mean(avgmadmem[, 1]))
        avgvarmempdiff <- sapply(1:ncol(avgvarmem), function(x) (mean(avgvarmem[, 
            x]) * 100)/mean(avgvarmem[, 1]))
    }
    print("Finished analysis, plotting....")
    {
        palette(c("red", "green", "blue", "orange", "darkgray", 
            "blueviolet", "darkslateblue", "darkviolet", "gray", 
            "bisque4", "brown", "cadetblue4", "darkgreen", "darkcyan", 
            "darkmagenta", "darkgoldenrod4", "coral1"))
        jobdir <- paste(getwd(), "/", currentjob[1], sep = "")
        pdf(file = paste(jobdir, "/Norm_report-", currentjob[1], 
            ".pdf", sep = ""), paper = "a4r", width = 0, height = 0)
        def.par <- par(no.readonly = T)
        par(mfrow = c(4, 1))
        textplot(paste("Project Name: ", currentjob[2], sep = ""), 
            halign = "center", valign = "center", cex = 2.5)
        plot(0:10, type = "n", xaxt = "n", yaxt = "n", bty = "n", 
            xlab = "", ylab = "")
        textplot("Normalization report by Normalyzer", halign = "center", 
            valign = "center", cex = 1.5)
        textplot(paste("Report created on: ", Sys.Date(), sep = ""), 
            halign = "center", valign = "center", cex = 1.5)
        pageno = 2
        {
            tout <- rbind(c(1, 2), c(3))
            layout(tout)
            par(mar = c(2, 2, 2, 1), oma = c(2, 2, 3, 2), xpd = NA)
            datacoltotal <- apply(filterrawdata, 2, function(x) {
                sum(x, na.rm = T)
            })
            barplot(datacoltotal, las = 2, main = "Total intensity", 
                cex.names = 0.5)
            datamissingcol <- apply(filterrawdata, 2, function(x) {
                sum(is.na(x))
            })
            barplot(datamissingcol, las = 2, main = "Total missing", 
                cex.names = 0.5)
            mtext("Data Summary", side = 3, adj = 0, outer = TRUE, 
                col = "gray", line = 1)
            mtext(paste("Page ", pageno, sep = ""), side = 1, 
                adj = 1, outer = T, col = "gray", line = 0)
            pageno = pageno + 1
            mtext(paste("Project: ", currentjob[2], sep = ""), 
                side = 1, outer = T, col = "gray", line = 0)
            mtext("Normalyzer Report", side = 1, adj = 0, outer = T, 
                col = "gray", line = 0)
        }
        {
            tout <- rbind(c(1, 2, 3), c(4))
            layout(tout)
            par(mar = c(2, 2, 2, 1), oma = c(2, 2, 3, 2), xpd = NA)
            boxplot(avgcvmem, main = "\%CV - Intragroup", names = c(methodnames), 
                border = "red", density = 20, cex = 0.3, cex.axis = 0.9, 
                las = 2, frame.plot = F)
            stripchart(as.data.frame(avgcvmem), vertical = T, 
                cex = 0.4, las = 2, pch = 20, add = T, col = "darkgray")
            boxplot(avgmadmem, main = "MAD - Intragroup", names = c(methodnames), 
                border = "red", density = 20, cex = 0.3, cex.axis = 0.9, 
                las = 2, frame.plot = F)
            stripchart(as.data.frame(avgmadmem), vertical = T, 
                cex = 0.4, las = 2, pch = 20, add = T, col = "darkgray")
            boxplot(avgvarmem, main = "PV - Intragroup", names = c(methodnames), 
                border = "red", density = 20, cex = 0.3, cex.axis = 0.9, 
                las = 2, frame.plot = F)
            stripchart(as.data.frame(avgvarmem), vertical = T, 
                cex = 0.4, las = 2, pch = 20, add = T, col = "darkgray")
            mtext("Data Summary", side = 3, adj = 0, outer = TRUE, 
                col = "gray", line = 1)
            mtext(paste("Page ", pageno, sep = ""), side = 1, 
                adj = 1, outer = T, col = "gray", line = 0)
            pageno = pageno + 1
            mtext(paste("Project: ", currentjob[2], sep = ""), 
                side = 1, outer = T, col = "gray", line = 0)
            mtext("Normalyzer Report", side = 1, adj = 0, outer = T, 
                col = "gray", line = 0)
        }
        {
            tout <- rbind(c(1, 2, 3), c(4))
            layout(tout)
            par(mar = c(2, 2, 2, 1), oma = c(2, 2, 3, 2), xpd = NA)
            abc <- barplot(avgcvmempdiff, main = "\%CV compared to log2 ", 
                names.arg = c(methodnames), border = "red", ylim = c(min(avgcvmempdiff) - 
                  10, (max(avgmadmempdiff)) + 5), density = 20, 
                cex = 0.9, cex.axis = 0.7, las = 2, xpd = F)
            axis(1, at = c(0.2, (max(abc) + 0.5)), labels = F, 
                lwd.ticks = 0)
            axis(1, at = abc, labels = F, lwd = 0, lwd.ticks = 1)
            text(abc, avgcvmempdiff, labels = round(avgcvmempdiff, 
                digits = 0), pos = 3, las = 2)
            abc <- barplot(avgmadmempdiff, main = "\%MAD compared to log2", 
                names.arg = c(methodnames), border = "red", ylim = c(min(avgmadmempdiff) - 
                  10, (max(avgmadmempdiff)) + 5), density = 20, 
                cex = 0.9, cex.axis = 0.7, las = 2, xpd = F)
            axis(1, at = c(0.2, (max(abc) + 0.5)), labels = F, 
                lwd.ticks = 0)
            axis(1, at = abc, labels = F, lwd = 0, lwd.ticks = 1)
            text(abc, avgmadmempdiff, labels = round(avgmadmempdiff, 
                digits = 0), pos = 3, las = 2)
            abc <- barplot(avgvarmempdiff, main = "\%PV - compared to log2", 
                names.arg = c(methodnames), border = "red", ylim = c(min(avgvarmempdiff) - 
                  10, (max(avgmadmempdiff)) + 5), density = 20, 
                cex = 0.9, cex.axis = 0.7, las = 2, xpd = F)
            axis(1, at = c(0.2, (max(abc) + 0.5)), labels = F, 
                lwd.ticks = 0)
            axis(1, at = abc, labels = F, lwd = 0, lwd.ticks = 1)
            text(abc, avgvarmempdiff, labels = round(avgvarmempdiff, 
                digits = 0), pos = 3, las = 2)
            plot(0:10, type = "n", xaxt = "n", yaxt = "n", bty = "n", 
                xlab = "", ylab = "")
            mtext("Data Summary", side = 3, adj = 0, outer = TRUE, 
                col = "gray", line = 1)
            mtext(paste("Page ", pageno, sep = ""), side = 1, 
                adj = 1, outer = T, col = "gray", line = 0)
            pageno = pageno + 1
            mtext(paste("Project: ", currentjob[2], sep = ""), 
                side = 1, outer = T, col = "gray", line = 0)
            mtext("Normalyzer Report", side = 1, adj = 0, outer = T, 
                col = "gray", line = 0)
        }
        {
            datastore <- methodlist[[1]]
            tempcvmat1 <- matrix(nrow = nrow(datastore), ncol = length(methodlist), 
                byrow = T)
            tempavgmat1 <- matrix(nrow = nrow(datastore), ncol = length(methodlist), 
                byrow = T)
            maxtempcv <- 0
            for (j in 1:length(methodlist)) {
                datastore <- methodlist[[j]]
                datastore <- datastore[, 1:sum(filterED == 1)]
                for (i in 1:nrow(datastore)) {
                  tempcv <- numSummary(datastore[i, ], statistics = c("cv"))
                  tempavg <- numSummary(filterrawdata[i, ], statistics = c("mean"))
                  tempcvmat1[i, j] <- 100 * tempcv$table
                  tempavgmat1[i, j] <- tempavg$table
                }
                if (maxtempcv < max(tempcvmat1, na.rm = T)) {
                  maxtempcv <- max(tempcvmat1, na.rm = T)
                }
            }
            tout <- rbind(c(1, 2, 3, 4), c(5, 6, 7, 8), c(9, 
                10, 11, 12))
            layout(tout)
            par(mar = c(4, 4, 2, 1), oma = c(2, 2, 3, 2), xpd = NA)
            for (i in 1:ncol(tempcvmat1)) {
                plot(tempavgmat1[, i], tempcvmat1[, i], main = methodnames[i], 
                  xlab = "Raw intensity", ylab = "CV", cex = 0.3, 
                  ylim = c(0, maxtempcv))
            }
            mtext("CV vs Raw Intensity plots", side = 3, adj = 0, 
                outer = TRUE, col = "gray", line = 1)
            mtext(paste("Page ", pageno, sep = ""), side = 1, 
                adj = 1, outer = T, col = "gray", line = 0)
            pageno = pageno + 1
            mtext(paste("Project: ", currentjob[2], sep = ""), 
                side = 1, outer = T, col = "gray", line = 0)
            mtext("Normalyzer Report", side = 1, adj = 0, outer = T, 
                col = "gray", line = 0)
        }
        {
            tout <- rbind(c(1, 2, 3, 4), c(5, 6, 7, 8), c(9, 
                10, 11, 12))
            layout(tout)
            par(mar = c(2, 2, 2, 1), oma = c(2, 2, 3, 2), xpd = NA)
            colt <- (c("red", "green", "blue", "orange", "darkgray", 
                "blueviolet", "darkslateblue", "darkviolet", 
                "gray", "bisque4", "brown", "cadetblue4", "darkgreen", 
                "darkcyan", "darkmagenta", "darkgoldenrod4", 
                "coral1"))
            for (j in 1:length(methodlist)) {
                temp <- scale(t(na.omit(methodlist[[j]])), center = TRUE, 
                  scale = TRUE)
                hc <- hclust(dist(temp), "ave")
                plot(as.phylo(hc), main = methodnames[j], cex = 0.5, 
                  tip.color = colt[filterED])
                axisPhylo(side = 1)
            }
            mtext(paste("Dendrograms - Built from", ncol(temp), 
                "variables containing non-missing data", sep = " "), 
                side = 3, adj = 0, outer = TRUE, col = "gray", 
                line = 1)
            mtext(paste("Page ", pageno, sep = ""), side = 1, 
                adj = 1, outer = T, col = "gray", line = 0)
            pageno = pageno + 1
            mtext(paste("Project: ", currentjob[2], sep = ""), 
                side = 1, outer = T, col = "gray", line = 0)
            mtext("Normalyzer Report", side = 1, adj = 0, outer = T, 
                col = "gray", line = 0)
        }
        {
            par(mfrow = c(1, 1))
            tempcorsum <- vector()
            tempcorsumspear <- vector()
            corsumspear <- vector()
            avgcorsum <- matrix(nrow = ((length(levels(as.factor(filterED))) * 
                2)) - 1, ncol = length(methodlist), byrow = T)
            avgcorsumspear <- matrix(nrow = ((length(levels(as.factor(filterED))) * 
                2)) - 1, ncol = length(methodlist), byrow = T)
            corsum <- vector()
            for (i in 1:length(methodlist)) {
                flag1 <- 1
                datastore <- as.matrix(methodlist[[i]])
                for (j in 1:length(filterED)) {
                  if (j != length(filterED) && flag1 == filterED[j] && 
                    flag1 == filterED[j + 1]) {
                    z <- j + 1
                    test <- cor(as.matrix(datastore[, j:z]), 
                      use = "complete.obs", method = "pearson")
                    tempcorsum <- c(tempcorsum, test[2])
                    test <- cor(as.matrix(datastore[, j:z]), 
                      use = "complete.obs", method = "spearman")
                    tempcorsumspear <- c(tempcorsumspear, test[2])
                  }
                  else {
                    corsum <- c(corsum, mean(tempcorsum))
                    tempcorsum <- vector()
                    corsumspear <- c(corsumspear, mean(tempcorsumspear))
                    tempcorsumspear <- vector()
                    if (j != length(filterED)) {
                      z <- j + 1
                      flag1 <- filterED[j]
                      test <- cor(as.matrix(datastore[, j:z]), 
                        use = "complete.obs", method = "pearson")
                      tempcorsum <- c(tempcorsum, test[2])
                      test <- cor(as.matrix(datastore[, j:z]), 
                        use = "complete.obs", method = "spearman")
                      tempcorsumspear <- c(tempcorsumspear, test[2])
                    }
                  }
                }
                avgcorsum[, i] <- corsum
                corsum <- vector()
                avgcorsumspear[, i] <- corsumspear
                corsumspear <- vector()
            }
            tout <- rbind(c(1, 2), c(3))
            layout(tout)
            par(mar = c(2, 2, 2, 1), oma = c(2, 2, 3, 2), xpd = NA)
            abc <- boxplot(avgcorsum, main = "Pearson correlation - Intragroup", 
                names = c(methodnames), border = "red", density = 20, 
                cex = 0.3, cex.axis = 0.9, las = 2)
            stripchart(as.data.frame(avgcorsum), vertical = T, 
                cex = 0.4, las = 2, pch = 20, add = T, col = "darkgreen")
            abc <- boxplot(avgcorsumspear, main = "Spearman correlation - Intragroup", 
                names = c(methodnames), border = "red", density = 20, 
                cex = 0.3, cex.axis = 0.9, las = 2)
            stripchart(as.data.frame(avgcorsumspear), vertical = T, 
                cex = 0.4, las = 2, pch = 20, add = T, col = "darkgreen")
            mtext("Correlation plots", side = 3, adj = 0, outer = TRUE, 
                col = "gray", line = 1)
            mtext(paste("Page ", pageno, sep = ""), side = 1, 
                adj = 1, outer = T, col = "gray", line = 0)
            pageno = pageno + 1
            mtext(paste("Project: ", currentjob[2], sep = ""), 
                side = 1, outer = T, col = "gray", line = 0)
            mtext("Normalyzer Report", side = 1, adj = 0, outer = T, 
                col = "gray", line = 0)
        }
        {
            tout <- rbind(c(1, 2, 3, 4), c(5, 6, 7, 8), c(9, 
                10, 11, 12))
            layout(tout)
            par(mar = c(2, 2, 2, 1), oma = c(3, 2, 3, 2), xpd = NA)
            pc = 0
            for (i in 1:length(methodlist)) {
                datastore <- as.matrix(methodlist[[i]])
                tempcolname <- colnames(datastore)
                datastore <- datastore[, 1:sum(filterED == 1)]
                datastore1 <- datastore[!is.na(datastore[, 1]), 
                  ]
                limma::plotMA(datastore1, array = 1, xlab = "", 
                  ylab = "", main = paste(tempcolname[1], methodnames[i]))
            }
            mtext("MA plots", side = 3, adj = 0, outer = TRUE, 
                col = "gray", line = 1)
            mtext(paste("Project: ", currentjob[2], sep = ""), 
                side = 1, outer = T, col = "gray", line = 0)
            mtext("Normalyzer Report", side = 1, adj = 0, outer = T, 
                col = "gray", line = 0)
            mtext(paste("Page ", pageno, sep = ""), side = 1, 
                adj = 1, outer = T, col = "gray", line = 0)
            pageno = pageno + 1
        }
        {
            tout <- rbind(c(1, 2, 3, 4), c(5, 6, 7, 8), c(9, 
                10, 11, 12))
            layout(tout)
            par(mar = c(2, 2, 2, 1), oma = c(3, 2, 3, 2), xpd = NA)
            pc = 0
            for (i in 1:length(methodlist)) {
                datastore <- methodlist[[i]]
                plot(datastore[, 1], datastore[, 2], xlab = "", 
                  ylab = "", main = methodnames[i], pch = 19, 
                  cex = 0.2)
                fit <- lm(datastore[, 1] ~ datastore[, 2])
                legend("topleft", bty = "n", legend = paste("R2 ", 
                  format(summary(fit)$adj.r.squared, digits = 2)), 
                  cex = 0.7)
                pc = pc + 1
            }
            mtext("Scatterplots", side = 3, adj = 0, outer = TRUE, 
                col = "gray", line = 1)
            mtext(paste("Project: ", currentjob[2], sep = ""), 
                side = 1, outer = T, col = "gray", line = 0)
            mtext("Normalyzer Report", side = 1, adj = 0, outer = T, 
                col = "gray", line = 0)
            mtext(paste("Page ", pageno, sep = ""), side = 1, 
                adj = 1, outer = T, col = "gray", line = 0)
            pageno = pageno + 1
        }
        {
            {
            }
            {
            }
        }
        {
            tout <- rbind(c(1, 2, 3, 4), c(5, 6, 7, 8), c(9, 
                10, 11, 12))
            layout(tout)
            par(mar = c(2, 2, 2, 1), oma = c(3, 2, 3, 2), xpd = NA)
            mindata <- 1000
            maxdata <- 0
            for (i in 1:length(methodlist)) {
                tempmin <- min(methodlist[[i]], na.rm = T)
                tempmax <- max(methodlist[[i]], na.rm = T)
                if (tempmin < mindata) {
                  mindata <- tempmin
                }
                if (tempmax > maxdata) {
                  maxdata <- tempmax
                }
            }
            for (i in 1:length(methodlist)) {
                par(mar = c(5, 1, 1, 1))
                boxplot(methodlist[[i]], cex = 0.1, cex.axis = 0.7, 
                  las = 2, names = colnames(methodlist[[i]]), 
                  main = methodnames[i], col = (filterED), outcol = "lightgray", 
                  ylim = c((mindata - 1), (maxdata + 1)))
            }
            mtext("Boxplots", side = 3, adj = 0, outer = TRUE, 
                col = "gray", line = 1)
            mtext(paste("Page ", pageno, sep = ""), side = 1, 
                adj = 1, outer = T, col = "gray", line = 0)
            pageno = pageno + 1
            mtext(paste("Project: ", currentjob[2], sep = ""), 
                side = 1, outer = T, col = "gray", line = 0)
            mtext("Normalyzer Report", side = 1, adj = 0, outer = T, 
                col = "gray", line = 0)
        }
        tout <- rbind(c(1, 2, 3, 4), c(5, 6, 7, 8), c(9, 10, 
            11, 12))
        layout(tout)
        par(mar = c(2, 2, 2, 1), oma = c(3, 2, 3, 2), xpd = NA)
        pc = 0
        for (i in 1:length(methodlist)) {
            deviations = methodlist[[i]] - rowMedians(methodlist[[i]], 
                na.rm = T)
            boxplot(deviations, outcol = "lightgray", cex = 0.1, 
                cex.axis = 0.7, las = 2, names = colnames(methodlist[[i]]), 
                main = methodnames[i], col = (filterED))
        }
        mtext("Relative Log Expression (RLE) plots", side = 3, 
            adj = 0, outer = TRUE, col = "gray", line = 1)
        mtext(paste("Project: ", currentjob[2], sep = ""), side = 1, 
            outer = T, col = "gray", line = 0)
        mtext("Normalyzer Report", side = 1, adj = 0, outer = T, 
            col = "gray", line = 0)
        mtext(paste("Page ", pageno, sep = ""), side = 1, adj = 1, 
            outer = T, col = "gray", line = 0)
        pageno = pageno + 1
        {
            tout <- rbind(c(1, 2, 3, 4), c(5, 6, 7, 8), c(9, 
                10, 11, 12))
            layout(tout)
            par(mar = c(3, 2, 3, 1), oma = c(3, 2, 3, 2), xpd = NA)
            for (i in 1:length(methodlist)) {
                datastore <- (methodlist[[i]])
                tempd <- density(datastore[, 1], na.rm = T)
                plot(density(datastore[, 1], na.rm = T), xlab = "", 
                  ylab = "", ylim = c(min(tempd$y), max(tempd$y) * 
                    1.5), main = methodnames[i], lty = 2, lwd = 1, 
                  col = "darkgray")
                for (j in 2:ncol(datastore)) {
                  lines(density(datastore[, j], na.rm = T), , 
                    lty = 2, lwd = 1, col = "darkgray")
                }
            }
            mtext("Density plots", side = 3, adj = 0, outer = TRUE, 
                col = "gray", line = 1)
            mtext(paste("Page ", pageno, sep = ""), side = 1, 
                adj = 1, outer = T, col = "gray", line = 0)
            pageno = pageno + 1
            mtext(paste("Project: ", currentjob[2], sep = ""), 
                side = 1, outer = T, col = "gray", line = 0)
            mtext("Normalyzer Report", side = 1, adj = 0, outer = T, 
                col = "gray", line = 0)
        }
        {
            tout <- rbind(c(1, 2, 3, 4), c(5, 6, 7, 8), c(9, 
                10, 11, 12))
            layout(tout)
            par(mar = c(2, 2, 2, 1), oma = c(3, 2, 3, 2), xpd = NA)
            pc = 0
            for (i in 1:length(methodlist)) {
                datastore <- (methodlist[[i]])
                tempcolname <- colnames(datastore)
                qqnorm(datastore[, j], main = paste(tempcolname[j], 
                  methodnames[i]), xlab = "", ylab = "")
            }
            mtext("Q-Q plots", side = 3, adj = 0, outer = TRUE, 
                col = "gray", line = 1)
            mtext(paste("Project: ", currentjob[2], sep = ""), 
                side = 1, outer = T, col = "gray", line = 0)
            mtext("Normalyzer Report", side = 1, adj = 0, outer = T, 
                col = "gray", line = 0)
            mtext(paste("Page ", pageno, sep = ""), side = 1, 
                adj = 1, outer = T, col = "gray", line = 0)
            pageno = pageno + 1
        }
        tout <- rbind(c(1, 2, 3, 4), c(5, 6, 7, 8), c(9, 10, 
            11, 12))
        layout(tout)
        par(mar = c(2, 2, 2, 1), oma = c(3, 2, 3, 2), xpd = NA)
        for (i in 1:length(methodlist)) {
            datastore <- (methodlist[[i]])
            d <- dist(scale(t(na.omit(datastore)), center = TRUE, 
                scale = TRUE))
            fit <- cmdscale(d, eig = TRUE, k = 2)
            x <- fit$points[, 1]
            y <- fit$points[, 2]
            plot(x, y, type = "n", main = methodnames[i], xlab = "", 
                ylab = "")
            text((fit$points[, 1]), (fit$points[, 2]), col = filterED)
        }
        mtext(paste("MDS plots - Built from", ncol(temp), "variables with non-missing data", 
            sep = " "), side = 3, adj = 0, outer = TRUE, col = "gray", 
            line = 1)
        mtext(paste("Project: ", currentjob[2], sep = ""), side = 1, 
            outer = T, col = "gray", line = 0)
        mtext("Normalyzer Report", side = 1, adj = 0, outer = T, 
            col = "gray", line = 0)
        mtext(paste("Page ", pageno, sep = ""), side = 1, adj = 1, 
            outer = T, col = "gray", line = 0)
        pageno = pageno + 1
        tout <- rbind(c(1, 2, 3, 4), c(5, 6, 7, 8), c(9, 10, 
            11, 12))
        layout(tout)
        par(mar = c(2, 2, 2, 1), oma = c(3, 2, 3, 2), xpd = NA)
        for (i in 1:length(methodlist)) {
            datastore <- (methodlist[[i]])
            meanSdPlot(datastore, xlab = "", ylab = "", main = methodnames[i])
        }
        mtext("MeanSDplots", side = 3, adj = 0, outer = TRUE, 
            col = "gray", line = 1)
        mtext(paste("Project: ", currentjob[2], sep = ""), side = 1, 
            outer = T, col = "gray", line = 0)
        mtext("Normalyzer Report", side = 1, adj = 0, outer = T, 
            col = "gray", line = 0)
        mtext(paste("Page ", pageno, sep = ""), side = 1, adj = 1, 
            outer = T, col = "gray", line = 0)
        pageno = pageno + 1
        {
            tout <- rbind(c(1, 2, 3), c(4))
            layout(tout)
            par(mar = c(2, 2, 2, 1), oma = c(2, 2, 3, 2), xpd = NA)
            barplot(colSums(anfdr < 0.05), main = "ANOVA", names = c(methodnames), 
                border = "red", density = 20, cex = 0.5, cex.axis = 0.9, 
                las = 2, ylab = "No. of Variables with FDR<0.05")
            barplot(colSums(kwfdr < 0.05), main = "Kruskal Wallis", 
                names = c(methodnames), border = "red", density = 20, 
                cex = 0.5, cex.axis = 0.9, las = 2, ylab = "No. of Variables with FDR<0.05")
            mtext("Differential Expression", side = 3, adj = 0, 
                outer = TRUE, col = "gray", line = 1)
            mtext(paste("Page ", pageno, sep = ""), side = 1, 
                adj = 1, outer = T, col = "gray", line = 0)
            pageno = pageno + 1
            mtext(paste("Project: ", currentjob[2], sep = ""), 
                side = 1, outer = T, col = "gray", line = 0)
            mtext("Normalyzer Report", side = 1, adj = 0, outer = T, 
                col = "gray", line = 0)
        }
        dev.off()
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
