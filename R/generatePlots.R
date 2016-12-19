#' Generate full output report plot document
#' 
#' @param nr Normalyzer results object.
#' @param jobdir Path to output directory for run.
#' @export
#' @examples
#' generatePlots
#' generatePlots(my_normalyzer_results, "path/to/my_run")
generatePlots <- function(nr, jobdir) {
    
    nds <- nr@nds
    currentjob <- nds@jobName

    currentLayout <- grid.layout(nrow=5, ncol=6, heights=c(0.1, 1, 1, 1, 0.1), 
                                 widths=c(0.1, 1, 1, 1, 1, 0.1), 
                                 default.units=c('null', 'null'))
    currentFont <- "Helvetica"

    setupPlotting(currentjob, jobdir)
    plotFrontPage(currentjob, currentFont)
    
    # TI
    print("DEBUG: Plotting page 2")
    pageno <- 2
    plotSampleOutlierSummary(nr, currentLayout, pageno)
    
    # CV
    print("DEBUG: Plotting page 3")
    pageno <- pageno + 1
    plotReplicateVariance(nr, currentLayout, pageno)

    # Stable variables plot and CV in percent difference
    print("DEBUG: Plotting page 4")
    pageno <- pageno + 1
    plotReplicateVarAndStableVariables(nr, currentLayout, pageno)

    # CVvsintensityplot 
    print("DEBUG: Plotting page 5")
    pageno <- pageno + 1
    plotCVvsIntensity(nr, currentLayout, pageno)

    #MA plots
    print("DEBUG: Plotting page 6!")
    pageno <- pageno + 1
    plotMA(nr, currentLayout, pageno)
    
    # Scatterplots
    print("DEBUG: Plotting page 7")
    pageno <- pageno + 1
    plotScatter(nr, currentLayout, pageno)

    # QQplot
    print("DEBUG: Plotting page 8")
    pageno <- pageno + 1
    plotQQ(nr, currentLayout, pageno)

    # Boxplot
    print("DEBUG: Plotting page 9")
    pageno <- pageno + 1
    plotBoxPlot(nr, currentLayout, pageno)

    # RLE plots
    print("DEBUG: Plotting page 10")
    pageno <- pageno + 1
    plotRLE(nr, currentLayout, pageno)
    
    # Density plots
    print("DEBUG: Plotting page 11")
    pageno <- pageno + 1
    plotDensity(nr, currentLayout, pageno)

    # MDS plot
    print("DEBUG: Plotting page 12")
    pageno <- pageno + 1
    plotMDS(nr, currentLayout, pageno)
    
    # meanSDplot
    print("DEBUG: Plotting page 13")
    pageno <- pageno + 1
    plotMeanSD(nr, currentLayout, pageno)
    
    # Correlation
    print("DEBUG: Plotting page 15")
    pageno <- pageno + 1
    plotCorrelation(nr, currentLayout, pageno)

    # Dendrograms
    print("DEBUG: Plotting page 16")
    pageno <- pageno + 1
    plotDendrograms(nr, currentLayout, pageno)

    # DE plots
    print("DEBUG: Plotting page 17")
    pageno <- pageno + 1
    plotDEPlots(nr, currentLayout, pageno)
    
    dev.off()
}

#' Setup PDF report settings
#' 
#' @param currentjob Name of current run.
#' @param jobdir Path to output directory for run.
setupPlotting <- function(currentjob, jobdir) {

    palette(c("red", "green", "blue", "orange", "darkgray", "blueviolet", "darkslateblue", "darkviolet", "gray", "bisque4",
              "brown", "cadetblue4", "darkgreen", "darkcyan", "darkmagenta", "darkgoldenrod4", "coral1"))
    
    pdf(file=paste(jobdir, "/Norm_report-", currentjob, ".pdf", sep=""), paper="a4r", width=0, height=0)    
    theme_norm <- theme_set(theme_bw())
    theme_norm <- theme_update(panel.grid.minor=element_blank(), axis.text=element_text(size=7), axis.title=element_text(size=8), 
                               plot.title=element_text(size=8), plot.margin=unit(c(1, 1, 1, 1), "mm"))
    def.par <- par(no.readonly=T)
}

#' Generate first page in output report
#' 
#' @param currentjob Name of current run.
#' @param currentFont Font used for output document.
plotFrontPage <- function(currentjob, currentFont) {
    
    par(mfrow=c(4, 1))
    # TODO: Re-insert nice illustration (figure?)
    # data(data4pdftitle)
    plot(1, type="n", axes=F, xlab="", ylab="")
    
    # TODO: Re-insert nice illustration (figure?)
    # boxplot(data4pdftitle, axes=F, col=c("green","green","red","red"))
    la1 <- grid.layout(nrow=7, ncol=1, heights=c(0.2, 1, 0.1, 0.2, 0.1, 0.2, 0.2), default.units=c('null','null'))
    gpfill <- gpar(fill="gray90", lwd=0, lty=0)
    pushViewport(viewport(layout=la1))
    grid.rect(vp=viewport(layout.pos.row=1), gp=gpfill)
    grid.rect(vp=viewport(layout.pos.row=7), gp=gpfill)
    
    grid.text(paste("Project Name: ", currentjob, sep=""), vp=viewport(layout.pos.row=3), just=c("center","center"), 
              gp=gpar(fontsize=12, fontfamily=currentFont, col="black"))
    grid.text(paste("Normalyzer (ver 1.1.1)"), vp=viewport(layout.pos.row=4), just=c("center", "center"), 
              gp=gpar(fontface="bold", fontsize=32, fontfamily=currentFont, col="darkblue"))
    
    grid.text(paste("Report created on: ", Sys.Date(), sep=""), vp=viewport(layout.pos.row=5), just=c("center","center"),
              gp=gpar(fontsize=12, fontfamily=currentFont, col="black"))
    
    citationText <- "Citation: Chawade, A., Alexandersson, E., Levander, F. (2014). Normalyzer: a tool for rapid evaluation of normalization methods for omics data sets. J Proteome Res.,13 (6)"
    
    grid.text(citationText, vp=viewport(layout.pos.row=6), just=c("center","center"), gp=gpar(fontsize=10, fontfamily=currentFont, col="black"))
    
    grid.text("Documentation for analyzing this report can be found at http://quantitativeproteomics.org/normalyzer/help.php",
              vp=viewport(layout.pos.row=7), just=c("center", "center"), gp=gpar(fontsize=10, fontfamily=currentFont, col="black"))
}

#' Generate sample summary of intensities, missing values and MDS plot
#' 
#' @param nr Normalyzer results object.
#' @param currentLayout Layout used for document.
#' @param pageno Current page number.
plotSampleOutlierSummary <- function(nr, currentLayout, pageno) {

    nds <- nr@nds
    methodlist <- getNormalizationMatrices(nr)
    filterED <- nds@sampleReplicateGroups
    filterrawdata <- nds@filterrawdata
    currentjob <- nds@jobName
    
    tout <- rbind(c(1, 2), c(3, 4))
    layout(tout)
    par(mar=c(4, 4, 2, 1), oma=c(2, 2, 3, 2), xpd=NA)
    datacoltotal <- apply(filterrawdata, 2, function(x) { sum(x, na.rm=T) })
    
    barplot(datacoltotal, las=2, main="Total intensity", cex.names=0.5, names.arg=substr(names(datacoltotal), 1, 10))
    datamissingcol <- apply(filterrawdata, 2, function(x) { sum(is.na(x)) })
    barplot(datamissingcol, las=2, main="Total missing", cex.names=0.5, names.arg=substr(names(datamissingcol), 1, 10))
    datastore <- methodlist[[1]]
    
    d <- dist(scale(t(na.omit(datastore)), center=TRUE, scale=TRUE))
    fit <- cmdscale(d, eig=TRUE, k=2)
    x <- fit$points[, 1]
    y <- fit$points[, 2]
    
    plot(x, y, type="n", main="Log2-MDS plot", xlab="", ylab="")
    text(fit$points[, 1], fit$points[, 2], col=filterED, labels=filterED)
    pushViewport(viewport(layout=currentLayout))
    printMeta("Data Summary - Outlier detection", pageno, currentjob)
}

#' Generate normalization replicate variance summary
#' 
#' @param nr Normalyzer results object.
#' @param currentLayout Layout used for document.
#' @param pageno Current page number.
plotReplicateVariance <- function(nr, currentLayout, pageno) {

    nds <- nr@nds
    methodnames <- getUsedMethodNames(nr)
    currentjob <- nds@jobName
    
    ner <- nr@ner
    avgcvmem <- ner@avgcvmem
    avgmadmem <- ner@avgmadmem
    avgvarmem <- ner@avgvarmem
    
    tout <- rbind(c(1, 2, 3), c(4))
    layout(tout)
    par(mar=c(2, 2, 2, 1), oma=c(2, 2, 3, 2), xpd=NA)
    
    boxplot(avgcvmem, main="PCV - Intragroup", names=c(methodnames), border="red", density=20, cex=0.3, cex.axis=0.9, las=2, frame.plot=F)
    stripchart(as.data.frame(avgcvmem), vertical=T, cex=0.4, las=2, pch=20, add=T, col="darkgray")
    boxplot(avgmadmem, main="PMAD - Intragroup", names=c(methodnames), border="red", density=20, cex=0.3, cex.axis=0.9, las=2, frame.plot=F)
    stripchart(as.data.frame(avgmadmem), vertical=T, cex=0.4, las=2, pch=20, add=T, col="darkgray")
    boxplot(avgvarmem, main="PEV - Intragroup", names=c(methodnames), border="red", density=20, cex=0.3, cex.axis=0.9, las=2, frame.plot=F)
    stripchart(as.data.frame(avgvarmem), vertical=T, cex=0.4, las=2, pch=20, add=T, col="darkgray")
    
    pushViewport(viewport(layout=currentLayout))
    printMeta("Replicate variation", pageno, currentjob)
}

#' Generate replicate variance overview for normalization methods and CV
#' plot for stable variables
#' 
#' @param nr Normalyzer results object.
#' @param currentLayout Layout used for document.
#' @param pageno Current page number.
plotReplicateVarAndStableVariables <- function(nr, currentLayout, pageno) {

    nds <- nr@nds
    methodnames <- getUsedMethodNames(nr)
    currentjob <- nds@jobName
    
    ner <- nr@ner
    nonsiganfdrlistcvpdiff <- ner@nonsiganfdrlistcvpdiff
    avgcvmempdiff <- ner@avgcvmempdiff
    avgmadmempdiff <- ner@avgmadmempdiff
    avgvarmempdiff <- ner@avgvarmempdiff
    
    tout <- rbind(c(1, 2, 3), c(4, 5, 5))
    layout(tout)
    par(mar=c(6, 6, 3, 1), oma=c(2, 3, 3, 2), xpd=NA)
    abc <- barplot(avgcvmempdiff, main="PCV compared to log2 ", names.arg=c(methodnames), border="red", 
                   ylim=c(min(avgcvmempdiff) - 10, (max(avgmadmempdiff)) + 5), density=20, cex=0.9, cex.axis=0.7, las=2, xpd=F)
    axis(1, at=c(0.2, (max(abc) + 0.5)), labels=F, lwd.ticks=0)
    axis(1, at=abc, labels=F, lwd=0, lwd.ticks=1)
    text(abc,avgcvmempdiff, labels=round(avgcvmempdiff, digits=0), pos=3, las=2)  
    abc<-barplot(avgmadmempdiff, main="PMAD compared to log2", names.arg=c(methodnames), border="red", 
                 ylim=c(min(avgmadmempdiff) - 10, max(avgmadmempdiff) + 5), density=20, cex=0.9, cex.axis=0.7, las=2, xpd=F)
    axis(1, at=c(0.2, max(abc) + 0.5), labels=F, lwd.ticks=0)
    axis(1, at=abc, labels=F, lwd=0, lwd.ticks=1)
    text(abc, avgmadmempdiff, labels=round(avgmadmempdiff, digits=0), pos=3, las=2)
    abc <- barplot(avgvarmempdiff, main="%PEV - compared to log2", names.arg=c(methodnames), border="red", 
                   ylim=c(min(avgvarmempdiff) - 10, max(avgmadmempdiff) + 5), density=20, cex=0.9, cex.axis=0.7, las=2, xpd=F)
    axis(1, at=c(0.2, max(abc) + 0.5), labels=F, lwd.ticks=0)
    axis(1, at=abc, labels=F, lwd=0, lwd.ticks=1)
    text(abc, avgvarmempdiff, labels=round(avgvarmempdiff, digits=0), pos=3, las=2)
    
    if (!is.na(nonsiganfdrlistcvpdiff)) {
        
        if (min(avgcvmempdiff) < 0 || max(avgcvmempdiff) > 100 || min(nonsiganfdrlistcvpdiff) < 0 || max(nonsiganfdrlistcvpdiff) > 100) {
            plot(avgcvmempdiff, nonsiganfdrlistcvpdiff, pch=18, xlim=c(0, 100), ylim=c(0, 100), main="Stable variables plot", 
                 xlab="PCV (intragroup) compared to Log2", ylab="% Global CV of stable variables compared to Log2")
            showLabels(avgcvmempdiff, nonsiganfdrlistcvpdiff, labels=methodnames, id.method="mahal", id.cex=0.7, id.col="black")
        }
        else {
            plot(avgcvmempdiff, nonsiganfdrlistcvpdiff, pch=18, main="Stable variables plot", 
                 xlab="PCV (Intragroup) compared to Log2", ylab="% Global CV of stable variables compared to Log2")
            showLabels(avgcvmempdiff, nonsiganfdrlistcvpdiff, labels=methodnames, id.method="mahal", id.cex=0.7, id.col="black")
        }
    }
    
    pushViewport(viewport(layout=currentLayout))
    printMeta("Replicate variation (Relative to Log2)", pageno, currentjob)
}

#' Overview of coefficient of variation compared to intensity for variables
#'  for different normalization methods
#' 
#' @param nr Normalyzer results object.
#' @param currentLayout Layout used for document.
#' @param pageno Current page number.
plotCVvsIntensity <- function(nr, currentLayout, pageno) {

    nds <- nr@nds
    methodnames <- getUsedMethodNames(nr)
    methodlist <- getNormalizationMatrices(nr)
    currentjob <- nds@jobName
    sampleReplicateGroups <- nds@sampleReplicateGroups
    filterrawdata <- nds@filterrawdata
    
    datastore <- methodlist[[1]]
    tempcvmat1 <- matrix(nrow=nrow(datastore), ncol=length(methodlist), byrow=T)
    tempavgmat1 <- matrix(nrow=nrow(datastore), ncol=length(methodlist), byrow=T)
    maxtempcv <- 0
    
    for (j in 1:length(methodlist)) {
        
        datastore <- methodlist[[j]]
        datastore <- datastore[, 1:sum(sampleReplicateGroups == 1)]
        
        for (i in 1:nrow(datastore)) {
            
            tempcv <- numSummary(datastore[i, ], statistics=c("cv"))
            tempavg <- numSummary(filterrawdata[i, ], statistics=c("mean"))
            tempcvmat1[i, j] <- 100 * tempcv$table
            tempavgmat1[i, j] <- tempavg$table 
        }
        
        if (maxtempcv < max(tempcvmat1, na.rm=T)) {
            maxtempcv <- max(tempcvmat1, na.rm=T)
        }
    }
    
    tout <- rbind(c(1, 2, 3, 4), c(5, 6, 7, 8), c(9, 10, 11, 12))
    layout(tout)
    par(mar=c(4, 4, 2, 1), oma=c(2, 2, 3, 2), xpd=NA)
    
    for (i in 1:ncol(tempcvmat1)) {
        plot(tempavgmat1[, i], tempcvmat1[, i], main=methodnames[i], xlab="Raw intensity", ylab="CV", cex=0.3, ylim=c(0, maxtempcv))
    }
    
    pushViewport(viewport(layout=currentLayout))
    printMeta("CV vs Raw Intensity plots", pageno, currentjob)
}

#' Expression vs. fold-change for variables for each normalization method
#' ! TODO: Need to investigate what the actual fold change is here
#' 
#' @param nr Normalyzer results object.
#' @param currentLayout Layout used for document.
#' @param pageno Current page number.
plotMA <- function(nr, currentLayout, pageno) {

    nds <- nr@nds
    methodnames <- getUsedMethodNames(nr)
    methodlist <- getNormalizationMatrices(nr)
    currentjob <- nds@jobName
    filterED <- nds@sampleReplicateGroups
    filterrawdata <- nds@filterrawdata
    
    Malist <- list()
    for (i in 1:length(methodlist)) {
        
        datastore <- as.matrix(methodlist[[i]])
        tempcolname <- colnames(datastore)
        datastore <- datastore[, 1:sum(filterED == 1)]
        datastore1 <- datastore[!is.na(datastore[, 1]), ]
        avg <- rowMeans(datastore1)
        fc <- apply(cbind(datastore1[, 1], avg), 1, function(x) x[1] - x[2])
        df <- as.data.frame(cbind(avg, fc))
        
        Malist[[i]] <- ggplot(df, aes(avg, fc)) + 
            geom_point(color="darkgray", size=0.7) + 
            labs(x=("Replicate group mean"), y=("Replicate-1 Fold Change"), 
                 title=(paste(tempcolname[1], methodnames[i]))) + 
            stat_smooth(method="loess", se=F, colour="red") + 
            geom_abline(intercept=0, slope=0, size=0.3)
    } 
    
    grid.newpage()
    pushViewport(viewport(layout=currentLayout))
    
    printPlots(Malist, "MA plots", pageno, currentjob)
}

#' Scatter plot comparing replicate expression between samples
#' TODO: Investigate this further
#' 
#' @param nr Normalyzer results object.
#' @param currentLayout Layout used for document.
#' @param pageno Current page number.
plotScatter <- function(nr, currentLayout, pageno) {

    nds <- nr@nds
    methodnames <- getUsedMethodNames(nr)
    methodlist <- getNormalizationMatrices(nr)
    currentjob <- nds@jobName
    
    tout <- rbind(c(1, 2, 3, 4), c(5, 6, 7, 8), c(9, 10, 11, 12))
    layout(tout)
    par(mar=c(2, 2, 2, 1), oma=c(3, 2, 3, 2), xpd=NA)
    
    for (i in 1:length(methodlist)) {
        datastore <- methodlist[[i]]
        fit <- lm(datastore[, 1]~datastore[, 2])
        plot(datastore[, 1], datastore[, 2], xlab="", ylab="", main=methodnames[i], pch=19, cex=0.2)
        
        legend("topleft", bty="n", legend=paste("R2 ", format(summary(fit)$adj.r.squared, digits=2)), cex=0.7)
    }
    pushViewport(viewport(layout=currentLayout))
    printMeta("Scatterplots", pageno, currentjob)
}

#' QQ-plots for variable values
#' TODO: Investigate this further - What is actually extracted here?
#' 
#' @param nr Normalyzer results object.
#' @param currentLayout Layout used for document.
#' @param pageno Current page number.
plotQQ <- function(nr, currentLayout, pageno) {

    nds <- nr@nds
    methodnames <- getUsedMethodNames(nr)
    methodlist <- getNormalizationMatrices(nr)
    currentjob <- nds@jobName

    qqlist <- list()
    
    for (i in 1:length(methodlist)) {  
        datastore <- methodlist[[i]]
        tempcolname <- colnames(datastore)
        #qqnorm(datastore[,1],main=paste(tempcolname[1],methodnames[i]),xlab="",ylab="")
        qqlist[[i]] <- qplot(sample=datastore[, 1], stat="qq") + labs(x=(""), y=(""), title=(paste(tempcolname[1], methodnames[i])))
    }
    
    grid.newpage()
    pushViewport(viewport(layout=currentLayout))
    printPlots(qqlist, "Q-Q plots", pageno, currentjob)
}

#' Boxplots showing distribution of values after different normalizations
#' TODO: Investigate this further - What is actually extracted here?
#' 
#' @param nr Normalyzer results object.
#' @param currentLayout Layout used for document.
#' @param pageno Current page number.
plotBoxPlot <- function(nr, currentLayout, pageno) {

    nds <- nr@nds
    methodnames <- getUsedMethodNames(nr)
    methodlist <- getNormalizationMatrices(nr)
    currentjob <- nds@jobName
    filterED <- nds@sampleReplicateGroups
    filterrawdata <- nds@filterrawdata
    
    tout<-rbind(c(1, 2, 3, 4), c(5, 6, 7, 8), c(9, 10, 11, 12))
    layout(tout)
    par(mar=c(2, 2, 2, 1), oma=c(3, 2, 3, 2), xpd=NA)
    mindata <- 1000
    maxdata <- 0
    
    for (i in 1:length(methodlist)) { 
        tempmin <- min(methodlist[[i]], na.rm=T)
        tempmax <- max(methodlist[[i]], na.rm=T)
        if (tempmin < mindata) {
            mindata <- tempmin
        }
        if (tempmax > maxdata) {
            maxdata <- tempmax
        }
    }
    
    for (i in 1:length(methodlist)) {   
        par(mar=c(5, 1, 1, 1))
        boxplot(methodlist[[i]], cex=0.1, cex.axis=0.7, las=2, main=methodnames[i], col=(filterED), outcol="lightgray", 
                ylim=c(mindata - 1, maxdata + 1), names=substr(colnames(methodlist[[i]]), 1, 10))
    }
    
    pushViewport(viewport(layout=currentLayout))
    printMeta("Boxplots", pageno, currentjob)
}

#' Boxplots showing relative log expression after normalizations
#' TODO: Investigate this further - What is actually extracted here?
#' 
#' @param nr Normalyzer results object.
#' @param currentLayout Layout used for document.
#' @param pageno Current page number.
plotRLE <- function(nr, currentLayout, pageno) {

    nds <- nr@nds
    methodnames <- getUsedMethodNames(nr)
    methodlist <- getNormalizationMatrices(nr)
    currentjob <- nds@jobName
    filterED <- nds@sampleReplicateGroups
    
    tout <- rbind(c(1, 2, 3, 4), c(5, 6, 7, 8), c(9, 10, 11, 12))
    layout(tout)
    par(mar=c(2, 2, 2, 1), oma=c(3, 2, 3, 2), xpd=NA)
    
    for (i in 1:length(methodlist)) {
        deviations = methodlist[[i]] - rowMedians(methodlist[[i]], na.rm=T)
        boxplot(deviations, outcol="lightgray", cex=0.1, cex.axis=0.7, las=2, main=methodnames[i], col=(filterED), names=substr(colnames(methodlist[[i]]), 1, 6))
    }
    
    pushViewport(viewport(layout=currentLayout))
    printMeta("Relative Log Expression (RLE) plots", pageno, currentjob)
}

#' Density plots showing value distributions after normalizations
#' 
#' @param nr Normalyzer results object.
#' @param currentLayout Layout used for document.
#' @param pageno Current page number.
plotDensity <- function(nr, currentLayout, pageno) {

    nds <- nr@nds
    methodnames <- getUsedMethodNames(nr)
    methodlist <- getNormalizationMatrices(nr)
    currentjob <- nds@jobName

    tout <- rbind(c(1, 2, 3, 4), c(5, 6, 7, 8), c(9, 10, 11, 12))
    layout(tout)
    par(mar=c(3, 2, 3, 1), oma=c(3, 2, 3, 2), xpd=NA)
    
    for (i in 1:length(methodlist)) {
        datastore <- (methodlist[[i]])
        tempd <- density(datastore[, 1], na.rm=T)
        plot(density(datastore[, 1], na.rm=T), xlab="", ylab="", ylim=c(min(tempd$y), max(tempd$y) * 1.5), main=methodnames[i], lty=2, lwd=1, col="darkgray")
        
        for (j in 2:ncol(datastore)) {
            lines(density(datastore[, j], na.rm=T), , lty=2, lwd=1, col="darkgray")
        }
    }
    
    pushViewport(viewport(layout=currentLayout))
    printMeta("Density plots", pageno, currentjob)
}

#' MDS plots showing grouping of samples after normalizations
#' TODO: Investigate this further - What is actually extracted here?
#' TODO: Do the d-value with number non-missing data vanish?
#' 
#' @param nr Normalyzer results object.
#' @param currentLayout Layout used for document.
#' @param pageno Current page number.
plotMDS <- function(nr, currentLayout, pageno) {

    nds <- nr@nds
    methodnames <- getUsedMethodNames(nr)
    methodlist <- getNormalizationMatrices(nr)
    currentjob <- nds@jobName
    filterED <- nds@sampleReplicateGroups

    tout <- rbind(c(1, 2, 3, 4), c(5, 6, 7, 8), c(9, 10, 11, 12))
    layout(tout)
    par(mar=c(2, 2, 2, 1), oma=c(3, 2, 3, 2), xpd=NA)

    for (i in 1:length(methodlist)) {
        datastore <- (methodlist[[i]])
        d <- dist(scale(t(na.omit(datastore)), center=TRUE, scale=TRUE))
        fit <- cmdscale(d, eig=TRUE, k=2)
        x <- fit$points[, 1]
        y <- fit$points[, 2]
        plot(x, y, type="n", main=methodnames[i], xlab="", ylab="")
        text(fit$points[, 1], fit$points[, 2], col=filterED, labels=filterED)
    }
    
    pushViewport(viewport(layout=currentLayout))
    printMeta(paste("MDS plots - Built from", ncol(d), "variables with non-missing data", sep=" "), pageno, currentjob)
}

#' Visualize standard deviation over (expression?) for different values
#' TODO: What is on the x-axis here?
#' 
#' @param nr Normalyzer results object.
#' @param currentLayout Layout used for document.
#' @param pageno Current page number.
plotMeanSD <- function(nr, currentLayout, pageno) {

    nds <- nr@nds
    methodnames <- getUsedMethodNames(nr)
    methodlist <- getNormalizationMatrices(nr)
    currentjob <- nds@jobName

    sdPlots <- list()
    
    for (i in 1:length(methodlist)) {
        
        datastore <- methodlist[[i]]
        print(paste("Feeding method name: ", methodnames[i]))
        
        msd <- meanSdPlot(datastore, xlab="", ylab="", plot=F)
        
        sdPlots[[i]] <- msd$gg + 
            ggtitle(methodnames[i]) + 
            theme(legend.position="none", plot.margin=unit(c(1,0,0,0), "cm"))

        # TODO: The main=methodnames[i] seemed to cause crash here // Jakob
        # meanSdPlot(datastore, xlab="", ylab="", main=methodnames[i])
    }

    
    do.call("grid.arrange", c(sdPlots, nrow=3, ncol=4, padding=1, top=1))

    pushViewport(viewport(layout=currentLayout))
    printMeta("MeanSDplots", pageno, currentjob)
}

#' Visualize correlations for plots
#' TODO: What exactly is going on here?
#' 
#' @param nr Normalyzer results object.
#' @param currentLayout Layout used for document.
#' @param pageno Current page number.
plotCorrelation <- function(nr, currentLayout, pageno) {

    ner <- nr@ner
    avgpercorsum <- ner@avgpercorsum
    avgspecorsum <- ner@avgspecorsum
    
    nds <- nr@nds
    methodnames <- getUsedMethodNames(nr)
    methodlist <- getNormalizationMatrices(nr)
    currentjob <- nds@jobName
    filterED <- nds@sampleReplicateGroups
    filterrawdata <- nds@filterrawdata
    
    tout <- rbind(c(1, 2), c(3))
    layout(tout)
    par(mar=c(2, 2, 2, 1), oma=c(2, 2, 3, 2), xpd=NA)
    
    perdf <- data.frame(matrix(unlist(avgpercorsum), 
                               nrow=as.numeric(max(summary(avgpercorsum)[1])), 
                               byrow=T))
    
    abc <- boxplot(perdf, main="Pearson correlation - Intragroup", 
                   names=c(methodnames), border="red", 
                   density=20, cex=0.3, cex.axis=0.9, las=2)
    
    stripchart(as.data.frame(perdf), vertical=T, cex=0.4, las=2, pch=20, add=T, col="darkgreen")
    
    spedf <- data.frame(matrix(unlist(avgspecorsum), 
                               nrow=as.numeric(max(summary(avgspecorsum)[1])), 
                               byrow=T))
    
    abc <- boxplot(spedf, main="Spearman correlation - Intragroup", 
                   names=c(methodnames), border="red", 
                   density=20, cex=0.3, cex.axis=0.9, las=2)
    
    stripchart(as.data.frame(spedf), 
               vertical=T, cex=0.4, las=2, pch=20, add=T, col="darkgreen")
    pushViewport(viewport(layout=currentLayout))
    printMeta("Correlation plots", pageno, currentjob)
}

#' Visualize dendrogram grouping of samples
#' TODO: Why is the Quantile one crazy? (Numbers instead of replicate names)
#' 
#' @param nr Normalyzer results object.
#' @param currentLayout Layout used for document.
#' @param pageno Current page number.
plotDendrograms <- function(nr, currentLayout, pageno) {

    nds <- nr@nds
    methodnames <- getUsedMethodNames(nr)
    methodlist <- getNormalizationMatrices(nr)
    currentjob <- nds@jobName
    filterED <- nds@sampleReplicateGroups

    print(length(methodnames))
    print(length(filterED))
    
    tout <- rbind(c(1, 2, 3, 4), c(5, 6, 7, 8), c(9, 10, 11, 12))
    layout(tout)
    par(mar=c(2, 2, 2, 1), oma=c(2, 2, 3, 2), xpd=NA)
    colt<-(c("red", "green", "blue", "orange", "darkgray", "blueviolet", 
             "darkslateblue", "darkviolet", "gray", "bisque4", "brown", 
             "cadetblue4", "darkgreen", "darkcyan", "darkmagenta", 
             "darkgoldenrod4", "coral1"))
    
    for (j in 1:length(methodlist)) {
        
        scaledTransposedMatrix <- scale(t(na.omit(methodlist[[j]])), center=TRUE, scale=TRUE)
        
        hc <- hclust(dist(scaledTransposedMatrix), "ave")
        
        plot(as.phylo(hc), main=methodnames[j], cex=0.5, tip.color=colt[filterED])
        axisPhylo(side=1)
    }
    
    pushViewport(viewport(layout=currentLayout))
    printMeta(paste("Dendrograms - Built from", ncol(scaledTransposedMatrix), 
                    "variables containing non-missing data", sep=" "), 
              pageno, currentjob)
}

#' Visualize number of DE variables for ANOVA and Kruskal Wallis
#' 
#' @param nr Normalyzer results object.
#' @param currentLayout Layout used for document.
#' @param pageno Current page number.
plotDEPlots <- function(nr, currentLayout, pageno) {

    nds <- nr@nds
    methodnames <- getUsedMethodNames(nr)
    methodlist <- getNormalizationMatrices(nr)
    currentjob <- nds@jobName
    
    ner <- nr@ner
    anfdr <- ner@anfdr
    kwfdr <- ner@kwfdr

    tout <- rbind(c(1, 2, 3), c(4))
    layout(tout)
    par(mar=c(2, 2, 2, 1), oma=c(2, 2, 3, 2), xpd=NA)
    
    barplot(colSums(anfdr < 0.05), main="ANOVA", names=c(methodnames), 
            border="red", density=20, cex=0.5, cex.axis=0.9, las=2,
            ylab="No. of Variables with FDR < 0.05")
    
    barplot(colSums(kwfdr < 0.05), main="Kruskal Wallis", names=c(methodnames), 
            border="red", density=20, cex=0.5, cex.axis=0.9, 
            las=2, ylab="No. of Variables with FDR < 0.05")
    
    pushViewport(viewport(layout=currentLayout))
    printMeta("Differential Expression", pageno, currentjob)
}

