#' Generate full output report plot document
#' 
#' @param nr Normalyzer results object.
#' @param jobdir Path to output directory for run.
#' @param plotRows Number of plot rows.
#' @param plotCols Number of plot columns.
#' @return None
#' @export
#' @examples
#' normObj <- getVerifiedNormalyzerObject("data.tsv", "job_name", "design.tsv")
#' normResults <- normMethods(normObj)
#' normResultsWithEval <- analyzeNormalizations(normResults)
#' generatePlots(normResultsWithEval, "path/to/output")
generatePlots <- function(nr, jobdir, plotRows=3, plotCols=4) {
    
    nds <- nr@nds
    currentjob <- nds@jobName
    nrows <- plotRows + 2
    ncols <- plotCols + 2
    
    currentLayout <- grid::grid.layout(nrow=nrows, ncol=ncols,
                                       heights=c(0.1, rep(3/(nrows-2), (nrows-2)), 0.1), 
                                       widths=c(0.1, rep(4/(ncols-2), (ncols-2)), 0.1), 
                                       default.units=c('null', 'null'))
        
    currentFont <- "Helvetica"
    setupPlotting(currentjob, jobdir, "Norm-report")
    plotFrontPage(currentjob, currentFont)
    isLimitedRun <- nr@nds@singleReplicateRun
    
    # TI
    pageno <- 2
    plotSampleOutlierSummary(nr, currentLayout, pageno)
    
    # CV
    if (!isLimitedRun) {
        pageno <- pageno + 1
        plotReplicateVariance(nr, currentLayout, pageno)

        # Stable variables plot and CV in percent difference
        pageno <- pageno + 1
        plotReplicateVarAndStableVariables(nr, currentLayout, pageno)

        # CVvsintensityplot
        pageno <- pageno + 1
        plotCVvsIntensity(nr, currentLayout, pageno)
    }

    # MA plots
    pageno <- pageno + 1
    plotMA(nr, currentLayout, pageno)

    # Scatterplots
    pageno <- pageno + 1
    plotScatter(nr, currentLayout, pageno)

    # QQplot
    pageno <- pageno + 1
    plotQQ(nr, currentLayout, pageno)

    # Boxplot
    pageno <- pageno + 1
    plotBoxPlot(nr, currentLayout, pageno)

    # RLE plots
    pageno <- pageno + 1
    plotRLE(nr, currentLayout, pageno)

    # Density plots
    pageno <- pageno + 1
    plotDensity(nr, currentLayout, pageno)
    
    # MDS plot
    pageno <- pageno + 1
    plotMDS(nr, currentLayout, pageno)

    if (!isLimitedRun) {   
        
        # meanSDplot
        pageno <- pageno + 1
        plotMeanSD(nr, currentLayout, pageno)
        
        # Correlation
        pageno <- pageno + 1
        plotCorrelation(nr, currentLayout, pageno)
    }
     
    # Dendrograms
    pageno <- pageno + 1
    plotDendrograms(nr, currentLayout, pageno)

    if (!isLimitedRun) {
        # DE plots
        pageno <- pageno + 1
        # plotDEPlots(nr, currentLayout, pageno)
        plotPHist(nr, currentLayout, pageno)
    }
    
    grDevices::dev.off()
}

#' Setup PDF report settings
#' 
#' @param currentjob Name of current run.
#' @param jobdir Path to output directory for run.
#' @return None
setupPlotting <- function(currentjob, jobdir, suffix) {
    
    grDevices::palette(c("red", "green", "blue", "orange", "darkgray", 
                         "blueviolet", "darkslateblue", "darkviolet", "gray", 
                         "bisque4", "brown", "cadetblue4", "darkgreen", 
                         "darkcyan", "darkmagenta", "darkgoldenrod4", "coral1"))

    grDevices::pdf(file=paste(jobdir, "/", suffix, "-", currentjob, ".pdf", sep=""), 
                   paper="a4r", width=0, height=0)    

    themeNorm <- ggplot2::theme_set(ggplot2::theme_bw())
    themeNorm <- ggplot2::theme_update(panel.grid.minor=ggplot2::element_blank(), 
                                        axis.text=ggplot2::element_text(size=7), 
                                        axis.title=ggplot2::element_text(size=8), 
                                        plot.title=ggplot2::element_text(size=8), 
                                        plot.margin=ggplot2::unit(c(1, 1, 1, 1), "mm"))
    def.par <- graphics::par(no.readonly=TRUE)
}

#' Generate first page in output report
#' 
#' @param currentjob Name of current run.
#' @param currentFont Font used for output document.
#' @return None
plotFrontPage <- function(currentjob, currentFont) {
    
    graphics::par(mfrow=c(4, 1))
    # TODO: Re-insert nice illustration (figure?)
    # data(data4pdftitle)
    graphics::plot(1, type="n", axes=FALSE, xlab="", ylab="")
    
    if ("NormalyzerDE" %in% rownames(installed.packages())) {
      version <- packageVersion("NormalyzerDE")
    }
    else {
      version <- "(version not found)"
    }
    
    # TODO: Re-insert nice illustration (figure?)
    # boxplot(data4pdftitle, axes=FALSE, col=c("green","green","red","red"))
    la1 <- grid::grid.layout(nrow=7, ncol=1, heights=c(0.2, 1, 0.1, 0.2, 0.1, 0.2, 0.2), default.units=c('null','null'))
    gpfill <- grid::gpar(fill="gray90", lwd=0, lty=0)
    grid::pushViewport(grid::viewport(layout=la1))
    grid::grid.rect(vp=grid::viewport(layout.pos.row=1), gp=gpfill)
    grid::grid.rect(vp=grid::viewport(layout.pos.row=7), gp=gpfill)
    
    grid::grid.text(paste("Project Name: ", currentjob, sep=""), vp=grid::viewport(layout.pos.row=3), just=c("center","center"), 
                    gp=grid::gpar(fontsize=12, fontfamily=currentFont, col="black"))
    grid::grid.text(paste0("Normalyzer (ver ", version, " )"), vp=grid::viewport(layout.pos.row=4), just=c("center", "center"), 
                    gp=grid::gpar(fontface="bold", fontsize=32, fontfamily=currentFont, col="darkblue"))
    
    grid::grid.text(paste("Report created on: ", Sys.Date(), sep=""), vp=grid::viewport(layout.pos.row=5), just=c("center","center"),
                    gp=grid::gpar(fontsize=12, fontfamily=currentFont, col="black"))
    
    citationText <- "Citation: Chawade, A., Alexandersson, E., Levander, F. (2014). Normalyzer: a tool for rapid evaluation of normalization methods for omics data sets. J Proteome Res.,13 (6)"
    
    grid::grid.text(citationText, vp=grid::viewport(layout.pos.row=6), just=c("center","center"), gp=grid::gpar(fontsize=10, fontfamily=currentFont, col="black"))
    
    grid::grid.text("Documentation for analyzing this report can be found at http://quantitativeproteomics.org/normalyzer/help.php",
                    vp=grid::viewport(layout.pos.row=7), just=c("center", "center"), gp=grid::gpar(fontsize=10, fontfamily=currentFont, col="black"))
}

#' Generate sample summary of intensities, missing values and MDS plot
#' 
#' @param nr Normalyzer results object.
#' @param currentLayout Layout used for document.
#' @param pageno Current page number.
#' @return None
plotSampleOutlierSummary <- function(nr, currentLayout, pageno) {
    
    nds <- nr@nds
    methodlist <- getNormalizationMatrices(nr)
    filterED <- nds@sampleReplicateGroups
    filterrawdata <- nds@filterrawdata
    currentjob <- nds@jobName
    
    tout <- rbind(c(1, 2), c(3, 4))
    graphics::layout(tout)
    graphics::par(mar=c(4, 4, 2, 1), oma=c(2, 2, 3, 2), xpd=NA)        
    
    datacoltotal <- apply(filterrawdata, 2, function(x) { sum(x, na.rm=TRUE) })
    
    graphics::barplot(datacoltotal, las=2, main="Total intensity", cex.names=0.5, names.arg=substr(names(datacoltotal), 1, 10))
    datamissingcol <- apply(filterrawdata, 2, function(x) { sum(is.na(x)) })
    graphics::barplot(datamissingcol, las=2, main="Total missing", cex.names=0.5, names.arg=substr(names(datamissingcol), 1, 10))
    datastore <- methodlist[[1]]
    
    d <- stats::dist(scale(t(stats::na.omit(datastore)), center=TRUE, scale=TRUE))
    fit <- stats::cmdscale(d, eig=TRUE, k=2)
    x <- fit$points[, 1]
    y <- fit$points[, 2]
    
    graphics::plot(x, y, type="n", main="Log2-MDS plot", xlab="", ylab="")
    graphics::text(fit$points[, 1], fit$points[, 2], col=filterED, labels=filterED)
    grid::pushViewport(grid::viewport(layout=currentLayout))
    
    printMeta("Data Summary - Outlier detection", pageno, currentjob, currentLayout)
}

#' Generate normalization replicate variance summary
#' 
#' @param nr Normalyzer results object.
#' @param currentLayout Layout used for document.
#' @param pageno Current page number.
#' @return None
plotReplicateVariance <- function(nr, currentLayout, pageno) {
    
    nds <- nr@nds
    methodnames <- getUsedMethodNames(nr)
    currentjob <- nds@jobName
    
    ner <- nr@ner
    avgcvmem <- ner@avgcvmem
    avgmadmem <- ner@avgmadmem
    avgvarmem <- ner@avgvarmem
    
    tout <- rbind(c(1, 2, 3), c(4))
    graphics::layout(tout)
    graphics::par(mar=c(2, 2, 2, 1), oma=c(2, 2, 3, 2), xpd=NA)
    
    graphics::boxplot(avgcvmem, main="PCV - Intragroup", names=c(methodnames), 
                      border="red", density=20, cex=0.3, cex.axis=0.9, las=2, 
                      frame.plot=FALSE)
    
    graphics::stripchart(as.data.frame(avgcvmem), 
                         vertical=TRUE, cex=0.4, las=2, 
                         pch=20, add=TRUE, col="darkgray")
    
    graphics::boxplot(avgmadmem, main="PMAD - Intragroup", names=c(methodnames), 
                      border="red", density=20, cex=0.3, cex.axis=0.9, las=2, 
                      frame.plot=FALSE)
    
    graphics::stripchart(as.data.frame(avgmadmem), 
                         vertical=TRUE, cex=0.4, las=2, 
                         pch=20, add=TRUE, col="darkgray")
    
    graphics::boxplot(avgvarmem, main="PEV - Intragroup", names=c(methodnames), 
                      border="red", density=20, cex=0.3, cex.axis=0.9, las=2, 
                      frame.plot=FALSE)
    
    graphics::stripchart(as.data.frame(avgvarmem), 
                         vertical=TRUE, cex=0.4, las=2, 
                         pch=20, add=TRUE, col="darkgray")
    
    grid::pushViewport(grid::viewport(layout=currentLayout))
    printMeta("Replicate variation", pageno, currentjob, currentLayout)
}

#' Generate replicate variance overview for normalization methods and CV
#' plot for stable variables
#' 
#' @param nr Normalyzer results object.
#' @param currentLayout Layout used for document.
#' @param pageno Current page number.
#' @return None
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
    graphics::layout(tout)
    graphics::par(mar=c(6, 6, 3, 1), oma=c(2, 3, 3, 2), xpd=NA)
    
    abc <- graphics::barplot(avgcvmempdiff, main="PCV compared to log2 ", 
                             names.arg=c(methodnames), border="red", 
                             ylim=c(min(avgcvmempdiff) - 10, max(avgmadmempdiff) + 5), 
                             density=20, cex=0.9, cex.axis=0.7, las=2, xpd=FALSE)
    
    graphics::axis(1, at=c(0.2, (max(abc) + 0.5)), labels=FALSE, lwd.ticks=0)
    graphics::axis(1, at=abc, labels=FALSE, lwd=0, lwd.ticks=1)
    graphics::text(abc,avgcvmempdiff, labels=round(avgcvmempdiff, digits=0), 
                   pos=3, las=2)
    
    abc<-graphics::barplot(avgmadmempdiff, main="PMAD compared to log2", 
                           names.arg=c(methodnames), border="red", 
                           ylim=c(min(avgmadmempdiff) - 10, max(avgmadmempdiff) + 5), 
                           density=20, cex=0.9, cex.axis=0.7, las=2, xpd=FALSE)
    
    graphics::axis(1, at=c(0.2, max(abc) + 0.5), labels=FALSE, lwd.ticks=0)
    graphics::axis(1, at=abc, labels=FALSE, lwd=0, lwd.ticks=1)
    graphics::text(abc, avgmadmempdiff, 
                   labels=round(avgmadmempdiff, digits=0), 
                   pos=3, las=2)
    
    abc <- graphics::barplot(avgvarmempdiff, main="%PEV - compared to log2", 
                             names.arg=c(methodnames), border="red", 
                             ylim=c(min(avgvarmempdiff) - 10, max(avgmadmempdiff) + 5),
                             density=20, cex=0.9, cex.axis=0.7, las=2, xpd=FALSE)
    
    graphics::axis(1, at=c(0.2, max(abc) + 0.5), labels=FALSE, lwd.ticks=0)
    graphics::axis(1, at=abc, labels=FALSE, lwd=0, lwd.ticks=1)
    graphics::text(abc, avgvarmempdiff, labels=round(avgvarmempdiff, digits=0), pos=3, las=2)
    
    if (!all(is.na(nonsiganfdrlistcvpdiff))) {
        
        if (min(avgcvmempdiff) < 0 || max(avgcvmempdiff) > 100 || 
            min(nonsiganfdrlistcvpdiff) < 0 || max(nonsiganfdrlistcvpdiff) > 100) {
            
            graphics::plot(avgcvmempdiff, nonsiganfdrlistcvpdiff, pch=18, 
                           xlim=c(0, 100), ylim=c(0, 100), 
                           main="Stable variables plot", 
                           xlab="PCV (intragroup) compared to Log2", 
                           ylab="% Global CV of stable variables compared to Log2")
            car::showLabels(avgcvmempdiff, nonsiganfdrlistcvpdiff, 
                            labels=methodnames, id.method="mahal", 
                            id.cex=0.7, id.col="black")
        }
        else {
            graphics::plot(avgcvmempdiff, nonsiganfdrlistcvpdiff, pch=18, 
                           main="Stable variables plot", 
                           xlab="PCV (Intragroup) compared to Log2", 
                           ylab="% Global CV of stable variables compared to Log2")
            car::showLabels(avgcvmempdiff, nonsiganfdrlistcvpdiff, 
                            labels=methodnames, id.method="mahal", 
                            id.cex=0.7, id.col="black")
        }
    }
    
    grid::pushViewport(grid::viewport(layout=currentLayout))
    printMeta("Replicate variation (Relative to Log2)", pageno, currentjob, currentLayout)
}

#' Overview of coefficient of variation compared to intensity for variables
#'  for different normalization methods
#' 
#' @param nr Normalyzer results object.
#' @param currentLayout Layout used for document.
#' @param pageno Current page number.
#' @return None
plotCVvsIntensity <- function(nr, currentLayout, pageno) {
    
    nds <- nr@nds
    methodnames <- getUsedMethodNames(nr)
    methodlist <- getNormalizationMatrices(nr)
    currentjob <- nds@jobName
    sampleReplicateGroups <- nds@sampleReplicateGroups
    filterrawdata <- nds@filterrawdata
    
    datastore <- methodlist[[1]]
    tempcvmat1 <- matrix(nrow=nrow(datastore), 
                         ncol=length(methodlist), byrow=TRUE)
    tempavgmat1 <- matrix(nrow=nrow(datastore), 
                          ncol=length(methodlist), byrow=TRUE)
    maxtempcv <- 0
    
    for (j in 1:length(methodlist)) {
        
        datastore <- methodlist[[j]]
        datastore <- datastore[, sampleReplicateGroups == min(sampleReplicateGroups)]

        # browser()
        
        for (i in 1:nrow(datastore)) {
            
            # tempcv <- rcmdrNumSummary(datastore[i, ], statistics=c("cv"))
            # tempavg <- rcmdrNumSummary(filterrawdata[i, ], statistics=c("mean"))
            tempcv <- RcmdrMisc::numSummary(datastore[i, ], statistics=c("cv"))
            tempavg <- RcmdrMisc::numSummary(filterrawdata[i, ], statistics=c("mean"))
            
            tempcvmat1[i, j] <- 100 * tempcv$table
            tempavgmat1[i, j] <- tempavg$table 
        }
        
        if (maxtempcv < max(tempcvmat1, na.rm=TRUE)) {
            maxtempcv <- max(tempcvmat1, na.rm=TRUE)
        }
    }
    
    tout <- matrix(1:((currentLayout$nrow-2)*(currentLayout$ncol-2)), ncol=(currentLayout$ncol-2), byrow=TRUE)
    graphics::layout(tout)
    graphics::par(mar=c(4, 4, 2, 1), oma=c(2, 2, 3, 2), xpd=NA)
    
    for (i in 1:ncol(tempcvmat1)) {
        graphics::plot(tempavgmat1[, i], tempcvmat1[, i], main=methodnames[i], 
                       xlab="Raw intensity", ylab="CV", cex=0.3, 
                       ylim=c(0, maxtempcv))
    }
    
    grid::pushViewport(grid::viewport(layout=currentLayout))
    printMeta("CV vs Raw Intensity plots", pageno, currentjob, currentLayout)
}

#' Expression vs. fold-change for variables for each normalization method
#' ! TODO: Need to investigate what the actual fold change is here
#' 
#' @param nr Normalyzer results object.
#' @param currentLayout Layout used for document.
#' @param pageno Current page number.
#'
#' @return None
plotMA <- function(nr, currentLayout, pageno) {
    
    nds <- nr@nds
    methodnames <- getUsedMethodNames(nr)
    methodlist <- getNormalizationMatrices(nr)
    currentjob <- nds@jobName
    sampleReplicateGroups <- nds@sampleReplicateGroups
    filterrawdata <- nds@filterrawdata
    
    Malist <- list()
    for (i in 1:length(methodlist)) {
        
        datastore <- as.matrix(methodlist[[i]])
        tempcolname <- colnames(datastore)
        datastore <- datastore[, sampleReplicateGroups == min(sampleReplicateGroups)]
        datastore1 <- datastore[!is.na(datastore[, 1]), ]
        avg <- rowMeans(datastore1)
        fc <- apply(cbind(datastore1[, 1], avg), 1, function(x) x[1] - x[2])
        df <- as.data.frame(cbind(avg, fc))
        
        Malist[[i]] <- ggplot2::ggplot(df, ggplot2::aes(avg, fc)) + 
            ggplot2::geom_point(color="darkgray", size=0.7, na.rm=TRUE) + 
            ggplot2::labs(x=("Replicate group mean"), 
                          y=("Replicate-1 Fold Change"),
                          title=methodnames[i]) + 
            ggplot2::stat_smooth(method="loess", se=FALSE, colour="red", na.rm=TRUE) + 
            ggplot2::geom_abline(intercept=0, slope=0, size=0.3)
    } 
    
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout=currentLayout))
    
    printPlots(Malist, "MA plots", pageno, currentjob, currentLayout)
}

#' Scatter plot comparing replicate expression between samples
#' TODO: Investigate this further
#' 
#' @param nr Normalyzer results object.
#' @param currentLayout Layout used for document.
#' @param pageno Current page number.
#' @return None
plotScatter <- function(nr, currentLayout, pageno) {
    
    nds <- nr@nds
    methodnames <- getUsedMethodNames(nr)
    methodlist <- getNormalizationMatrices(nr)
    currentjob <- nds@jobName
    
    tout <- matrix(1:((currentLayout$nrow-2)*(currentLayout$ncol-2)), ncol=(currentLayout$ncol-2), byrow=TRUE)
    graphics::layout(tout)
    graphics::par(mar=c(2, 2, 2, 1), oma=c(3, 2, 3, 2), xpd=NA)
    
    for (i in 1:length(methodlist)) {
        datastore <- methodlist[[i]]
        fit <- stats::lm(datastore[, 1]~datastore[, 2])
        graphics::plot(datastore[, 1], datastore[, 2], xlab="", ylab="", main=methodnames[i], pch=19, cex=0.2)
        
        graphics::legend("topleft", bty="n", legend=paste("R2 ", format(summary(fit)$adj.r.squared, digits=2)), cex=0.7)
    }
    
    grid::pushViewport(grid::viewport(layout=currentLayout))
    printMeta("Scatterplots", pageno, currentjob, currentLayout)
}

#' QQ-plots for variable values
#' TODO: Investigate this further - What is actually extracted here?
#' 
#' @param nr Normalyzer results object.
#' @param currentLayout Layout used for document.
#' @param pageno Current page number.
#' @return None
plotQQ <- function(nr, currentLayout, pageno) {
    
    nds <- nr@nds
    methodnames <- getUsedMethodNames(nr)
    methodlist <- getNormalizationMatrices(nr)
    currentjob <- nds@jobName
    
    qqlist <- list()
    
    for (i in 1:length(methodlist)) {  
        datastore <- methodlist[[i]]
        tempcolname <- colnames(datastore)
        qqlist[[i]] <- ggplot2::qplot(sample=datastore[, 1], na.rm=TRUE) + 
            ggplot2::labs(x="", y="", title=methodnames[i])
    }
    
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout=currentLayout))
    printPlots(qqlist, "Q-Q plots", pageno, currentjob, currentLayout)
}

#' Boxplots showing distribution of values after different normalizations
#' TODO: Investigate this further - What is actually extracted here?
#' 
#' @param nr Normalyzer results object.
#' @param currentLayout Layout used for document.
#' @param pageno Current page number.
#' @return None
plotBoxPlot <- function(nr, currentLayout, pageno) {
    
    nds <- nr@nds
    methodnames <- getUsedMethodNames(nr)
    methodlist <- getNormalizationMatrices(nr)
    currentjob <- nds@jobName
    filterED <- nds@sampleReplicateGroups
    filterrawdata <- nds@filterrawdata
    
    tout <- matrix(1:((currentLayout$nrow-2)*(currentLayout$ncol-2)), ncol=(currentLayout$ncol-2), byrow=TRUE)
    graphics::layout(tout)
    graphics::par(mar=c(2, 2, 2, 1), oma=c(3, 2, 3, 2), xpd=NA)
    mindata <- 1000
    maxdata <- 0
    
    for (i in 1:length(methodlist)) { 
        tempmin <- min(methodlist[[i]], na.rm=TRUE)
        tempmax <- max(methodlist[[i]], na.rm=TRUE)
        if (tempmin < mindata) {
            mindata <- tempmin
        }
        if (tempmax > maxdata) {
            maxdata <- tempmax
        }
    }
    
    for (i in 1:length(methodlist)) {   
        graphics::par(mar=c(5, 1, 1, 1))
        graphics::boxplot(methodlist[[i]], 
                          cex=0.1, 
                          cex.axis=0.7, 
                          las=2, 
                          main=methodnames[i], 
                          col=(filterED), 
                          outcol="lightgray", 
                          ylim=c(mindata - 1, maxdata + 1), 
                          names=substr(colnames(methodlist[[i]]), 1, 10))
    }
    
    grid::pushViewport(grid::viewport(layout=currentLayout))
    printMeta("Boxplots", pageno, currentjob, currentLayout)
}

#' Boxplots showing relative log expression after normalizations
#' TODO: Investigate this further - What is actually extracted here?
#' 
#' @param nr Normalyzer results object.
#' @param currentLayout Layout used for document.
#' @param pageno Current page number.
#' @return None
plotRLE <- function(nr, currentLayout, pageno) {
    
    nds <- nr@nds
    methodnames <- getUsedMethodNames(nr)
    methodlist <- getNormalizationMatrices(nr)
    currentjob <- nds@jobName
    filterED <- nds@sampleReplicateGroups
    
    tout <- matrix(1:((currentLayout$nrow-2)*(currentLayout$ncol-2)), ncol=(currentLayout$ncol-2), byrow=TRUE)
    graphics::layout(tout)
    graphics::par(mar=c(2, 2, 2, 1), oma=c(3, 2, 3, 2), xpd=NA)
    
    for (i in 1:length(methodlist)) {
        deviations = methodlist[[i]] - Biobase::rowMedians(methodlist[[i]], na.rm=TRUE)
        graphics::boxplot(deviations, outcol="lightgray", cex=0.1, cex.axis=0.7, 
                          las=2, main=methodnames[i], col=(filterED), 
                          names=substr(colnames(methodlist[[i]]), 1, 6))
    }
    
    grid::pushViewport(grid::viewport(layout=currentLayout))
    printMeta("Relative Log Expression (RLE) plots", pageno, currentjob, currentLayout)
}

#' Density plots showing value distributions after normalizations
#' 
#' @param nr Normalyzer results object.
#' @param currentLayout Layout used for document.
#' @param pageno Current page number.
#' @return None
plotDensity <- function(nr, currentLayout, pageno) {
    
    nds <- nr@nds
    methodnames <- getUsedMethodNames(nr)
    methodlist <- getNormalizationMatrices(nr)
    currentjob <- nds@jobName
    
    tout <- matrix(1:((currentLayout$nrow-2)*(currentLayout$ncol-2)), ncol=(currentLayout$ncol-2), byrow=TRUE)
    graphics::layout(tout)
    graphics::par(mar=c(3, 2, 3, 1), oma=c(3, 2, 3, 2), xpd=NA)
    
    for (i in 1:length(methodlist)) {
        
        datastore <- methodlist[[i]]
        tempd <- stats::density(datastore[, 1], na.rm=TRUE)
        graphics::plot(stats::density(datastore[, 1], na.rm=TRUE), xlab="", 
                       ylab="", ylim=c(min(tempd$y), max(tempd$y) * 1.5), 
                       main=methodnames[i], lty=2, lwd=1, col="darkgray")
        
        for (j in 2:ncol(datastore)) {
            
            graphics::lines(stats::density(datastore[, j], na.rm=TRUE), 
                            , lty=2, lwd=1, col="darkgray")
        }
    }
    
    grid::pushViewport(grid::viewport(layout=currentLayout))
    printMeta("Density plots", pageno, currentjob, currentLayout)
}

#' MDS plots showing grouping of samples after normalizations
#' TODO: Investigate this further - What is actually extracted here?
#' TODO: Do the d-value with number non-missing data vanish?
#' 
#' @param nr Normalyzer results object.
#' @param currentLayout Layout used for document.
#' @param pageno Current page number.
#' @return None
plotMDS <- function(nr, currentLayout, pageno) {
    
    nds <- nr@nds
    methodnames <- getUsedMethodNames(nr)
    methodlist <- getNormalizationMatrices(nr)
    currentjob <- nds@jobName
    filterED <- nds@sampleReplicateGroups
    
    tout <- matrix(1:((currentLayout$nrow-2)*(currentLayout$ncol-2)), ncol=(currentLayout$ncol-2), byrow=TRUE)
    graphics::layout(tout)
    graphics::par(mar=c(2, 2, 2, 1), oma=c(3, 2, 3, 2), xpd=NA)
    
    for (i in 1:length(methodlist)) {
        
        datastore <- methodlist[[i]]
        d <- stats::dist(scale(t(stats::na.omit(datastore)), center=TRUE, scale=TRUE))
        fit <- stats::cmdscale(d, eig=TRUE, k=2)
        x <- fit$points[, 1]
        y <- fit$points[, 2]
        graphics::plot(x, y, type="n", main=methodnames[i], xlab="", ylab="")
        graphics::text(fit$points[, 1], fit$points[, 2], col=filterED, labels=filterED)
    }
    
    grid::pushViewport(grid::viewport(layout=currentLayout))
    printMeta(paste("MDS plots - Built from", ncol(d), "variables with non-missing data", sep=" "), pageno, currentjob, currentLayout)
}

#' Visualize standard deviation over (expression?) for different values
#' 
#' @param nr Normalyzer results object.
#' @param currentLayout Layout used for document.
#' @param pageno Current page number.
#' @return None
plotMeanSD <- function(nr, currentLayout, pageno) {
    
    nds <- nr@nds
    methodnames <- getUsedMethodNames(nr)
    methodlist <- getNormalizationMatrices(nr)
    currentjob <- nds@jobName
    
    sdPlots <- list()
    
    for (i in 1:length(methodlist)) {
        
        datastore <- methodlist[[i]]
        msd <- vsn::meanSdPlot(datastore, xlab="", ylab="", plot=FALSE, na.rm=TRUE)
        
        sdPlots[[i]] <- msd$gg + ggplot2::ggtitle(methodnames[i]) +
            ggplot2::theme(legend.position="none", plot.margin=ggplot2::unit(c(1,0,0,0), "cm"))
    }
    
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout=currentLayout))
    printPlots(sdPlots, "MeanSDplots", pageno, currentjob, currentLayout)  
}

#' Visualize correlations for plots
#' TODO: What exactly is going on here?
#' 
#' @param nr Normalyzer results object.
#' @param currentLayout Layout used for document.
#' @param pageno Current page number.
#' @return None
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
    graphics::layout(tout)
    graphics::par(mar=c(2, 2, 2, 1), oma=c(2, 2, 3, 2), xpd=NA)
    
    perdf <- data.frame(matrix(unlist(avgpercorsum), 
                               nrow=as.numeric(max(summary(avgpercorsum)[1])), 
                               byrow=TRUE))
    
    abc <- graphics::boxplot(perdf, main="Pearson correlation - Intragroup", 
                             names=c(methodnames), border="red", 
                             density=20, cex=0.3, cex.axis=0.9, las=2)
    
    graphics::stripchart(as.data.frame(perdf), vertical=TRUE, cex=0.4, las=2, pch=20, add=TRUE, col="darkgreen")
    
    spedf <- data.frame(matrix(unlist(avgspecorsum), 
                               nrow=as.numeric(max(summary(avgspecorsum)[1])), 
                               byrow=TRUE))
    
    abc <- graphics::boxplot(spedf, main="Spearman correlation - Intragroup", 
                             names=c(methodnames), border="red", 
                             density=20, cex=0.3, cex.axis=0.9, las=2)
    
    graphics::stripchart(as.data.frame(spedf), 
                         vertical=TRUE, cex=0.4, las=2, pch=20, add=TRUE, col="darkgreen")
    grid::pushViewport(grid::viewport(layout=currentLayout))
    printMeta("Correlation plots", pageno, currentjob, currentLayout)
}

#' Visualize dendrogram grouping of samples
#' TODO: Why is the Quantile one crazy? (Numbers instead of replicate names)
#' 
#' @param nr Normalyzer results object.
#' @param currentLayout Layout used for document.
#' @param pageno Current page number.
#' @return None
plotDendrograms <- function(nr, currentLayout, pageno) {
    
    nds <- nr@nds
    methodnames <- getUsedMethodNames(nr)
    methodlist <- getNormalizationMatrices(nr)
    currentjob <- nds@jobName
    filterED <- nds@sampleReplicateGroups
    
    tout <- matrix(1:((currentLayout$nrow-2)*(currentLayout$ncol-2)), ncol=(currentLayout$ncol-2), byrow=TRUE)
    graphics::layout(tout)
    graphics::par(mar=c(2, 2, 2, 1), oma=c(2, 2, 3, 2), xpd=NA)
    
    colorVect <- c("red", "green", "blue", "orange", "darkgray", "blueviolet", 
                    "darkslateblue", "darkviolet", "gray", "bisque4", "brown", 
                    "cadetblue4", "darkgreen", "darkcyan", "darkmagenta", 
                    "darkgoldenrod4", "coral1")
    colt <- rep(colorVect, ceiling(length(filterED) / length(colorVect)))
    
    for (j in 1:length(methodlist)) {
        
        dataMatrix <- stats::na.omit(methodlist[[j]])
        colnames(dataMatrix) <- filterED
        scaledTransposedMatrix <- scale(t(dataMatrix), center=TRUE, scale=TRUE)

        hc <- stats::hclust(stats::dist(scaledTransposedMatrix), "ave")
        
        graphics::plot(ape::as.phylo(hc), main=methodnames[j], cex=0.5, tip.color=colt[filterED])
        ape::axisPhylo(side=1)
    }
    
    grid::pushViewport(grid::viewport(layout=currentLayout))
    printMeta(paste("Dendrograms - Built from", ncol(scaledTransposedMatrix), 
                    "variables containing non-missing data", sep=" "), 
              pageno, currentjob, currentLayout)
}

#' Visualize number of DE variables for ANOVA and Kruskal Wallis
#' 
#' @param nr Normalyzer results object.
#' @param currentLayout Layout used for document.
#' @param pageno Current page number.
#' @return None
plotDEPlots <- function(nr, currentLayout, pageno) {
    
    fdrThreshold <- 0.05
    
    nds <- nr@nds
    methodnames <- getUsedMethodNames(nr)
    methodlist <- getNormalizationMatrices(nr)
    currentjob <- nds@jobName
    
    ner <- nr@ner
    anfdr <- ner@anfdr
    kwfdr <- ner@kwfdr
    
    tout <- rbind(c(1, 2, 3), c(4))
    graphics::layout(tout)
    graphics::par(mar=c(2, 2, 2, 1), oma=c(2, 2, 3, 2), xpd=NA)
    
    graphics::barplot(sapply(anfdr, function(col) { length(col[col < fdrThreshold]) }), 
                      main="ANOVA", names=c(methodnames), 
                      border="red", density=20, cex=0.5, cex.axis=0.9, las=2,
                      ylab=paste("No. of Variables with FDR <", fdrThreshold))
    
    graphics::barplot(sapply(kwfdr, function(col) { length(col[col < fdrThreshold]) }), 
                      main="Kruskal Wallis", names=c(methodnames), 
                      border="red", density=20, cex=0.5, cex.axis=0.9, las=2, 
                      ylab=paste("No. of Variables with FDR <", fdrThreshold
                             ))
    
    grid::pushViewport(grid::viewport(layout=currentLayout))
    printMeta("Differential Expression", pageno, currentjob, currentLayout)
}

#' Generate P-histograms for ANOVA calculated after each normalization
#' 
#' @param nr Normalyzer results object.
#' @param currentLayout Layout used for document.
#' @param pageno Current page number.
#' @return None
plotPHist <- function(nr, currentLayout, pageno) {
    
    nds <- nr@nds
    ner <- nr@ner
    methodnames <- getUsedMethodNames(nr)
    currentjob <- nds@jobName
    anovaP <- ner@anovaP
    histPlots <- list()
    
    for (i in 1:length(methodnames)) {
        
        anovaPVals <- anovaP[, i]
        df <- data.frame(anovaPVals=anovaPVals)
        histPlots[[i]] <- ggplot2::ggplot(df) + 
            ggplot2::geom_histogram(ggplot2::aes(anovaPVals), na.rm=TRUE, binwidth=0.01) +
            ggplot2::ggtitle(methodnames[i]) +
            ggplot2::xlim(0,1)
    }
    
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout=currentLayout))
    printPlots(histPlots, "HistPlots", pageno, currentjob, currentLayout)  
}

