#' Generates a number of visualizations for the performance measures calculated
#' for the normalized matrices. These contain both general measures and
#' direct comparisons for different normalization approaches.
#' 
#' They include:
#' 
#' "Total intensity" 
#' Barplot showing the summed intensity in each sample for thelog2-transformed 
#' data
#' 
#' "Total missing" 
#' Barplot showing the number of missing values found in each sample for the 
#' log2-tranformed data
#' 
#' Log2-MDS plot: MDS plot where data is reduced to two dimensions allowing 
#' inspection of the main global changes in the data
#' 
#' PCV - Intragroup: Mean of intragroup CV of all replicate groups
#' 
#' PMAD - Intragroup: Mean of intragroup median absolute deviation across 
#' replicate groups
#' 
#' PEV - Intragroup: Mean of intragroup pooled estimate of variance across the 
#' replicate groups
#' 
#' Relative PCV, PMAD and PEV compared to log2: The results from PCV, PMAD
#' and PEV from all normalized data compared to the log2 data
#' 
#' Stable variables plot:
#' 5% of least differentially expressed variables are identified by ANOVA
#' analysis of log2 transformed data. Thereafter, global CV of these variables
#' is estimated from different normalized datasets. A plot of global CV of the
#' stable variables from all datsets on the y-axis and PCV-compared to log2 on
#' the x-axis is generated.
#' 
#' CV vs Raw Intensity plots:
#' For the first replicate group in each of the normalized dataset, a plot of
#' PCV of each variable compared to the average intensity of the variable in
#' the replicate group is plotted.
#' 
#' MA plots:
#' Plotted using the plotMA function of the limma package. The first sample in
#' each dataset is plotted against the average of the replicate group that
#' sample belong to.
#' 
#' Scatterplots:
#' The first two samples from each dataset are plotted.
#' 
#' Q-Q plots:
#' QQ-plots are plotted for the first sample in each normalized dataset.
#' 
#' Boxplots:
#' Boxplots for all samples are plotted and colored according to the replicate
#' grouping.
#' 
#' Relative Log Expression (RLE) plots:
#' Relative log expression value plots. Ratio between the expression of the
#' variable and the median expression of this variable across all samples.
#' The samples should be aligned around zero. Any deviation would indicate
#' discrepancies in the data.
#' 
#' Density plots:
#' Density distributions for each sample using the density function. Can 
#' capture outliers (if single densities lies far from the others) and see
#' if there is batch effects in the dataset (if for instance there is two
#' clear collections of lines in the data).
#' 
#' MDS plots
#' Multidimensional scaling plot using the cmdscale() function from the stats
#' package. Is often able to show whether replicates group together, and
#' whether there are any clear outliers in the data.
#' 
#' MeanSDplots
#' Displays the standard deviation values against values ordered according
#' to mean. If no dependency on mean is present (as is desired) a flat red
#' line is shown.
#' 
#' Pearson and Spearman correlation
#' Mean of intragroup Pearson and Spearman correlation values for each method.
#' 
#' Dendograms
#' Generated using the hclust function. Data is centered and scaled prior to
#' analysis. Coloring of replicates is done using as.phylo from the ape package.
#' 
#' P-value histograms
#' Histogram plots of p-values after calculating an ANOVA between different
#' condition groups. If no effect is present in the data a flat distribution
#' is expected. If an effect is present a flat distribution is still expected,
#' but with a sharp peak close to zero. If other effects are present it might
#' indicate that the data doesn't support the assumptions of ANOVA, for
#' instance if there are batch effects present in the data.
#' 
#' @param nr Normalyzer results object.
#' @param jobdir Path to output directory for run.
#' @param plotRows Number of plot rows.
#' @param plotCols Number of plot columns.
#' @return None
#' @export
#' @examples
#' data(example_summarized_experiment)
#' normObj <- getVerifiedNormalyzerObject("job_name", example_summarized_experiment)
#' normResults <- normMethods(normObj)
#' normResultsWithEval <- analyzeNormalizations(normResults)
#' outputDir <- tempdir()
#' generatePlots(normResultsWithEval, outputDir)
generatePlots <- function(nr, jobdir, plotRows=3, plotCols=4) {
    
    nds <- nds(nr)
    currentjob <- jobName(nds)
    nrows <- plotRows + 2
    ncols <- plotCols + 2
    
    currentLayout <- grid::grid.layout(
        nrow=nrows, 
        ncol=ncols,
        heights=c(0.1, rep(3 / (nrows - 2), (nrows-  2)), 0.1), 
        widths=c(0.1, rep(4 / (ncols - 2), (ncols - 2)), 0.1), 
        default.units=c('null', 'null')
    )
        
    currentFont <- "Helvetica"
    setupPlotting(currentjob, jobdir, "Norm-report")
    plotFrontPage(currentjob, currentFont)
    isLimitedRun <- singleReplicateRun(nds)
    
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
        plotPHist(nr, currentLayout, pageno)
    }
    
    grDevices::dev.off()
}

#' Setup PDF report settings by initializing the color palette, format
#' for the PDF report and the graphical device
#' 
#' @param currentJob Name of current run.
#' @param jobDir Path to output directory for run.
#' @param suffix Text to add to output filename.
#' @return None
#' @keywords internal
setupPlotting <- function(currentJob, jobDir, suffix) {
    
    grDevices::palette(c(
        "red", "green", "blue", "orange", "darkgray", "blueviolet", 
        "darkslateblue", "darkviolet", "gray", "bisque4", "brown", 
        "cadetblue4", "darkgreen", "darkcyan", "darkmagenta", "darkgoldenrod4", 
        "coral1"
    ))

    grDevices::pdf(
        file=paste(jobDir, "/", suffix, "-", currentJob, ".pdf", sep=""), 
        paper="a4r", width=0, height=0
    )

    themeNorm <- ggplot2::theme_set(ggplot2::theme_bw())
    themeNorm <- ggplot2::theme_update(
        panel.grid.minor=ggplot2::element_blank(), 
        axis.text=ggplot2::element_text(size=7), 
        axis.title=ggplot2::element_text(size=8), 
        plot.title=ggplot2::element_text(size=8), 
        plot.margin=ggplot2::unit(c(1, 1, 1, 1), "mm")
    )
    def.par <- graphics::par(no.readonly=TRUE)
}

#' Generate first page in output report and write to viewport
#' 
#' @param currentjob Name of current run.
#' @param currentFont Font used for output document.
#' @return None
#' @keywords internal
plotFrontPage <- function(currentjob, currentFont) {
    
    graphics::par(mfrow=c(4, 1))
    # TODO: Re-insert nice illustration (figure?)
    # data(data4pdftitle)
    graphics::plot(1, type="n", axes=FALSE, xlab="", ylab="")
    
    if ("NormalyzerDE" %in% rownames(utils::installed.packages())) {
      version <- utils::packageVersion("NormalyzerDE")
    }
    else {
      version <- "(version not found)"
    }
    
    la1 <- grid::grid.layout(
        nrow=7, 
        ncol=1, 
        heights=c(0.2, 1, 0.1, 0.2, 0.1, 0.2, 0.2), 
        default.units=c('null','null')
    )
    gpfill <- grid::gpar(fill="gray90", lwd=0, lty=0)
    grid::pushViewport(grid::viewport(layout=la1))
    grid::grid.rect(vp=grid::viewport(layout.pos.row=1), gp=gpfill)
    grid::grid.rect(vp=grid::viewport(layout.pos.row=7), gp=gpfill)
    
    grid::grid.text(
        paste("Project Name: ", currentjob, sep=""), 
        vp=grid::viewport(layout.pos.row=3), 
        just=c("center","center"), 
        gp=grid::gpar(fontsize=12, fontfamily=currentFont, col="black"))
    
    grid::grid.text(
        paste0("Normalyzer (ver ", version, " )"), 
        vp=grid::viewport(layout.pos.row=4), 
        just=c("center", "center"), 
        gp=grid::gpar(
            fontface="bold", 
            fontsize=32, 
            fontfamily=currentFont, 
            col="darkblue")
        )
    
    grid::grid.text(
        paste("Report created on: ", Sys.Date(), sep=""), 
        vp=grid::viewport(layout.pos.row=5), 
        just=c("center","center"),
        gp=grid::gpar(fontsize=12, fontfamily=currentFont, col="black"))
    
    citationText <- paste(
        "Citation: Willforss, J., Chawade, A. and Levander, F.",
        "Submitted"
    )
    
    grid::grid.text(
        citationText, 
        vp=grid::viewport(layout.pos.row=6), 
        just=c("center","center"), 
        gp=grid::gpar(fontsize=10, fontfamily=currentFont, col="black")
    )
    
    grid::grid.text(
        paste("Documentation for analyzing this report can be found at",
              "http://quantitativeproteomics.org/normalyzer/help.php"),
        vp=grid::viewport(layout.pos.row=7), 
        just=c("center", "center"), 
        gp=grid::gpar(fontsize=10, fontfamily=currentFont, col="black")
    )
}

#' Write page containing sample summary of intensities, missing values and 
#' MDS plot to the viewport
#' 
#' @param nr Normalyzer results object.
#' @param currentLayout Layout used for document.
#' @param pageno Current page number.
#' @return None
#' @keywords internal
plotSampleOutlierSummary <- function(nr, currentLayout, pageno) {
    
    nds <- nds(nr)
    methodlist <- normalizations(nr)
    filterED <- sampleReplicateGroups(nds)
    filterrawdata <- filterrawdata(nds)
    currentjob <- jobName(nds)
    
    tout <- rbind(c(1, 2), c(3, 4))
    graphics::layout(tout)
    graphics::par(mar=c(4, 4, 2, 1), oma=c(2, 2, 3, 2), xpd=NA)        

    datacoltotal <- colSums(filterrawdata, na.rm=TRUE)
    
    graphics::barplot(
        datacoltotal, 
        las=2, 
        main="Total intensity", 
        cex.names=0.5, 
        names.arg=substr(names(datacoltotal), 1, 10)
    )
    
    datamissingcol <- apply(filterrawdata, 2, function(x) { sum(is.na(x)) })
    graphics::barplot(
        datamissingcol, 
        las=2, 
        main="Total missing", 
        cex.names=0.5, 
        names.arg=substr(names(datamissingcol), 1, 10)
    )
    
    datastore <- methodlist[[1]]
    
    d <- stats::dist(
        scale(t(stats::na.omit(datastore)), center=TRUE, scale=TRUE)
    )
    fit <- stats::cmdscale(d, eig=TRUE, k=2)
    x <- fit$points[, 1]
    y <- fit$points[, 2]
    
    graphics::plot(x, y, type="n", main="Log2-MDS plot", xlab="", ylab="")
    graphics::text(
        fit$points[, 1], fit$points[, 2], col=filterED, labels=filterED
    )
    grid::pushViewport(grid::viewport(layout=currentLayout))
    
    printMeta(
        "Data Summary - Outlier detection", 
        pageno, 
        currentjob, 
        currentLayout
    )
}

#' Generate normalization replicate variance summary by displaying
#' CV (coefficient of variance), MAD (mean of intragroup median absolute
#' deviation) and PEV (Pooled Estimate of Variance) as mean of intragroups
#' 
#' @param nr Normalyzer results object.
#' @param currentLayout Layout used for document.
#' @param pageno Current page number.
#' @return None
#' @keywords internal
plotReplicateVariance <- function(nr, currentLayout, pageno) {
    
    nds <- nds(nr)
    methodnames <- names(normalizations(nr))
    currentjob <- jobName(nds)
    
    ner <- ner(nr)
    avgCVMem <- avgcvmem(ner)
    avgMADMem <- avgmadmem(ner)
    avgVarMem <- avgvarmem(ner)
    
    tout <- rbind(c(1, 2, 3), c(4))
    graphics::layout(tout)
    graphics::par(mar=c(2, 2, 2, 1), oma=c(2, 2, 3, 2), xpd=NA)
    
    graphics::boxplot(
        avgCVMem, 
        main="PCV - Intragroup", 
        names=c(methodnames), 
        border="red", 
        density=20, 
        cex=0.3, 
        cex.axis=0.9, 
        las=2, 
        frame.plot=FALSE
    )
    
    graphics::stripchart(
        as.data.frame(avgCVMem), 
        vertical=TRUE, 
        cex=0.4, 
        las=2, 
        pch=20, 
        add=TRUE, 
        col="darkgray"
    )
    
    graphics::boxplot(
        avgMADMem, 
        main="PMAD - Intragroup", 
        names=c(methodnames), 
        border="red", 
        density=20, 
        cex=0.3, 
        cex.axis=0.9, 
        las=2, 
        frame.plot=FALSE
    )
    
    graphics::stripchart(
        as.data.frame(avgMADMem), 
        vertical=TRUE, 
        cex=0.4, 
        las=2, 
        pch=20, 
        add=TRUE, 
        col="darkgray"
    )
    
    graphics::boxplot(
        avgVarMem, 
        main="PEV - Intragroup", 
        names=c(methodnames), 
        border="red", 
        density=20, 
        cex=0.3, 
        cex.axis=0.9, 
        las=2, 
        frame.plot=FALSE
    )
    
    graphics::stripchart(
        as.data.frame(avgVarMem), 
        vertical=TRUE, 
        cex=0.4, 
        las=2, 
        pch=20, 
        add=TRUE, 
        col="darkgray"
    )
    
    grid::pushViewport(grid::viewport(layout=currentLayout))
    printMeta("Replicate variation", pageno, currentjob, currentLayout)
}

#' Write figures displaying pooled coefficient of variance, median absolute 
#' deviation and pooled estimate of variance percentage compared to 
#' log2-transformed and stable variables plot displaying CV of stable variables
#' against pooled CV measure. The stable variables are calculated by an 
#' ANOVA comparison across sample conditions and selecting features with the 
#' least clear difference.
#'  
#' @param nr Normalyzer results object.
#' @param currentLayout Layout used for document.
#' @param pageno Current page number.
#' @return None
#' @keywords internal
plotReplicateVarAndStableVariables <- function(nr, currentLayout, pageno) {
    
    nds <- nds(nr)
    methodnames <- names(normalizations(nr))
    currentjob <- jobName(nds)
    
    ner <- ner(nr)
    lowVarFeaturesCVsPercDiff <- lowVarFeaturesCVsPercDiff(ner)
    avgcvmempdiff <- avgcvmempdiff(ner)
    avgmadmempdiff <- avgmadmempdiff(ner)
    avgvarmempdiff <- avgvarmempdiff(ner)
    
    tout <- rbind(c(1, 2, 3), c(4, 5, 5))
    graphics::layout(tout)
    graphics::par(mar=c(6, 6, 3, 1), oma=c(2, 3, 3, 2), xpd=NA)
    
    abc <- graphics::barplot(
        avgcvmempdiff, 
        main="PCV compared to log2 ", 
        names.arg=c(methodnames), 
        border="red", 
        ylim=c(min(avgcvmempdiff) - 10, 
               max(avgmadmempdiff) + 5), 
        density=20, 
        cex=0.9, 
        cex.axis=0.7, 
        las=2, 
        xpd=FALSE
    )
    
    graphics::axis(1, at=c(0.2, (max(abc) + 0.5)), labels=FALSE, lwd.ticks=0)
    graphics::axis(1, at=abc, labels=FALSE, lwd=0, lwd.ticks=1)
    graphics::text(abc,avgcvmempdiff, labels=round(avgcvmempdiff, digits=0), 
                   pos=3, las=2)
    
    abc <- graphics::barplot(
        avgmadmempdiff, 
        main="PMAD compared to log2", 
        names.arg=c(methodnames), 
        border="red", 
        ylim=c(min(avgmadmempdiff) - 10, 
               max(avgmadmempdiff) + 5), 
        density=20, 
        cex=0.9, 
        cex.axis=0.7, 
        las=2, 
        xpd=FALSE
    )
    
    graphics::axis(1, at=c(0.2, max(abc) + 0.5), labels=FALSE, lwd.ticks=0)
    graphics::axis(1, at=abc, labels=FALSE, lwd=0, lwd.ticks=1)
    graphics::text(
        abc, 
        avgmadmempdiff, 
        labels=round(avgmadmempdiff, digits=0), 
        pos=3, 
        las=2
    )
    
    abc <- graphics::barplot(
        avgvarmempdiff, 
        main="%PEV - compared to log2", 
        names.arg=c(methodnames), 
        border="red", 
        ylim=c(min(avgvarmempdiff) - 10, 
               max(avgmadmempdiff) + 5),
        density=20, 
        cex=0.9, 
        cex.axis=0.7, 
        las=2, 
        xpd=FALSE
    )
    
    graphics::axis(
        1, at = c(0.2, max(abc) + 0.5), labels = FALSE, lwd.ticks = 0
    )
    graphics::axis(
        1, at=abc, labels=FALSE, lwd=0, lwd.ticks=1
    )
    graphics::text(
        abc, 
        avgvarmempdiff, 
        labels=round(avgvarmempdiff, digits=0), 
        pos=3, 
        las=2
    )
    
    if (!all(is.na(lowVarFeaturesCVsPercDiff))) {
        
        if (min(avgcvmempdiff) < 0 || 
            max(avgcvmempdiff) > 100 || 
            min(lowVarFeaturesCVsPercDiff) < 0 || 
            max(lowVarFeaturesCVsPercDiff) > 100) {
            
            graphics::plot(
                avgcvmempdiff, lowVarFeaturesCVsPercDiff, pch=18, 
                xlim=c(0, 100), ylim=c(0, 100), 
                main="Stable variables plot", 
                xlab="PCV (intragroup) compared to Log2", 
                ylab="% Global CV of stable variables compared to Log2"
            )
            car::showLabels(
                avgcvmempdiff, 
                lowVarFeaturesCVsPercDiff, 
                labels=methodnames, 
                id.method="mahal", 
                id.cex=0.7, 
                id.col="black"
            )
        }
        else {
            graphics::plot(
                avgcvmempdiff, 
                lowVarFeaturesCVsPercDiff, 
                pch=18, 
                main="Stable variables plot", 
                xlab="PCV (Intragroup) compared to Log2", 
                ylab="% Global CV of stable variables compared to Log2"
            )
            car::showLabels(
                avgcvmempdiff, 
                lowVarFeaturesCVsPercDiff, 
                labels=methodnames, 
                id.method="mahal", 
                id.cex=0.7, 
                id.col="black"
            )
        }
    }
    
    grid::pushViewport(grid::viewport(layout=currentLayout))
    printMeta(
        "Replicate variation (Relative to Log2)", 
        pageno, 
        currentjob, 
        currentLayout
    )
}

#' Plots page displaying coefficient of variance (CV) against raw intensity
#' for features across the performed normalizations
#' 
#' @param nr Normalyzer results object.
#' @param currentLayout Layout used for document.
#' @param pageno Current page number.
#' @return None
#' @keywords internal
plotCVvsIntensity <- function(nr, currentLayout, pageno) {
    
    nds <- nds(nr)
    methodnames <- names(normalizations(nr))
    methodlist <- normalizations(nr)
    currentjob <- jobName(nds)
    sampleReplicateGroups <- sampleReplicateGroups(nds)
    filterrawdata <- filterrawdata(nds)
    
    log2Mat <- methodlist[[1]]
    tempcvmat1 <- matrix(
        nrow=nrow(log2Mat), 
        ncol=length(methodlist), 
        byrow=TRUE
    )
    tempavgmat1 <- matrix(
        nrow=nrow(log2Mat), 
        ncol=length(methodlist), 
        byrow=TRUE
    )
    maxtempcv <- 0
    
    for (j in seq_along(methodlist)) {
        
        log2Mat <- methodlist[[j]]
        log2Mat <- log2Mat[, sampleReplicateGroups == min(sampleReplicateGroups)]

        for (i in seq_len(nrow(log2Mat))) {
            
            tempcv <- RcmdrMisc::numSummary(
                log2Mat[i, ], 
                statistics=c("cv")
            )
            tempavg <- RcmdrMisc::numSummary(
                filterrawdata[i, ], 
                statistics=c("mean")
            )
            
            tempcvmat1[i, j] <- 100 * tempcv$table
            tempavgmat1[i, j] <- tempavg$table 
        }
        
        if (maxtempcv < max(tempcvmat1, na.rm=TRUE)) {
            maxtempcv <- max(tempcvmat1, na.rm=TRUE)
        }
    }
    
    tout <- matrix(
        seq_len((currentLayout$nrow - 2) * (currentLayout$ncol - 2)), 
        ncol=(currentLayout$ncol - 2), 
        byrow=TRUE
    )
    graphics::layout(tout)
    graphics::par(mar=c(4, 4, 2, 1), oma=c(2, 2, 3, 2), xpd=NA)
    
    for (i in seq_len(ncol(tempcvmat1))) {
        graphics::plot(
            tempavgmat1[, i], 
            tempcvmat1[, i], 
            main=methodnames[i], 
            xlab="Raw intensity", 
            ylab="CV", 
            cex=0.3, 
            ylim=c(0, maxtempcv)
        )
    }
    
    grid::pushViewport(grid::viewport(layout=currentLayout))
    printMeta("CV vs Raw Intensity plots", pageno, currentjob, currentLayout)
}

#' Produces a page containing expression vs. fold-change figures (MA plots)
#' The visualized fold is between the first sample in each group and the
#' average of the replicate to which that sample belongs
#' 
#' @param nr Normalyzer results object.
#' @param currentLayout Layout used for document.
#' @param pageno Current page number.
#'
#' @return None
#' @keywords internal
plotMA <- function(nr, currentLayout, pageno) {
    
    nds <- nds(nr)
    methodNames <- names(normalizations(nr))
    normalizedDataList <- normalizations(nr)
    currentjob <- jobName(nds)
    sampleReplicateGroups <- sampleReplicateGroups(nds)
    filterrawdata <- filterrawdata(nds)
    
    Malist <- list()
    for (i in seq_along(normalizedDataList)) {
        
        methodData <- as.matrix(normalizedDataList[[i]])
        methodDataFirstCond <- methodData[, sampleReplicateGroups == min(sampleReplicateGroups)]
        firstColWoNA <- methodDataFirstCond[!is.na(methodDataFirstCond[, 1]), ]
        avgExpr <- rowMeans(firstColWoNA)
        fold <- apply(cbind(firstColWoNA[, 1], avgExpr), 1, function(x) x[1] - x[2])
        plotDf <- as.data.frame(cbind(avgExpr, fold))
        
        Malist[[i]] <- ggplot2::ggplot(plotDf, ggplot2::aes(avgExpr, fold)) + 
            ggplot2::geom_point(color="darkgray", size=0.7, na.rm=TRUE) + 
            ggplot2::labs(x=("Replicate group mean"), 
                          y=("Replicate-1 Fold Change"),
                          title=methodNames[i]) + 
            ggplot2::stat_smooth(method="loess", se=FALSE, colour="red", na.rm=TRUE) + 
            ggplot2::geom_abline(intercept=0, slope=0, size=0.3)
    } 
    
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout=currentLayout))
    
    printPlots(Malist, "MA plots", pageno, currentjob, currentLayout)
}

#' Produces page containing scatter plot plotting the first two samples from 
#' each dataset against each other for each normalization method
#' 
#' @param nr Normalyzer results object.
#' @param currentLayout Layout used for document.
#' @param pageno Current page number.
#' @return None
#' @keywords internal
plotScatter <- function(nr, currentLayout, pageno) {
    
    nds <- nds(nr)
    methodnames <- names(normalizations(nr))
    methodlist <- normalizations(nr)
    currentjob <- jobName(nds)
    
    tout <- matrix(
        seq_len((currentLayout$nrow - 2) * (currentLayout$ncol - 2)), 
        ncol=(currentLayout$ncol - 2), 
        byrow=TRUE
    )
    graphics::layout(tout)
    graphics::par(mar=c(2, 2, 2, 1), oma=c(3, 2, 3, 2), xpd=NA)
    
    for (i in seq_along(methodlist)) {
        methodData <- methodlist[[i]]
        fit <- stats::lm(methodData[, 1]~methodData[, 2])
        graphics::plot(
            methodData[, 1], 
            methodData[, 2], 
            xlab="", 
            ylab="", 
            main=methodnames[i], 
            pch=19, 
            cex=0.2
        )
        
        graphics::legend(
            "topleft", 
            bty="n", 
            legend=paste("R2 ", format(summary(fit)$adj.r.squared, digits=2)), cex=0.7)
    }
    
    grid::pushViewport(grid::viewport(layout=currentLayout))
    printMeta("Scatterplots", pageno, currentjob, currentLayout)
}

#' Produces page showing QQ-plots for the first sample for each normalization
#' method. This plot can be used to assess whether the data follows a normal
#' distribution.
#' 
#' @param nr Normalyzer results object.
#' @param currentLayout Layout used for document.
#' @param pageno Current page number.
#' @return None
#' @keywords internal
plotQQ <- function(nr, currentLayout, pageno) {
    
    nds <- nds(nr)
    methodnames <- names(normalizations(nr))
    methodlist <- normalizations(nr)
    currentjob <- jobName(nds)
    
    qqlist <- list()
    
    for (i in seq_along(methodlist)) {  
        methodData <- methodlist[[i]]
        tempcolname <- colnames(methodData)
        qqlist[[i]] <- ggplot2::qplot(sample=methodData[, 1], na.rm=TRUE) + 
            ggplot2::labs(x="", y="", title=methodnames[i])
    }
    
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout=currentLayout))
    printPlots(qqlist, "Q-Q plots", pageno, currentjob, currentLayout)
}

#' Boxplots showing distribution of values after different normalizations
#' 
#' @param nr Normalyzer results object.
#' @param currentLayout Layout used for document.
#' @param pageno Current page number.
#' @return None
#' @keywords internal
plotBoxPlot <- function(nr, currentLayout, pageno) {
    
    nds <- nds(nr)
    methodnames <- names(normalizations(nr))
    methodlist <- normalizations(nr)
    currentjob <- jobName(nds)
    filterED <- sampleReplicateGroups(nds)
    filterrawdata <- filterrawdata(nds)
    
    tout <- matrix(
        seq_len((currentLayout$nrow - 2) * (currentLayout$ncol - 2)), 
        ncol=(currentLayout$ncol - 2), 
        byrow=TRUE
    )
    graphics::layout(tout)
    graphics::par(mar=c(2, 2, 2, 1), oma=c(3, 2, 3, 2), xpd=NA)
    mindata <- 1000
    maxdata <- 0
    
    for (i in seq_along(methodlist)) { 
        tempmin <- min(methodlist[[i]], na.rm=TRUE)
        tempmax <- max(methodlist[[i]], na.rm=TRUE)
        if (tempmin < mindata) {
            mindata <- tempmin
        }
        if (tempmax > maxdata) {
            maxdata <- tempmax
        }
    }
    
    for (i in seq_along(methodlist)) {   
        graphics::par(mar=c(5, 1, 1, 1))
        graphics::boxplot(
            methodlist[[i]], 
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
#' 
#' @param nr Normalyzer results object.
#' @param currentLayout Layout used for document.
#' @param pageno Current page number.
#' @return None
#' @keywords internal
plotRLE <- function(nr, currentLayout, pageno) {
    
    nds <- nds(nr)
    methodnames <- names(normalizations(nr))
    methodlist <- normalizations(nr)
    currentjob <- jobName(nds)
    filterED <- sampleReplicateGroups(nds)
    
    tout <- matrix(
        seq_len((currentLayout$nrow - 2) * (currentLayout$ncol - 2)), 
        ncol=(currentLayout$ncol - 2), 
        byrow=TRUE
    )
    graphics::layout(tout)
    graphics::par(mar=c(2, 2, 2, 1), oma=c(3, 2, 3, 2), xpd=NA)
    
    for (i in seq_along(methodlist)) {
        deviations = methodlist[[i]] - Biobase::rowMedians(methodlist[[i]], na.rm=TRUE)
        graphics::boxplot(
            deviations, 
            outcol="lightgray", 
            cex=0.1, 
            cex.axis=0.7, 
            las=2,
            main=methodnames[i], 
            col=(filterED), 
            names=substr(colnames(methodlist[[i]]), 
                         1, 
                         6))
    }
    
    grid::pushViewport(grid::viewport(layout=currentLayout))
    printMeta(
        "Relative Log Expression (RLE) plots", 
        pageno, 
        currentjob, 
        currentLayout
    )
}

#' Density plots showing value distributions after normalizations
#' 
#' @param nr Normalyzer results object.
#' @param currentLayout Layout used for document.
#' @param pageno Current page number.
#' @return None
#' @keywords internal
plotDensity <- function(nr, currentLayout, pageno) {
    
    nds <- nds(nr)
    methodnames <- names(normalizations(nr))
    methodlist <- normalizations(nr)
    currentjob <- jobName(nds)
    
    tout <- matrix(
        seq_len((currentLayout$nrow-2)*(currentLayout$ncol-2)), 
        ncol=(currentLayout$ncol-2), 
        byrow=TRUE
    )
    graphics::layout(tout)
    graphics::par(mar=c(3, 2, 3, 1), oma=c(3, 2, 3, 2), xpd=NA)
    
    for (i in seq_along(methodlist)) {
        
        methodData <- methodlist[[i]]
        tempd <- stats::density(methodData[, 1], na.rm=TRUE)
        graphics::plot(stats::density(methodData[, 1], na.rm=TRUE), xlab="", 
                       ylab="", ylim=c(min(tempd$y), max(tempd$y) * 1.5), 
                       main=methodnames[i], lty=2, lwd=1, col="darkgray")
        
        for (j in 2:ncol(methodData)) {
            
            graphics::lines(stats::density(methodData[, j], na.rm=TRUE), 
                            , lty=2, lwd=1, col="darkgray")
        }
    }
    
    grid::pushViewport(grid::viewport(layout=currentLayout))
    printMeta("Density plots", pageno, currentjob, currentLayout)
}

#' MDS plots showing grouping of samples after normalizations
#' 
#' @param nr Normalyzer results object.
#' @param currentLayout Layout used for document.
#' @param pageno Current page number.
#' @return None
#' @keywords internal
plotMDS <- function(nr, currentLayout, pageno) {
    
    nds <- nds(nr)
    methodnames <- names(normalizations(nr))
    methodlist <- normalizations(nr)
    currentjob <- jobName(nds)
    filterED <- sampleReplicateGroups(nds)
    
    tout <- matrix(
        seq_len((currentLayout$nrow-2)*(currentLayout$ncol-2)), 
        ncol=(currentLayout$ncol-2), 
        byrow=TRUE
    )
    graphics::layout(tout)
    graphics::par(mar=c(2, 2, 2, 1), oma=c(3, 2, 3, 2), xpd=NA)
    
    for (i in seq_along(methodlist)) {
        
        datastore <- methodlist[[i]]
        d <- stats::dist(
            scale(t(stats::na.omit(datastore)), 
                  center=TRUE, 
                  scale=TRUE)
            )
        fit <- stats::cmdscale(d, eig=TRUE, k=2)
        x <- fit$points[, 1]
        y <- fit$points[, 2]
        graphics::plot(x, y, type="n", main=methodnames[i], xlab="", ylab="")
        graphics::text(
            fit$points[, 1], 
            fit$points[, 2], 
            col=filterED, 
            labels=filterED
        )
    }
    
    grid::pushViewport(grid::viewport(layout=currentLayout))
    printMeta(
        paste("MDS plots - Built from", 
              ncol(d), 
              "variables with non-missing data", sep=" "), 
        pageno,
        currentjob, 
        currentLayout
    )
}

#' Visualize standard deviation over (expression?) for different values
#' 
#' @param nr Normalyzer results object.
#' @param currentLayout Layout used for document.
#' @param pageno Current page number.
#' @return None
#' @keywords internal
plotMeanSD <- function(nr, currentLayout, pageno) {
    
    nds <- nds(nr)
    methodnames <- names(normalizations(nr))
    methodlist <- normalizations(nr)
    currentjob <- jobName(nds)
    
    sdPlots <- list()
    
    for (i in seq_along(methodlist)) {
        
        methodData <- methodlist[[i]]
        msd <- vsn::meanSdPlot(
            methodData, 
            xlab="", 
            ylab="", 
            plot=FALSE, 
            na.rm=TRUE
        )
        
        sdPlots[[i]] <- msd$gg + ggplot2::ggtitle(methodnames[i]) +
            ggplot2::theme(legend.position="none", plot.margin=ggplot2::unit(c(1,0,0,0), "cm"))
    }
    
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout=currentLayout))
    printPlots(sdPlots, "MeanSDplots", pageno, currentjob, currentLayout)  
}

#' Visualize within-replicates correlations
#' 
#' @param nr Normalyzer results object.
#' @param currentLayout Layout used for document.
#' @param pageno Current page number.
#' @return None
#' @keywords internal
plotCorrelation <- function(nr, currentLayout, pageno) {
    
    ner <- ner(nr)
    repCorPear <- repCorPear(ner)
    repCorSpear <- repCorSpear(ner)
    
    nds <- nds(nr)
    methodnames <- names(normalizations(nr))
    methodlist <- normalizations(nr)
    currentjob <- jobName(nds)
    filterED <- sampleReplicateGroups(nds)
    filterrawdata <- filterrawdata(nds)
    
    tout <- rbind(c(1, 2), c(3))
    graphics::layout(tout)
    graphics::par(mar=c(2, 2, 2, 1), oma=c(2, 2, 3, 2), xpd=NA)
    
    # perdf <- data.frame(matrix(unlist(repCorPear), 
    #                            nrow=as.numeric(max(summary(repCorPear)[1])), 
    #                            byrow=TRUE))
    
    perdf <- data.frame(repCorPear)
    
    abc <- graphics::boxplot(perdf, main="Pearson correlation - Intragroup", 
                             names=c(methodnames), border="red", 
                             density=20, cex=0.3, cex.axis=0.9, las=2)
    
    graphics::stripchart(
        as.data.frame(perdf), 
        vertical=TRUE, 
        cex=0.4, 
        las=2, 
        pch=20, 
        add=TRUE, 
        col="darkgreen"
    )
    
    # spedf <- data.frame(matrix(unlist(repCorSpear), 
    #                            nrow=as.numeric(max(summary(repCorSpear)[1])), 
    #                            byrow=TRUE))
    
    spedf <- data.frame(repCorSpear)
    
    abc <- graphics::boxplot(spedf, main="Spearman correlation - Intragroup", 
                             names=c(methodnames), border="red", 
                             density=20, cex=0.3, cex.axis=0.9, las=2)
    
    graphics::stripchart(
        as.data.frame(spedf), 
        vertical=TRUE, 
        cex=0.4, 
        las=2, 
        pch=20, 
        add=TRUE, 
        col="darkgreen"
    )
    grid::pushViewport(grid::viewport(layout=currentLayout))
    printMeta("Correlation plots", pageno, currentjob, currentLayout)
}

#' Visualize dendrogram grouping of samples
#' 
#' @param nr Normalyzer results object.
#' @param currentLayout Layout used for document.
#' @param pageno Current page number.
#' @return None
#' @keywords internal
plotDendrograms <- function(nr, currentLayout, pageno) {
    
    nds <- nds(nr)
    methodnames <- names(normalizations(nr))
    methodlist <- normalizations(nr)
    currentjob <- jobName(nds)
    filterED <- sampleReplicateGroups(nds)
    
    tout <- matrix(
        seq_len((currentLayout$nrow - 2) * (currentLayout$ncol - 2)), 
        ncol=(currentLayout$ncol - 2), 
        byrow=TRUE
    )
    graphics::layout(tout)
    graphics::par(mar=c(2, 2, 2, 1), oma=c(2, 2, 3, 2), xpd=NA)
    
    colorVect <- c("red", "green", "blue", "orange", "darkgray", "blueviolet", 
                    "darkslateblue", "darkviolet", "gray", "bisque4", "brown", 
                    "cadetblue4", "darkgreen", "darkcyan", "darkmagenta", 
                    "darkgoldenrod4", "coral1")
    colt <- rep(colorVect, ceiling(length(filterED) / length(colorVect)))
    
    for (j in seq_along(methodlist)) {
        
        dataMatrix <- stats::na.omit(methodlist[[j]])
        colnames(dataMatrix) <- filterED
        scaledTransposedMatrix <- scale(t(dataMatrix), center=TRUE, scale=TRUE)

        hc <- stats::hclust(stats::dist(scaledTransposedMatrix), "ave")
        
        graphics::plot(
            ape::as.phylo(hc), 
            main=methodnames[j], 
            cex=0.5, 
            tip.color=colt[filterED]
        )
        ape::axisPhylo(side=1)
    }
    
    grid::pushViewport(grid::viewport(layout=currentLayout))
    printMeta(paste("Dendrograms - Built from", ncol(scaledTransposedMatrix), 
                    "variables containing non-missing data", sep=" "), 
              pageno, currentjob, currentLayout)
}

#' Generate P-histograms for ANOVA calculated after each normalization
#' 
#' @param nr Normalyzer results object.
#' @param currentLayout Layout used for document.
#' @param pageno Current page number.
#' @return None
#' @keywords internal
plotPHist <- function(nr, currentLayout, pageno) {
    
    nds <- nds(nr)
    ner <- ner(nr)
    methodnames <- names(normalizations(nr))
    currentjob <- jobName(nds)
    anovaP <- anovaP(ner)
    histPlots <- list()
    
    for (i in seq_along(methodnames)) {
        
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

