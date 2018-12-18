#' Remove technical replicates from data and design
#' 
#' Collapses sample values into their average. If only one value is present
#' due to NA-values in other technical replicates, then that value is used.
#' 
#' Takes a SummarizedExperiment where the data is present as the assay
#' and the colData contains the design conditions. In the design conditions
#' there should be one column with the technical replicate groups and one
#' column containing the sample names
#' 
#' @param se Summarized experiment where the assay contains the data to be
#'   reduced, and the colData the data frame
#' @param techRepColName Technical replicates column name in colData
#' @param sampleColName Sample names column name in colData
#' @return reducedSe Summarized experiment with reduced data
#' @export
#' @examples
#' testData <- as.matrix(data.frame(
#'     c(1,1,1), 
#'     c(1,2,1), 
#'     c(7,7,7), 
#'     c(7,9,7)))
#' colnames(testData) <- c("a1", "a2", "b1", "b2")
#' designDf <- data.frame(
#'     sample=c("a1", "a2", "b1", "b2"), 
#'     techrep=c("a", "a", "b", "b"))
#' se <- SummarizedExperiment::SummarizedExperiment(
#'     assay=testData,
#'     colData=designDf
#' )
#' statObj <- reduceTechnicalReplicates(se, "techrep", "sample")
reduceTechnicalReplicates <- function(se, techRepColName, sampleColName) {
    
    designDf <- data.frame(SummarizedExperiment::colData(se))
    dataMat <- SummarizedExperiment::assay(se)
        
    techRepGroups <- as.character(designDf[[techRepColName]])
    indices <- lapply(unique(techRepGroups), function(i) { which(techRepGroups %in% i) })
    
    # Reduce design
    collDesignMatList <- lapply(indices, function(inds) { designDf[inds[1],] })
    collDesignDf <- do.call(rbind.data.frame, collDesignMatList)
    collDesignDf <- droplevels(collDesignDf)
    rownames(collDesignDf) <- seq_len(nrow(collDesignDf))
    
    collapsedNames <- unlist(lapply(indices, function(inds) { paste(designDf[[sampleColName]][inds], collapse=".") }))
    collDesignDf[[sampleColName]] <- collapsedNames
    
    # Reduce data
    collDataMat <- as.matrix(data.frame(lapply(
        indices, 
        function(inds) {rowMeans(dataMat[, inds, drop=FALSE], na.rm = TRUE)})))
    colnames(collDataMat) <- collapsedNames
    
    reducedSe <- SummarizedExperiment::SummarizedExperiment(
        assay=collDataMat,
        colData=collDesignDf,
        rowData=SummarizedExperiment::rowData(se)
    )
    reducedSe
}

#' Generate an annotated data frame from statistics object
#' 
#' Extracts key values (p-value, adjusted p-value, log2-fold change and
#' average expression values) from an NormalyzerStatistics instance
#' and appends these to the annotation- and data-matrices
#' 
#' @param nst NormalyzerDE statistics object.
#' @return outDf Annotated statistics matrix
#' @export
#' @examples
#' data(example_stat_summarized_experiment)
#' statObj <- NormalyzerStatistics(example_stat_summarized_experiment)
#' statObj <- calculateContrasts(statObj, comparisons=c("1-2", "2-3"), condCol="group", type="limma")
#' annotDf <- generateAnnotatedMatrix(statObj)
generateAnnotatedMatrix <- function(nst) {
    
    pairwiseHead <- paste(
        names(pairwiseCompsP(nst)), "PValue", sep="_"
    )
    pMat <- data.frame(pairwiseCompsP(nst))
    colnames(pMat) <- pairwiseHead
    
    pairwiseHeadFdr <- paste(
        names(pairwiseCompsFdr(nst)), "AdjPVal", sep="_"
    )
    fdrMat <- data.frame(pairwiseCompsFdr(nst))
    colnames(fdrMat) <- pairwiseHeadFdr
    
    pairwiseHeadFold <- paste(
        names(pairwiseCompsFold(nst)), "log2FoldChange", sep="_"
    )
    foldMat <- data.frame(pairwiseCompsFold(nst))
    colnames(foldMat) <- pairwiseHeadFold
    
    aveMat <- data.frame(pairwiseCompsAve(nst))
    outDf <- data.frame(cbind(annotMat(nst), pMat, fdrMat, foldMat, 
                   featureAvg=aveMat[, 1], dataMat(nst)), check.names=FALSE)
    outDf
}

#' Generate full output report plot document. Plots p-value histograms
#' for each contrast in the NormalyzerStatistics instance and writes these
#' to a PDF report.
#' 
#' @param nst NormalyzerDE statistics object.
#' @param jobName Name of processing run.
#' @param jobDir Path to output directory.
#' @param sigThres Significance threshold for indicating as significant
#' @param sigThresType Type of significance threshold (FDR or p)
#' @param log2FoldThres log2 fold-change required for being counted as significant
#' @param plotRows Number of plot rows.
#' @param plotCols Number of plot columns.
#' @return None
#' @export
#' @examples
#' data(example_stat_summarized_experiment)
#' statObj <- NormalyzerStatistics(example_stat_summarized_experiment)
#' statObj <- calculateContrasts(statObj, comparisons=c("1-2", "2-3"), 
#'   condCol="group", type="limma")
#' outputDir <- tempdir()
#' generateStatsReport(statObj, "jobName", outputDir)
generateStatsReport <- function(nst, jobName, jobDir,
                                sigThres=0.1, sigThresType="fdr", log2FoldThres=0,
                                plotRows=3, plotCols=4) {
    
    nrows <- plotRows + 1
    ncols <- plotCols + 1
    
    currentLayout <- grid::grid.layout(
        nrow=nrows, 
        ncol=ncols,
        heights=c(0.1, rep(3 / (nrows - 2), (nrows - 2)), 0.1), 
        widths=c(0.1, rep(4 / (ncols - 2), (ncols - 2)), 0.1), 
        default.units=c("null", "null"))
    
    currentFont <- "Helvetica"
    setupPlotting(jobName, jobDir, "Norm-stats-report")
    plotFrontPage(jobName, currentFont)
    
    pageNo <- 2
    plotContrastPHists(nst, jobName, currentLayout, pageNo)
    
    pageNo <- pageNo + 1
    plotSigScatter(
        nst, jobName, currentLayout, pageNo, type="MA", 
        sigThres=sigThres, sigThresType=sigThresType, log2FoldThres=log2FoldThres
    )

    pageNo <- pageNo + 1
    plotSigScatter(
        nst, jobName, currentLayout, pageNo, type="Vulcano",
        sigThres=sigThres, sigThresType=sigThresType, log2FoldThres=log2FoldThres
    )
    
    pageNo <- pageNo + 1
    plotContrastPCA(
        nst, jobName, currentLayout, pageNo, pcs=c(1,2)
    )

    pageNo <- pageNo + 1
    plotContrastPCA(
        nst, jobName, currentLayout, pageNo, pcs=c(3,4)
    )
    
    if (length(comparisons(nst)) > 1) {
        pageNo <- pageNo + 1
        plotComparisonVenns(
            nst, jobName, currentLayout, pageNo,
            sigThres=sigThres, sigThresType=sigThresType, log2FoldThres=log2FoldThres
        )
    }
    
    grDevices::dev.off()
}

getSigs <- function(nst, sigThresType, sigThres, log2FoldThres) {

    contrastPLists <- pairwiseCompsP(nst)
    contrastFDRLists <- pairwiseCompsFdr(nst)
    contrastFoldLists <- pairwiseCompsFold(nst)

    sigLists <- list()
    
    getSigList <- function(contrast, contrastPLists, contrastFDRLists, contrastFoldLists) {
        
        pVals <- contrastPLists[[contrast]]
        fdrVals <- contrastFDRLists[[contrast]]
        fold <- contrastFoldLists[[contrast]]
        
        if (sigThresType == "fdr") {
            statSig <- (fdrVals < sigThres)
        }
        else if (sigThresType == "p") {
            statSig <- (pVals < sigThres)
        }
        else {
            stop("Unknown significance threshold type: ", sigThresType)
        }
        
        if (log2FoldThres != 0) {
            foldSig <- abs(fold) >= log2FoldThres
            sig <- statSig & foldSig
        }
        else {
            sig <- statSig
        }
        
        sig
    }

    sigLists <- lapply(names(contrastPLists), getSigList, 
           contrastPLists=contrastPLists,
           contrastFDRLists=contrastFDRLists,
           contrastFoldLists=contrastFoldLists
    )

    sigLists
}

#' If multiple comparisons - Show overlap in Venn diagrams
#' 
#' @param nst NormalyzerDE statistics object.
#' @param jobName Name of processing run.
#' @param currentLayout Layout used for document.
#' @param pageno Current page number.
#' @param sigThres Cutoff value for significance theshold
#' @param sigThresType Type of significance cutoff
#' @param log2FoldThres Log2-fold based cutoff threshold
#' @param maxContrasts Maximum contrasts to show pairwise comparisons for
#' @return None
#' @keywords internal
plotComparisonVenns <- function(nst, jobName, currentLayout, pageno,
                                sigThres=0.1, sigThresType="fdr", 
                                log2FoldThres=0, maxContrasts=4) {
    
    contrastSig <- getSigs(nst, sigThresType, sigThres, log2FoldThres)
    contrasts <- comparisons(nst)
    compCount <- min(length(contrasts), maxContrasts)
    folds <- pairwiseCompsFold(nst)
    plts <- list()
    index <- 0
    
    for (indOut in seq_len(compCount-1)) {
        for (indIn in (indOut+1):compCount) {

            # Credit goes to:
            # https://scriptsandstatistics.wordpress.com/2018/04/26/how-to-plot-venn-diagrams-using-r-ggplot2-and-ggforce/
            
            sigInds <- as.data.frame(contrastSig[c(indIn, indOut)])
            
            bothSigInds <- which(apply(sigInds, 1, all))
            contraCount <- length(which(sign(folds[[indIn]][bothSigInds]) != sign(folds[[indOut]][bothSigInds])))
            
            vdc <- limma::vennCounts(sigInds)
            class(vdc) <- 'matrix'
            df.vdc <- data.frame(vdc[-1, ])
            df.vdc <- rbind(df.vdc, c(-1, -1, contraCount))
            df.vdc$x <- c(-1.5, 1.5, 0, 0)
            df.vdc$y <- c(0, 0, 0, -0.3)
            df.vdc$label <- df.vdc$Counts
            df.vdc$label[3] <- paste0(df.vdc$label[3] - df.vdc$label[4], " same")
            df.vdc$label[4] <- paste0(df.vdc$label[4], " contra")
            
            df.venn <- data.frame(x = c(-0.5, 0.5),
                                  y = c(0, 0),
                                  labels = c(contrasts[indOut], contrasts[indIn]))
            
            plt <- ggplot2::ggplot(data=df.venn) +
                ggforce::geom_circle(
                    ggplot2::aes_string(x0 = "x", y0 = "y", r = 1.5, fill = "labels"), 
                    alpha = .3, size = 0.5, colour = 'darkgray') +
                ggplot2::coord_fixed() +
                ggplot2::theme_void() +
                ggplot2::theme(legend.position = 'bottom', legend.direction='vertical') +
                ggplot2::scale_fill_manual(values = c('cornflowerblue', 'firebrick',  'gold')) +
                ggplot2::scale_colour_manual(values = c('cornflowerblue', 'firebrick', 'gold'), guide = FALSE) +
                ggplot2::labs(fill = NULL) +
                ggplot2::annotate("text", x = df.vdc$x, y = df.vdc$y, label = df.vdc$label, size = 5) +
                ggplot2::ggtitle(paste0(sigThresType, " < ", sigThres, 
                           ", |log2 fold| >= ", log2FoldThres))
            
            index <- index + 1
            plts[[index]] <- plt
        }
    }
    
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout=currentLayout))
    printPlots(plts, "Venns", pageno, jobName, currentLayout)  
}


#' Takes an NormalyzerStatistics instance and generates and prints a p-value
#' histogram for each onto the viewport
#' 
#' @param nst NormalyzerDE statistics object.
#' @param jobName Name of processing run.
#' @param currentLayout Layout used for document.
#' @param pageno Current page number.
#' @return None
#' @keywords internal
plotContrastPHists <- function(nst, jobName, currentLayout, pageno) {

    contrastPLists <- pairwiseCompsP(nst)
    histPlots <- list()
    
    for (i in seq_along(contrastPLists)) {
        contrast <- names(contrastPLists)[i]
        pVals <- contrastPLists[[contrast]]
        df <- data.frame(pVals=pVals)
        histPlots[[i]] <- ggplot2::ggplot(df, ggplot2::aes(pVals)) + 
            ggplot2::theme_classic() + 
            ggplot2::geom_histogram(na.rm=TRUE, binwidth=0.01, color="#268BD2", fill="#268BD2") +
            ggplot2::ggtitle(paste("Contrast:", contrast)) + 
            ggplot2::xlab("P-values") + 
            ggplot2::ylab("Count") +
            ggplot2::theme(
                axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust=0.5), 
                legend.position="none"
            )
    }
    
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout=currentLayout))
    printPlots(histPlots, "HistPlots", pageno, jobName, currentLayout)  
}

#' Takes an NormalyzerStatistics instance and generates and prints a 
#' vulcano plot
#' 
#' @param nst NormalyzerDE statistics object.
#' @param jobName Name of processing run.
#' @param currentLayout Layout used for document.
#' @param pageno Current page number.
#' @param type Specify whether to plot 'Vulcano' or 'MA'.
#' @param sigThres FDR threshold for DE coloring.
#' @return None
#' @keywords internal
plotSigScatter <- function(nst, jobName, currentLayout, pageno, type="Vulcano",
                           sigThres=0.1, sigThresType="fdr", log2FoldThres=0) {
    
    contrastPLists <- pairwiseCompsP(nst)
    contrastFDRLists <- pairwiseCompsFdr(nst)
    contrastFoldLists <- pairwiseCompsFold(nst)
    contrastAveLists <- pairwiseCompsAve(nst)
    plots <- list()
    
    for (i in seq_along(contrastPLists)) {
        
        contrast <- names(contrastPLists)[i]
        pVals <- contrastPLists[[contrast]]
        fdrVals <- contrastFDRLists[[contrast]]
        fold <- contrastFoldLists[[contrast]]
        expression <- contrastAveLists[[contrast]]
        
        if (sigThresType == "fdr") {
            stat_sig <- (fdrVals < sigThres)
            legend_label <- paste("FDR <", sigThres)
        }
        else if (sigThresType == "p") {
            stat_sig <- (pVals < sigThres)
            legend_label <- paste("P <", sigThres)
        }
        else {
            stop("Unknown significance threshold type: ", sigThresType)
        }
        
        if (log2FoldThres != 0) {
            fold_sig <- abs(fold) >= log2FoldThres
            sig <- stat_sig & fold_sig
            legend_label <- paste0(legend_label, "\nFold > ", log2FoldThres)
        }
        else {
            sig <- stat_sig
        }
        
        df <- data.frame(pVals=pVals, fold=fold, sig=sig, expression=expression)
        df <- df[stats::complete.cases(df),]
        
        if (type == "Vulcano") {
            plots[[i]] <- ggplot2::ggplot(df, ggplot2::aes(fold, -log10(pVals), color=sig)) + 
                ggplot2::geom_point(size=1, alpha=0.5, na.rm=TRUE) +
                ggplot2::xlab("Fold (log2)") + 
                ggplot2::ylab("-log10(Pvalue)")
        }
        else if (type == "MA") {
            plots[[i]] <- ggplot2::ggplot(df, ggplot2::aes(expression, fold, color=sig)) + 
                ggplot2::geom_point(size=1, alpha=0.5, na.rm=TRUE) +
                ggplot2::xlab("Expression (log2)") + 
                ggplot2::ylab("Fold (log2)")
        }
        else {
            stop("Unknown plot type: ", type)
        }
        
        plots[[i]] <- plots[[i]] + 
            ggplot2::ggtitle(paste("Contrast:", contrast)) + 
            ggplot2::theme_classic() + 
            ggplot2::labs(color=legend_label)
    }
    
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout=currentLayout))
    printPlots(plots, type, pageno, jobName, currentLayout)  
}

#' Show in a PCA plot what samples are compared in statistical contrast
#' This is useful to understand what conditions are compared and for checking
#' for outliers in the contrast
#' 
#' @param nst NormalyzerDE statistics object.
#' @param jobName Name of processing run.
#' @param currentLayout Layout used for document.
#' @param pageno Current page number.
#' @param pcs Principal components to show.
#' @return None
#' @keywords internal
plotContrastPCA <- function(nst, jobName, currentLayout, pageno, pcs=c(1,2)) {
    
    activeLevelsFactor <- function(fullFactor, activeLevelsList) {
        
        activeIndsList <- list()
        for (levelName in names(activeLevelsList)) {
            activeLevel <- activeLevelsList[[levelName]]
            activeInds <- which(fullFactor %in% activeLevel)
            for (ind in activeInds) {
                activeIndsList[[as.character(ind)]] <- levelName
            }
        }
        getLevel <- function(index, activeIndsList) {
            if (index %in% names(activeIndsList)) {
                activeIndsList[[as.character(index)]]
            }
            else {
                "other"
            }
        }
        subLevel <- vapply(seq_len(length(fullFactor)), getLevel, 
                           "", activeIndsList = activeIndsList)
        subLevel
    }
    
    contrasts <- comparisons(nst)
    dataDf <- dataMat(nst)
    dfPCA <- stats::prcomp(t(dataDf[stats::complete.cases(dataDf), ]), scale=TRUE, center=TRUE)
    dfOut <- as.data.frame(dfPCA$x)
    percentageVar <- round(dfPCA$sdev^2 / sum(dfPCA$sdev^2) * 100, 2)
    percentageVar <- paste0(colnames(dfOut), " (", paste0(as.character(percentageVar), "%)"))
    groups <- condCol(nst)
    dfOut$sample <- colnames(dataDf)
    
    pc1 <- paste0("PC", pcs[1])
    pc2 <- paste0("PC", pcs[2])
    
    plots <- list()
    
    for (contrast in contrasts) {
        
        contrastLevels <- unlist(strsplit(contrast, "-"))
        
        dfOut$group <- activeLevelsFactor(
            condCol(nst), 
            list(high=contrastLevels[1], low=contrastLevels[2])
        )
        
        plt <- ggplot2::ggplot(
            dfOut, 
            ggplot2::aes_string(x=pc1, y=pc2, color="group", label="sample")) + 
            ggplot2::geom_text() + 
            ggplot2::theme_classic() +
            ggplot2::scale_color_manual(values=c("#00AAAA", "#AA0000", "#BBBBBB")) +
            ggplot2::ggtitle(paste(contrastLevels, collapse=" vs. ")) +
            ggplot2::xlab(percentageVar[pcs[1]]) +
            ggplot2::ylab(percentageVar[pcs[2]])  
        
        plots[[contrast]] <- plt
    }
    
    
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout=currentLayout))
    title <- paste0("PCA, contrast colored (factors", pcs[1], " and ", pcs[2], ")")
    printPlots(plots, title, pageno, jobName, currentLayout)  
}



