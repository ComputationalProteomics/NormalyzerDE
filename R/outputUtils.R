#' Write normalization matrices to file
#' 
#' Outputs each of the normalized datasets to the specified directory.
#' 
#' @param nr Results object.
#' @param jobdir Path to output directory.
#' @param includePairwiseComparisons Include limma-based pairwise comparisons
#'   between all groups. For each normalization output file, columns named like
#'   \code{comp_A-B_p} and \code{comp_A-B_fdr} are added.
#' @param includeCvCol Include CV column in output.
#' @param includeAnovaP Include ANOVA p-value in output.
#' @param normSuffix String used to name output together with normalization names.
#' @param rawdataName Name of output raw data file.
#' @return None
#' @export
#' @examples
#' data(example_summarized_experiment)
#' normObj <- getVerifiedNormalyzerObject("job_name", example_summarized_experiment)
#' normResults <- normMethods(normObj)
#' normResultsWithEval <- analyzeNormalizations(normResults)
#' outputDir <- tempdir()
#' writeNormalizedDatasets(normResultsWithEval, outputDir)
writeNormalizedDatasets <- function(nr, jobdir, includePairwiseComparisons=FALSE, 
                                    includeCvCol=FALSE, includeAnovaP=FALSE,
                                    normSuffix="-normalized.txt",
                                    rawdataName="submitted_rawdata.txt") {
    
    nds <- nds(nr)
    ner <- ner(nr)

    methodlist <- normalizations(nr)
    methodnames <- names(methodlist)
    annotationColumns <- annotationValues(nds)
    if (ncol(annotationColumns) == 0) {
        annotationColumns <- NULL
    }

    if (includePairwiseComparisons) {
        designDf <- designMatrix(nds)
        groupCol <- groupNameCol(nds)
        if (!(groupCol %in% colnames(designDf))) {
            stop("Group column '", groupCol, "' not found in design matrix")
        }

        groupFactor <- as.factor(as.character(designDf[[groupCol]]))
        groupLevels <- levels(base::droplevels(groupFactor))

        if (length(groupLevels) < 2) {
            stop("At least two groups are required for pairwise comparisons")
        }

        safeGroupLevels <- make.names(groupLevels, unique=TRUE)
        groupMap <- stats::setNames(safeGroupLevels, groupLevels)

        design <- stats::model.matrix(~0 + groupFactor)
        colnames(design) <- safeGroupLevels

        groupPairs <- utils::combn(groupLevels, 2, simplify=FALSE)
        comparisonLabels <- vapply(groupPairs, function(pair) paste(pair, collapse="-"), "", USE.NAMES=FALSE)

        contrastMatrix <- matrix(
            0,
            nrow=length(groupLevels),
            ncol=length(groupPairs),
            dimnames=list(safeGroupLevels, comparisonLabels)
        )

        for (idx in seq_along(groupPairs)) {
            high <- groupPairs[[idx]][1]
            low <- groupPairs[[idx]][2]
            contrastMatrix[groupMap[[high]], idx] <- 1
            contrastMatrix[groupMap[[low]], idx] <- -1
        }
    }
    
    for (sampleIndex in seq_along(methodnames)) {
        
        currentMethod <- methodnames[sampleIndex]
        filePath <- paste(jobdir, "/", currentMethod, normSuffix, sep="")
        outputTable <- cbind(annotationColumns, methodlist[[sampleIndex]])

        if (includeAnovaP) {
            anovaP <- anovaP(ner)[,sampleIndex]
            
            if (nrow(outputTable) != length(anovaP)) {
                stop("Table row count: ", nrow(outputTable), 
                     " must match p-value vector length for anova: ", 
                     length(anovaP))
            }
            
            outputTable <- cbind(outputTable, anovaP=anovaP)
        }
        
        if (includePairwiseComparisons) {
            fit <- limma::lmFit(methodlist[[sampleIndex]], design)
            fit <- limma::contrasts.fit(fit, contrastMatrix)
            fit <- limma::eBayes(fit)

            pMat <- fit$p.value
            fdrMat <- apply(pMat, 2, function(p) stats::p.adjust(p, method="BH"))
            if (is.null(dim(fdrMat))) {
                fdrMat <- matrix(fdrMat, ncol=1)
            }
            colnames(fdrMat) <- colnames(pMat)

            compNames <- colnames(pMat)
            pColNames <- paste("comp", compNames, "p", sep="_")
            fdrColNames <- paste("comp", compNames, "fdr", sep="_")

            pairCols <- matrix(NA_real_, nrow=nrow(pMat), ncol=2 * length(compNames))
            pairCols[, seq(1, ncol(pairCols), by=2)] <- pMat
            pairCols[, seq(2, ncol(pairCols), by=2)] <- fdrMat
            colnames(pairCols) <- as.vector(rbind(pColNames, fdrColNames))

            outputTable <- cbind(outputTable, pairCols)
        }

        if (includeCvCol) {
            cvCol <- featureCVPerMethod(ner)[, sampleIndex]
            outputTable <- cbind(outputTable, CV=cvCol)
        }

        utils::write.table(
            outputTable, file=filePath, sep="\t", row.names=FALSE, quote=FALSE)
    }
    
    rawFilePath <- paste(jobdir, "/", rawdataName, sep="")
    rawOutputTable <- cbind(annotationColumns, filterrawdata(nds))
    
    utils::write.table(
        rawOutputTable, 
        file=rawFilePath, 
        sep="\t", 
        row.names=FALSE, 
        quote=FALSE
    )
}
