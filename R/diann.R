diannIsParquet <- function(filePath) {
    grepl("\\.parquet$", filePath, ignore.case=TRUE)
}

diannReadHeader <- function(filePath, sep="\t") {
    headerLine <- readLines(filePath, n=1, warn=FALSE)
    if (length(headerLine) < 1) {
        stop("DIANN input file was empty: ", filePath)
    }
    strsplit(headerLine, sep, fixed=TRUE)[[1]]
}

diannIsReportHeader <- function(header) {
    any(c("Run", "File.Name") %in% header)
}

diannSelectFirstPresent <- function(candidates, available) {
    first <- candidates[candidates %in% available][1]
    if (is.na(first) || is.null(first)) {
        NULL
    }
    else {
        first
    }
}

diannResolveMinPositive <- function(diannMinPositive) {
    if (is.null(diannMinPositive)) {
        return(0)
    }

    diannMinPositive <- as.numeric(diannMinPositive)[1]
    if (is.na(diannMinPositive) || diannMinPositive < 0) {
        stop("diannMinPositive must be a single non-negative numeric value.")
    }

    diannMinPositive
}

diannStripPath <- function(paths) {
    sub("^.*[\\\\/]", "", paths)
}

diannStripKnownExtension <- function(fileNames) {
    sub("\\.(raw|mzml|d|wiff)$", "", fileNames, ignore.case=TRUE)
}

diannCleanSampleName <- function(sampleNames) {
    diannStripKnownExtension(diannStripPath(sampleNames))
}

diannCandidateSampleColumns <- function(columnNames) {
    grepl("[\\\\/]|\\.(raw|mzml|d|wiff)$", columnNames, ignore.case=TRUE)
}

diannRenameSampleColumnsForDesign <- function(dataFrame, designSampleNames) {
    if (is.null(designSampleNames) || length(designSampleNames) == 0) {
        return(dataFrame)
    }

    designSampleNames <- as.character(designSampleNames)
    if (all(designSampleNames %in% colnames(dataFrame))) {
        return(dataFrame)
    }

    candidateCols <- diannCandidateSampleColumns(colnames(dataFrame))
    if (!any(candidateCols)) {
        return(dataFrame)
    }

    cleaned <- colnames(dataFrame)
    cleaned[candidateCols] <- diannCleanSampleName(cleaned[candidateCols])

    dupNames <- unique(cleaned[candidateCols][duplicated(cleaned[candidateCols])])
    if (length(dupNames) > 0) {
        stop(
            "DIANN sample columns are not unique after stripping paths/extensions. ",
            "Duplicate sample names include: ", paste(utils::head(dupNames, 10), collapse=", "), "\n",
            "Provide unique sample names in DIA-NN export, or use full file paths in the design matrix."
        )
    }

    if (all(designSampleNames %in% cleaned)) {
        colnames(dataFrame) <- cleaned
    }

    dataFrame
}

diannChooseReportSpec <- function(reportColumns,
                                 diannLevel=c("auto", "protein", "precursor"),
                                 diannSampleCol=NULL,
                                 diannFeatureCol=NULL,
                                 diannQuantityCol=NULL) {

    diannLevel <- match.arg(diannLevel)

    sampleCol <- if (!is.null(diannSampleCol)) {
        if (!(diannSampleCol %in% reportColumns)) {
            stop("DIANN report file is missing requested sample column: ", diannSampleCol)
        }
        diannSampleCol
    }
    else if ("Run" %in% reportColumns) {
        "Run"
    }
    else if ("File.Name" %in% reportColumns) {
        "File.Name"
    }
    else {
        stop("DIANN report file is missing both 'Run' and 'File.Name' columns")
    }

    proteinFeatureCol <- "Protein.Group"
    precursorFeatureCol <- "Precursor.Id"
    proteinQuantityCandidates <- c("PG.Quantity", "PG.Normalised", "PG.MaxLFQ")
    precursorQuantityCandidates <- c("Precursor.Quantity", "Precursor.Normalised", "Precursor.Translated")

    inferQuantity <- function(featureCol) {
        if (identical(featureCol, proteinFeatureCol)) {
            diannSelectFirstPresent(proteinQuantityCandidates, reportColumns)
        }
        else if (identical(featureCol, precursorFeatureCol)) {
            diannSelectFirstPresent(precursorQuantityCandidates, reportColumns)
        }
        else {
            NULL
        }
    }

    featureCol <- NULL
    quantityCol <- NULL

    if (!is.null(diannFeatureCol)) {
        if (!(diannFeatureCol %in% reportColumns)) {
            stop("DIANN report file is missing requested feature column: ", diannFeatureCol)
        }
        featureCol <- diannFeatureCol
        if (!is.null(diannQuantityCol)) {
            if (!(diannQuantityCol %in% reportColumns)) {
                stop("DIANN report file is missing requested quantity column: ", diannQuantityCol)
            }
            quantityCol <- diannQuantityCol
        }
        else {
            quantityCol <- inferQuantity(featureCol)
            if (is.null(quantityCol)) {
                stop(
                    "Could not infer DIANN quantity column for feature column '", featureCol, "'. ",
                    "Provide diannQuantityCol explicitly."
                )
            }
        }
    }
    else if (!is.null(diannQuantityCol)) {
        if (!(diannQuantityCol %in% reportColumns)) {
            stop("DIANN report file is missing requested quantity column: ", diannQuantityCol)
        }

        quantityCol <- diannQuantityCol
        featureCol <- if (identical(diannLevel, "protein")) {
            proteinFeatureCol
        }
        else if (identical(diannLevel, "precursor")) {
            precursorFeatureCol
        }
        else if (proteinFeatureCol %in% reportColumns) {
            proteinFeatureCol
        }
        else {
            precursorFeatureCol
        }

        if (!(featureCol %in% reportColumns)) {
            stop("DIANN report file is missing requested feature column: ", featureCol)
        }
    }
    else {
        if (identical(diannLevel, "protein") || identical(diannLevel, "auto")) {
            if (proteinFeatureCol %in% reportColumns) {
                featureCol <- proteinFeatureCol
                quantityCol <- diannSelectFirstPresent(proteinQuantityCandidates, reportColumns)
            }
        }

        if (is.null(quantityCol) && (identical(diannLevel, "precursor") || identical(diannLevel, "auto"))) {
            if (precursorFeatureCol %in% reportColumns) {
                featureCol <- precursorFeatureCol
                quantityCol <- diannSelectFirstPresent(precursorQuantityCandidates, reportColumns)
            }
        }

        if (is.null(quantityCol) || !(featureCol %in% reportColumns)) {
            stop(
                "Could not infer DIANN report feature/quantity columns. ",
                "Expected protein-level columns like 'Protein.Group' + 'PG.Quantity'/'PG.MaxLFQ' ",
                "or precursor-level columns like 'Precursor.Id' + 'Precursor.Quantity'."
            )
        }
    }

    extraColsCandidates <- if (identical(featureCol, proteinFeatureCol)) {
        c("Protein.Group", "Protein.Ids", "Protein.Names", "Genes", "First.Protein.Description")
    }
    else if (identical(featureCol, precursorFeatureCol)) {
        c(
            "Precursor.Id",
            "Modified.Sequence",
            "Stripped.Sequence",
            "Precursor.Charge",
            "Protein.Group",
            "Protein.Names",
            "Genes",
            "Proteotypic"
        )
    }
    else {
        c(featureCol)
    }

    decoyCol <- if ("Decoy" %in% reportColumns) "Decoy" else NULL

    list(
        sampleCol=sampleCol,
        featureCol=featureCol,
        quantityCol=quantityCol,
        extraCols=intersect(extraColsCandidates, reportColumns),
        decoyCol=decoyCol
    )
}

diannReadReportTSV <- function(filePath, sep="\t", selectCols) {
    if (requireNamespace("data.table", quietly=TRUE)) {
        suppressWarnings(
            data.table::fread(
                filePath,
                sep=sep,
                select=selectCols,
                data.table=FALSE,
                showProgress=FALSE
            )
        )
    }
    else {
        header <- diannReadHeader(filePath, sep=sep)
        colClasses <- rep("NULL", length(header))
        keep <- header %in% selectCols
        colClasses[keep] <- "character"
        utils::read.table(
            filePath,
            sep=sep,
            header=TRUE,
            quote="",
            comment.char="",
            check.names=FALSE,
            stringsAsFactors=FALSE,
            colClasses=colClasses,
            na.strings=c("NA", "null", "")
        )
    }
}

diannReadReportParquet <- function(filePath, selectCols) {
    if (!requireNamespace("arrow", quietly=TRUE)) {
        stop(
            "Reading DIANN parquet files requires the optional 'arrow' package. ",
            "Install it, or export DIANN output as TSV instead."
        )
    }
    reader <- arrow::ParquetFileReader$create(filePath)
    schemaNames <- reader$GetSchema()$names

    selectCols <- unique(selectCols)
    selectIndices <- match(selectCols, schemaNames)
    missing <- selectCols[is.na(selectIndices)]
    if (length(missing) > 0) {
        stop(
            "DIANN parquet file is missing expected columns: ",
            paste(utils::head(missing, 10), collapse=", ")
        )
    }

    arrow::read_parquet(filePath, col_select=as.integer(selectIndices), as_data_frame=TRUE)
}

diannFilterDecoys <- function(reportDf, decoyCol="Decoy") {
    if (is.null(decoyCol) || !(decoyCol %in% colnames(reportDf))) {
        return(reportDf)
    }

    decoy <- reportDf[[decoyCol]]
    keep <- rep(TRUE, length(decoy))

    if (is.logical(decoy)) {
        keep <- is.na(decoy) | !decoy
    }
    else if (is.numeric(decoy)) {
        keep <- is.na(decoy) | decoy == 0
    }
    else {
        decoyStr <- tolower(as.character(decoy))
        keep <- is.na(decoyStr) | !(decoyStr %in% c("1", "true", "t", "yes", "y"))
    }

    reportDf[keep, , drop=FALSE]
}

diannDefaultQValueCols <- function(featureCol) {
    if (identical(featureCol, "Precursor.Id")) {
        return(c("Q.Value", "Lib.Q.Value", "Lib.PG.Q.Value"))
    }
    if (identical(featureCol, "Protein.Group")) {
        return(c("PG.Q.Value", "Global.PG.Q.Value", "Lib.PG.Q.Value"))
    }
    c("Q.Value", "Lib.Q.Value")
}

diannResolveQValueCols <- function(diannQValueCols, reportColumns, featureCol) {
    if (is.null(diannQValueCols)) {
        candidates <- diannDefaultQValueCols(featureCol)
        return(intersect(candidates, reportColumns))
    }

    diannQValueCols <- as.character(diannQValueCols)
    if (length(diannQValueCols) == 0) {
        return(character())
    }
    if (length(diannQValueCols) == 1 && identical(tolower(diannQValueCols), "auto")) {
        candidates <- diannDefaultQValueCols(featureCol)
        return(intersect(candidates, reportColumns))
    }
    if (length(diannQValueCols) == 1 && identical(tolower(diannQValueCols), "none")) {
        return(character())
    }

    missing <- setdiff(diannQValueCols, reportColumns)
    if (length(missing) > 0) {
        stop(
            "DIANN report file is missing requested q-value columns: ",
            paste(utils::head(missing, 10), collapse=", ")
        )
    }
    diannQValueCols
}

diannResolveQValueCutoffs <- function(diannQValueCols, diannQValueCutoffs) {
    if (length(diannQValueCols) == 0) {
        return(numeric())
    }
    if (is.null(diannQValueCutoffs) || length(diannQValueCutoffs) == 0) {
        return(rep(0.01, length(diannQValueCols)))
    }
    diannQValueCutoffs <- as.numeric(diannQValueCutoffs)
    diannQValueCutoffs <- diannQValueCutoffs[!is.na(diannQValueCutoffs)]
    if (length(diannQValueCutoffs) == 0) {
        return(rep(0.01, length(diannQValueCols)))
    }
    if (length(diannQValueCutoffs) != length(diannQValueCols)) {
        diannQValueCutoffs <- rep_len(diannQValueCutoffs[1], length(diannQValueCols))
    }
    diannQValueCutoffs
}

diannFilterByQValue <- function(reportDf, qCols, qCutoffs) {
    if (length(qCols) == 0) {
        return(reportDf)
    }

    keep <- rep(TRUE, nrow(reportDf))
    for (i in seq_along(qCols)) {
        q <- reportDf[[qCols[i]]]
        if (!is.numeric(q)) {
            q <- suppressWarnings(as.numeric(q))
        }
        keep[q > qCutoffs[i]] <- FALSE
    }

    reportDf[keep, , drop=FALSE]
}

diannResolveRTCol <- function(diannRTCol, reportColumns, featureCol) {
    if (!identical(featureCol, "Precursor.Id")) {
        return(NULL)
    }
    if (is.null(diannRTCol) || length(diannRTCol) == 0) {
        return(NULL)
    }
    diannRTCol <- as.character(diannRTCol)[1]
    if (is.na(diannRTCol) || diannRTCol == "") {
        return(NULL)
    }
    if (diannRTCol %in% reportColumns) {
        return(diannRTCol)
    }
    NULL
}

diannReportToWide <- function(reportDf,
                              sampleCol,
                              featureCol,
                              quantityCol,
                              extraCols,
                              designSampleNames=NULL,
                              diannMinPositive=0,
                              rtCol=NULL) {
    reportDf[[sampleCol]] <- as.character(reportDf[[sampleCol]])
    reportDf[[featureCol]] <- as.character(reportDf[[featureCol]])

    diannMinPositive <- diannResolveMinPositive(diannMinPositive)

    quantities <- reportDf[[quantityCol]]
    if (!is.numeric(quantities)) {
        quantities <- suppressWarnings(as.numeric(quantities))
    }
    quantities[quantities == 0] <- NA_real_
    if (diannMinPositive > 0) {
        quantities[!is.na(quantities) & quantities < diannMinPositive] <- NA_real_
    }

    samples <- unique(reportDf[[sampleCol]])
    if (!is.null(designSampleNames) && !all(designSampleNames %in% samples)) {
        cleanedSamples <- diannCleanSampleName(samples)
        sampleMap <- stats::setNames(cleanedSamples, samples)
        mapped <- sampleMap[reportDf[[sampleCol]]]
        mapped <- as.character(mapped)
        mapped[is.na(mapped)] <- reportDf[[sampleCol]][is.na(mapped)]

        mappedUnique <- unique(mapped)
        if (all(designSampleNames %in% mappedUnique)) {
            reportDf[[sampleCol]] <- mapped
            samples <- unique(reportDf[[sampleCol]])
        }
    }

    features <- unique(reportDf[[featureCol]])
    wide <- matrix(NA_real_, nrow=length(features), ncol=length(samples))
    colnames(wide) <- samples

    sampleIndex <- match(reportDf[[sampleCol]], samples)
    featureIndex <- match(reportDf[[featureCol]], features)
    linearIndex <- featureIndex + (sampleIndex - 1L) * length(features)
    keep <- !is.na(linearIndex) & !is.na(quantities)
    linearIndex <- linearIndex[keep]
    quantities <- quantities[keep]

    if (length(linearIndex) > 0) {
        if (anyDuplicated(linearIndex)) {
            maxByIndex <- tapply(quantities, linearIndex, max, na.rm=TRUE)
            wide[as.integer(names(maxByIndex))] <- maxByIndex
        }
        else {
            wide[linearIndex] <- quantities
        }
    }

    dedup <- !duplicated(reportDf[[featureCol]])
    annotation <- reportDf[dedup, extraCols, drop=FALSE]
    annotation <- annotation[match(features, annotation[[featureCol]]), , drop=FALSE]

    if (!is.null(rtCol) && (rtCol %in% colnames(reportDf))) {
        rtValues <- reportDf[[rtCol]]
        if (!is.numeric(rtValues)) {
            rtValues <- suppressWarnings(as.numeric(rtValues))
        }
        rtByFeature <- tapply(rtValues, reportDf[[featureCol]], stats::median, na.rm=TRUE)
        annotation[["RT"]] <- as.numeric(rtByFeature[features])
    }

    data.frame(annotation, as.data.frame(wide, check.names=FALSE), check.names=FALSE)
}

readDiannToDataFrame <- function(filePath,
                                 sep="\t",
                                 designSampleNames=NULL,
                                 diannLevel=c("auto", "protein", "precursor"),
                                 diannSampleCol=NULL,
                                 diannFeatureCol=NULL,
                                 diannQuantityCol=NULL,
                                 diannFilterDecoy=TRUE,
                                 diannFilterQValue=TRUE,
                                 diannQValueCols=NULL,
                                 diannQValueCutoffs=0.01,
                                 diannMinPositive=0,
                                 diannRTCol="RT") {

    diannLevel <- match.arg(diannLevel)

    if (diannIsParquet(filePath)) {
        if (!requireNamespace("arrow", quietly=TRUE)) {
            stop(
                "Reading DIANN parquet files requires the optional 'arrow' package. ",
                "Install it, or export DIANN output as TSV instead."
            )
        }

        reportColumns <- arrow::ParquetFileReader$create(filePath)$GetSchema()$names
        spec <- diannChooseReportSpec(
            reportColumns,
            diannLevel=diannLevel,
            diannSampleCol=diannSampleCol,
            diannFeatureCol=diannFeatureCol,
            diannQuantityCol=diannQuantityCol
        )

        qCols <- if (isTRUE(diannFilterQValue)) {
            diannResolveQValueCols(diannQValueCols, reportColumns, spec$featureCol)
        } else {
            character()
        }
        qCutoffs <- diannResolveQValueCutoffs(qCols, diannQValueCutoffs)
        rtCol <- diannResolveRTCol(diannRTCol, reportColumns, spec$featureCol)

        selectCols <- unique(c(spec$sampleCol, spec$featureCol, spec$quantityCol, spec$extraCols, spec$decoyCol, qCols, rtCol))
        reportDf <- diannReadReportParquet(filePath, selectCols=selectCols)
        if (isTRUE(diannFilterDecoy)) {
            reportDf <- diannFilterDecoys(reportDf, spec$decoyCol)
        }
        if (length(qCols) > 0) {
            reportDf <- diannFilterByQValue(reportDf, qCols=qCols, qCutoffs=qCutoffs)
        }
        return(
            diannReportToWide(
                reportDf,
                sampleCol=spec$sampleCol,
                featureCol=spec$featureCol,
                quantityCol=spec$quantityCol,
                extraCols=spec$extraCols,
                designSampleNames=designSampleNames,
                diannMinPositive=diannMinPositive,
                rtCol=rtCol
            )
        )
    }

    header <- diannReadHeader(filePath, sep=sep)
    if (diannIsReportHeader(header)) {
        spec <- diannChooseReportSpec(
            header,
            diannLevel=diannLevel,
            diannSampleCol=diannSampleCol,
            diannFeatureCol=diannFeatureCol,
            diannQuantityCol=diannQuantityCol
        )

        qCols <- if (isTRUE(diannFilterQValue)) {
            diannResolveQValueCols(diannQValueCols, header, spec$featureCol)
        } else {
            character()
        }
        qCutoffs <- diannResolveQValueCutoffs(qCols, diannQValueCutoffs)
        rtCol <- diannResolveRTCol(diannRTCol, header, spec$featureCol)

        selectCols <- unique(c(spec$sampleCol, spec$featureCol, spec$quantityCol, spec$extraCols, spec$decoyCol, qCols, rtCol))
        reportDf <- diannReadReportTSV(filePath, sep=sep, selectCols=selectCols)
        if (isTRUE(diannFilterDecoy)) {
            reportDf <- diannFilterDecoys(reportDf, spec$decoyCol)
        }
        if (length(qCols) > 0) {
            reportDf <- diannFilterByQValue(reportDf, qCols=qCols, qCutoffs=qCutoffs)
        }
        return(
            diannReportToWide(
                reportDf,
                sampleCol=spec$sampleCol,
                featureCol=spec$featureCol,
                quantityCol=spec$quantityCol,
                extraCols=spec$extraCols,
                designSampleNames=designSampleNames,
                diannMinPositive=diannMinPositive,
                rtCol=rtCol
            )
        )
    }

    matrixDf <- utils::read.table(
        filePath,
        sep=sep,
        header=TRUE,
        quote="",
        comment.char="",
        check.names=FALSE,
        stringsAsFactors=FALSE
    )

    diannMinPositive <- diannResolveMinPositive(diannMinPositive)

    matrixDf <- diannRenameSampleColumnsForDesign(matrixDf, designSampleNames)
    sampleCols <- intersect(as.character(designSampleNames), colnames(matrixDf))
    if (length(sampleCols) > 0) {
        for (colName in sampleCols) {
            if (is.numeric(matrixDf[[colName]])) {
                matrixDf[[colName]][matrixDf[[colName]] == 0] <- NA_real_
                if (diannMinPositive > 0) {
                    matrixDf[[colName]][!is.na(matrixDf[[colName]]) & matrixDf[[colName]] < diannMinPositive] <- NA_real_
                }
            }
        }
    }

    matrixDf
}
