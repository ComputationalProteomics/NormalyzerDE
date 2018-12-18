designPath <- system.file(package="NormalyzerDE", "extdata", "tiny_design.tsv")
fourColDesignPath <- system.file(package="NormalyzerDE", "extdata", "four_col_design.tsv")
designSingleRepPath <- system.file(package="NormalyzerDE", "extdata", "tiny_design_singlerep.tsv")
designSingleCondPath <- system.file(package="NormalyzerDE", "extdata", "tiny_design_singlecond.tsv")

dataPath <- system.file(package="NormalyzerDE", "extdata", "tiny_data.tsv")
dataPathlog2 <- system.file(package="NormalyzerDE", "extdata", "tiny_data_log2.tsv")
dataPathNonNA <- system.file(package="NormalyzerDE", "extdata", "tiny_data_nonna.tsv")
dataPathMaxQuantPep <- system.file(package="NormalyzerDE", "extdata", "mq_peptides_100.txt")
dataPathMaxQuantProt <- system.file(package="NormalyzerDE", "extdata", "mq_proteinGroups_100.txt")
dataPathProteios <- system.file(package="NormalyzerDE", "extdata", "tiny_data_proteios.tsv")
dataPathSingleAnnot <- system.file(package="NormalyzerDE", "extdata", "tiny_data_single_annot.tsv")
dataPathNoAnnot <- system.file(package="NormalyzerDE", "extdata", "tiny_data_no_annot.tsv")
fourColDataPath <- system.file(package="NormalyzerDE", "extdata", "four_col_data.tsv")

tempOut <- tempdir()

referenceNormResultsDir <- system.file(package="NormalyzerDE", "extdata", "unit_test_run_norm_reference")

statResultsDirNormal <- system.file(package="NormalyzerDE", "extdata", "stat_runs/normal")
statResultsTechRep <- system.file(package="NormalyzerDE", "extdata", "stat_runs/techrep")
statResultsDirBatch <- system.file(package="NormalyzerDE", "extdata", "stat_runs/batch")
statResultsDirLog2 <- system.file(package="NormalyzerDE", "extdata", "stat_runs/log2")
statResultsDirMultipleComps <- system.file(package="NormalyzerDE", "extdata", "stat_runs/multiple_comps")
statResultsDirWelch <- system.file(package="NormalyzerDE", "extdata", "stat_runs/welch")

# Note: During testing the expect_silent uses sink() to capture all output
# For debugging using the browser() statement while in sink - run
# closeAllConnections()

# When assigned TRUE forceAll will run all tests
forceAll <- TRUE

### Normal runs ###
runAllNormalization <- FALSE
normalRun <- FALSE || runAllNormalization || forceAll
singleRepRun <- FALSE || runAllNormalization || forceAll
singleCondRun <- FALSE || runAllNormalization || forceAll
nonNAEmptyRun <- FALSE || runAllNormalization || forceAll
logTransformedRun <- FALSE || runAllNormalization || forceAll
maxQuantPepRun <- FALSE || runAllNormalization || forceAll
maxQuantProtRun <- FALSE || runAllNormalization || forceAll
proteiosRun <- FALSE || runAllNormalization || forceAll
singleAnnotColRun <- FALSE || runAllNormalization || forceAll
noAnnotRun <- FALSE || runAllNormalization || forceAll

### Stats report runs ###
runAllStats <- FALSE
statisticsRunNormal <- FALSE || runAllStats || forceAll
statisticsRunSingleAnnot <- FALSE || runAllStats || forceAll
statisticsRunNoAnnot <- FALSE || runAllStats || forceAll
statisticsRunCustomThreshold <- FALSE || runAllStats || forceAll
statisticsRunWelch <- FALSE || runAllStats || forceAll
statisticsRunBatch <- FALSE || runAllStats || forceAll
statisticsRunTechRepRed <- FALSE || runAllStats || forceAll
statisticsLogTransform <- FALSE || runAllStats || forceAll
statisticsMultipleComparisons <- TRUE || runAllStats || forceAll

### SummarizedExperiments runs ###
summarizedExperimentsRun <- FALSE || forceAll

are_matrices_identical <- function(label, samples, expected_path, found_path, custom_annot=NULL, stat_cols=NULL, decimals=5) {
    
    exp_df <- read.csv(expected_path, sep="\t", stringsAsFactors=FALSE, check.names=FALSE)
    found_df <- read.csv(found_path, sep="\t", stringsAsFactors=FALSE, check.names=FALSE)
    
    if (is.null(custom_annot)) {
        expect_true(
            all(dim(exp_df) == dim(found_df)), 
            paste("Expected identical dimensions for", label, 
                  "\nFound:", paste(dim(exp_df), collapse=", "), "and", paste(dim(found_df), collapse=", "))
        )
        
        expect_true(
            all(colnames(exp_df) == colnames(found_df)),
            paste("Expected identical column names for", label,
                  "\nFound:\n", paste(colnames(found_df), collapse=", "), 
                  "\nexpected:\n", 
                  paste(colnames(exp_df), collapse=", "))
        )
    }
    else {
        
        expect_true(
            nrow(exp_df) == nrow(found_df), 
            paste("Expected same row count for", label, 
                  "\nFound:", nrow(exp_df), "and", nrow(found_df))
        )
        
        expect_true(
            all(colnames(custom_annot) %in% colnames(found_df)),
            paste("Custom annotation not present in custom file:",
                  "\nFound: ", paste(colnames(found_df), collapse=", "),
                  "\nCustom:", paste(colnames(custom_annot), collapse=", ")
            )
        )
    }

    expDfColSums <- round(colSums(exp_df[, samples], na.rm=TRUE), decimals)
    foundDfColSums <- round(colSums(found_df[, samples], na.rm=TRUE), decimals)

    expect_true(
        all(expDfColSums == foundDfColSums), 
        paste("Expected identical column sums for", label, 
              "\nFound:\n", paste(expDfColSums, collapse="\t"), "\n", paste(foundDfColSums, collapse="\t"))
    )
    
    foundAnnots <- found_df[, !(colnames(found_df) %in% c(samples, stat_cols)), drop=FALSE]
    if (is.null(custom_annot)) {
        expAnnots <- exp_df[, !(colnames(exp_df) %in% c(samples, stat_cols)), drop=FALSE]
    }
    else {
        expAnnots <- custom_annot
    }
    
    expect_true(
        all(expAnnots == foundAnnots), 
        paste("Expected annotations for", label, 
              "\nFound:", paste(head(expAnnots), collapse=", "), "and", paste(foundAnnots, collapse=", "))
    )
    
    if (!is.null(stat_cols)) {
        
        exp_stat_cols <- round(exp_df[, stat_cols], decimals)
        found_stat_cols <- round(found_df[, stat_cols], decimals)
        
        expect_true(
            all(exp_stat_cols == found_stat_cols, na.rm=TRUE),
            paste("Expected identical statistics for", label, 
                  "\nFound colsums:\n", 
                  paste(colSums(exp_stat_cols, na.rm=TRUE), collapse="\t"), 
                  "\nExpected colsums:\n", 
                  paste(colSums(foundDfColSums, na.rm=TRUE), collapse="\t")
            )
        )
    }
}

compare_output_directories <- function(label, samples, dir1, dir2, expected_count=NULL, 
                                       ignores=c(), custom_annot=NULL, stat_cols=NULL) {

    currFiles <- list.files(dir1 , pattern="*(.txt|.tsv)", full.names=TRUE)
    refFiles <- list.files(dir2, pattern="*(.txt|.tsv)", full.names=TRUE)

    if (length(which(basename(currFiles) %in% basename(refFiles))) == 0) {
        stop("Expected observed files also present in reference, found: \n",
             paste(currFiles, collapse=", "),
             "\nReference:\n",
             paste(refFiles, collapse=", "))
    }
    
    names(currFiles) <- basename(currFiles)
    names(refFiles) <- basename(refFiles)

    if (is.null(expected_count)) {
        expect_true(
            length(currFiles) == length(refFiles), 
            paste0("Number of normalized matrices (", length(currFiles), 
                   ") should be same as reference matrices (", length(refFiles), ")")
        )
    }
    else {
        expect_true(
            length(currFiles) == expected_count, 
            paste0("Number of normalized matrices (", length(currFiles), 
                   ") should be same as provided count (", expected_count, ")")
        )
    }
    
    targetCurrFiles <- currFiles[!(names(currFiles) %in% ignores)]
    for (name in names(targetCurrFiles)) {
        are_matrices_identical(name, samples, found_path=currFiles[[name]], expected_path=refFiles[[name]], 
                               custom_annot=custom_annot, stat_cols=stat_cols)
    }
}

context("Normalization runs: Normal run")
if (normalRun) {
    test_that("Normalization run succeeds quietly", {
        
        expect_silent(
            normalyzer(
                jobName="unit_test_run_norm",
                dataPath=dataPath,
                designPath=designPath,
                outputDir=tempOut,
                quiet=TRUE
            )
        )
    })
    
    test_that("Normalization results are identical to previous", {
        
        designDf <- read.csv(designPath, sep="\t")
        samples <- as.character(designDf$sample)
        currOutDir <- paste0(tempOut, "/unit_test_run_norm")
        compare_output_directories("NormalRun", samples, currOutDir, referenceNormResultsDir)
    })
}

context("Normalization runs: Single replicate run")
if (singleRepRun) {
    test_that("Normalization run succeeds without errors (single replicates)", {
        
        expect_silent(
            normalyzer(
                jobName="unit_test_run_norm_singlerep",
                dataPath=dataPath,
                designPath=designSingleRepPath,
                outputDir=tempOut,
                quiet=TRUE,
                requireReplicates=FALSE
            )
        )
    })
    
    test_that("Single-replicate run output has identical normalizations", {
        
        designDf <- read.csv(designSingleRepPath, sep="\t")
        samples <- as.character(designDf$sample)
        currOutDir <- paste0(tempOut, "/unit_test_run_norm_singlerep")
        compare_output_directories("SingleRepRun", samples, currOutDir, referenceNormResultsDir)
    })
}

context("Normalization runs: Single condition run")
if (singleCondRun) {
    test_that("Normalization run succeeds without errors (single condition)", {
        
        expect_silent(
            normalyzer(
                jobName="unit_test_run_norm_singlecond",
                dataPath=dataPath,
                designPath=designSingleCondPath,
                outputDir=tempOut,
                quiet=TRUE,
                requireReplicates=FALSE
            )
        )
    })
    
    test_that("Single-condition run output has identical normalizations", {
        
        designDf <- read.csv(designSingleCondPath, sep="\t")
        samples <- as.character(designDf$sample)
        currOutDir <- paste0(tempOut, "/unit_test_run_norm_singlecond")
        compare_output_directories("SingleCondRun", samples, currOutDir, referenceNormResultsDir)
    })
    
}

context("Normalization runs: Zeroes and empty cells instead of NAs")
if (nonNAEmptyRun) {
    test_that("Normalization run succeeds without errors (Non-NA)", {
        
        expect_silent(
            normalyzer(
                jobName="unit_test_run_norm_nonna",
                dataPath=dataPathNonNA,
                designPath=designPath,
                outputDir=tempOut,
                quiet=TRUE
            )
        )
    })
    
    test_that("Non-NA run output has identical normalizations", {
        
        designDf <- read.csv(designPath, sep="\t")
        samples <- as.character(designDf$sample)
        currOutDir <- paste0(tempOut, "/unit_test_run_norm_nonna")
        compare_output_directories("NonNARun", samples, currOutDir, referenceNormResultsDir)
    })
    
}

context("Normalization runs: Log2-transformed input")
if (logTransformedRun) {
    test_that("Normalization run succeeds without errors (already log transformed)", {
        expect_silent(
            normalyzer(
                jobName="unit_test_run_norm_log_transformed",
                dataPath=dataPathlog2,
                designPath=designPath,
                outputDir=tempOut,
                quiet=TRUE,
                noLogTransform=TRUE
            )
        )
    })
    
    test_that("Log2-run output has identical normalizations", {
        
        designDf <- read.csv(designPath, sep="\t")
        samples <- as.character(designDf$sample)
        currOutDir <- paste0(tempOut, "/unit_test_run_norm_log_transformed")
        compare_output_directories(
            "Log2-transformed run", 
            samples, 
            currOutDir, 
            referenceNormResultsDir, 
            expected_count=11,
            ignores="submitted_rawdata.txt"
        )
    })
}

context("Normalization runs: MaxQuant peptide")
if (maxQuantPepRun) {
    
    test_that("Normalization run succeeds without errors (MaxQuant pep)", {
        expect_silent(
            normalyzer(
                jobName="unit_test_run_norm_maxquantpep",
                dataPath=dataPathMaxQuantPep,
                designPath=designPath,
                outputDir=tempOut,
                quiet=TRUE,
                inputFormat = "maxquantpep",
                normalizeRetentionTime = FALSE
            )
        )
    })
    
    test_that("MaxQuant peptide output has identical normalizations", {
        
        designDf <- read.csv(designPath, sep="\t")
        samples <- as.character(designDf$sample)
        
        rawDf <- read.csv(dataPathMaxQuantPep, sep="\t", stringsAsFactors=FALSE, comment.char="", quote="", header=TRUE)
        annotDf <- rawDf[, c("Sequence", "Mass", "Proteins", "Leading.razor.protein", "PEP", "Charges")]
        
        currOutDir <- paste0(tempOut, "/unit_test_run_norm_maxquantpep")
        compare_output_directories(
            "MaxQuant peptide run", 
            samples, 
            currOutDir, 
            referenceNormResultsDir,
            expected_count=9, 
            custom_annot=annotDf
        )
    })
}

context("Normalization runs: MaxQuant protein")
if (maxQuantProtRun) {
    test_that("Normalization run succeeds without errors (MaxQuant prot)", {
        expect_silent(
            normalyzer(
                jobName="unit_test_run_norm_maxquantprot",
                dataPath=dataPathMaxQuantProt,
                designPath=designPath,
                outputDir=tempOut,
                quiet=TRUE,
                inputFormat = "maxquantprot",
                normalizeRetentionTime = FALSE
            )
        )
    })
    
    test_that("MaxQuant protein output has identical normalizations", {
        
        designDf <- read.csv(designPath, sep="\t")
        samples <- as.character(designDf$sample)
        
        rawDf <- read.csv(dataPathMaxQuantProt, sep="\t", stringsAsFactors=FALSE, comment.char="", quote="", header=TRUE)
        annotDf <- rawDf[, c("Protein.IDs", "Majority.protein.IDs", "Fasta.headers")]
        
        currOutDir <- paste0(tempOut, "/unit_test_run_norm_maxquantprot")
        compare_output_directories(
            "MaxQuant protein run", 
            samples, 
            currOutDir, 
            referenceNormResultsDir,
            expected_count=9, 
            custom_annot=annotDf
        )
    })
}

context("Normalization runs: Proteios")
if (proteiosRun) {
    test_that("Normalization run succeeds without errors (Proteios)", {
        
        expect_silent(
            normalyzer(
                jobName="unit_test_run_norm_proteios",
                dataPath=dataPathProteios,
                designPath=designPath,
                outputDir=tempOut,
                quiet=TRUE,
                inputFormat="proteios"
            )
        )
    })
    
    test_that("Proteios protein output has identical normalizations", {
        
        designDf <- read.csv(designPath, sep="\t")
        samples <- as.character(designDf$sample)
        
        currOutDir <- paste0(tempOut, "/unit_test_run_norm_proteios")
        compare_output_directories(
            "Proteios run", 
            samples, 
            currOutDir, 
            referenceNormResultsDir
        )
    })
}

context("Normalization runs: Single annotation column")
if (singleAnnotColRun) {
    
    test_that("Normalization run succeeds without errors (Single annotation column)", {
        expect_silent(
            normalyzer(
                jobName="unit_test_run_norm_singleannot",
                dataPath=dataPathSingleAnnot,
                designPath=designPath,
                outputDir=tempOut,
                quiet=TRUE
            )
        )
    })
    
    test_that("Single annotation column output has identical output", {
        
        designDf <- read.csv(designPath, sep="\t")
        samples <- as.character(designDf$sample)
        rawDf <- read.csv(dataPathSingleAnnot, sep="\t", stringsAsFactors=FALSE, 
                          comment.char="", quote="", header=TRUE, check.names=FALSE)
        annotDf <- rawDf[, "Average RT", drop=FALSE]
        currOutDir <- paste0(tempOut, "/unit_test_run_norm_singleannot")
        
        compare_output_directories(
            "Single annotation run", 
            samples, 
            currOutDir, 
            referenceNormResultsDir,
            custom_annot=annotDf
        )
    })
}

context("Normalization runs: No annotation column")
if (noAnnotRun) {
    test_that("Normalization run succeeds without errors (no annotation column)", {
        
        expect_silent(
            normalyzer(
                jobName="unit_test_run_norm_noannot",
                dataPath=dataPathNoAnnot,
                designPath=designPath,
                outputDir=tempOut,
                quiet=TRUE
            )
        )
    })
    
    test_that("No annotation column output has identical output", {
        
        designDf <- read.csv(designPath, sep="\t")
        samples <- as.character(designDf$sample)
        currOutDir <- paste0(tempOut, "/unit_test_run_norm_noannot")
        
        compare_output_directories(
            "No annotation run", 
            samples, 
            currOutDir, 
            referenceNormResultsDir,
            custom_annot=as.data.frame(matrix(nrow=100, ncol=0)),
            expected_count=9
        )
    })
}

context("Statistics runs: Normal")
if (statisticsRunNormal) {
    
    currName <- "normal"
    currDesignPath <- designPath
    currDataPath <- dataPath
    currRefDir <- statResultsDirNormal
    
    designDf <- read.csv(currDesignPath, sep="\t")
    samples <- as.character(designDf$sample)
    currOutDir <- paste0(tempOut, "/", currName)
    statCols <- c("4-5_PValue", "4-5_AdjPVal", "4-5_log2FoldChange", "featureAvg")
    
    test_that("Statistics run succeeds without errors", {

        expect_silent(
            normalyzerDE(
                jobName=currName,
                dataPath=currDataPath,
                designPath=currDesignPath,
                outputDir=tempOut,
                comparisons=c("4-5"),
                logTrans=TRUE,
                quiet=TRUE
            )   
        )
    })
    
    test_that("Statistics results are identical to previous", {
        
        compare_output_directories(
            "Stats: Normal run",
            samples,
            currOutDir,
            currRefDir,
            stat_cols=statCols
        )
    })
}

context("Statistics runs: Single annotation column")
if (statisticsRunSingleAnnot) {
    
    currName <- "normal"
    currDesignPath <- designPath
    currDataPath <- dataPathSingleAnnot
    currRefDir <- statResultsDirNormal
    
    designDf <- read.csv(currDesignPath, sep="\t")
    samples <- as.character(designDf$sample)
    currOutDir <- paste0(tempOut, "/", currName)
    statCols <- c("4-5_PValue", "4-5_AdjPVal", "4-5_log2FoldChange", "featureAvg")
    
    rawDf <- read.csv(currDataPath, sep="\t", stringsAsFactors=FALSE, 
                      comment.char="", quote="", header=TRUE, check.names=FALSE)
    annotDf <- rawDf[, "Average RT", drop=FALSE]
    statCols <- c("4-5_PValue", "4-5_AdjPVal", "4-5_log2FoldChange", "featureAvg")
    
    test_that("Statistics (single annot) run succeeds without errors", {
        
        expect_silent(
            normalyzerDE(
                jobName=currName,
                dataPath=currDataPath,
                designPath=currDesignPath,
                outputDir=tempOut,
                comparisons=c("4-5"),
                logTrans=TRUE,
                quiet=TRUE
            )   
        )
    })
    
    test_that("Statistics (single annot) results are identical to previous", {
        
        compare_output_directories(
            "Stats: Single annotation column run",
            samples,
            currOutDir,
            currRefDir,
            custom_annot=annotDf,
            stat_cols=statCols
        )
    })
}

context("Statistics runs: No annotation column")
if (statisticsRunNoAnnot) {
    
    currName <- "normal"
    currDesignPath <- designPath
    currDataPath <- dataPathNoAnnot
    currRefDir <- statResultsDirNormal
    
    designDf <- read.csv(currDesignPath, sep="\t")
    samples <- as.character(designDf$sample)
    currOutDir <- paste0(tempOut, "/", currName)
    statCols <- c("4-5_PValue", "4-5_AdjPVal", "4-5_log2FoldChange", "featureAvg")
    
    test_that("Statistics (no annot) run succeeds without errors", {
        
        expect_silent(
            normalyzerDE(
                jobName=currName,
                dataPath=currDataPath,
                designPath=currDesignPath,
                outputDir=tempOut,
                comparisons=c("4-5"),
                logTrans=TRUE,
                quiet=TRUE
            )   
        )
    })
    
    test_that("Statistics (no annot) results are identical to previous", {
        
        compare_output_directories(
            "Stats: No annotation column run",
            samples,
            currOutDir,
            statResultsDirNormal,
            custom_annot=as.data.frame(matrix(ncol=0, nrow=100)),
            stat_cols=statCols
        )
    })
}

context("Statistics runs: Custom threshold (no identity check)")
if (statisticsRunCustomThreshold) {
    
    test_that("Statistics with custom thresholding", {
        
        expect_silent(
            normalyzerDE(
                jobName="unit_test_run_stat_custom_thres",
                dataPath=dataPath,
                designPath=designPath,
                outputDir=tempOut,
                comparisons=c("4-5"),
                logTrans=TRUE,
                quiet=TRUE,
                sigThres=0.05,
                sigThresType="p",
                log2FoldThres="1"
            )   
        )
    })
}

context("Statistics runs: Welch")
if (statisticsRunWelch) {
    
    currName <- "welch"
    currDesignPath <- designPath
    currDataPath <- dataPath
    currRefDir <- statResultsDirNormal
    
    designDf <- read.csv(currDesignPath, sep="\t")
    samples <- as.character(designDf$sample)
    currOutDir <- paste0(tempOut, "/", currName)
    statCols <- c("4-5_PValue", "4-5_AdjPVal", "4-5_log2FoldChange", "featureAvg")
    
    test_that("Statistics run succeeds without errors", {
        expect_silent(
            normalyzerDE(
                jobName=currName, 
                designPath=currDesignPath, 
                dataPath=currDataPath, 
                outputDir=tempOut, 
                logTrans=TRUE, 
                comparisons=c("4-5"), 
                leastRepCount=2, 
                sampleCol="sample", 
                condCol="group", 
                type="welch",
                quiet=TRUE
            )
        )
    })
    
    test_that("Statistics (welch) results are identical to previous", {
        
        compare_output_directories(
            "Stats: No annotation column run",
            samples,
            currOutDir,
            statResultsDirWelch,
            stat_cols=statCols
        )
    })
}

context("Statistics runs: Batch")
if (statisticsRunBatch) {
    
    currName <- "batch"
    currDesignPath <- designPath
    currDataPath <- dataPath
    currRefDir <- statResultsDirBatch
    
    designDf <- read.csv(currDesignPath, sep="\t")
    samples <- as.character(designDf$sample)
    currOutDir <- paste0(tempOut, "/", currName)
    statCols <- c("4-5_PValue", "4-5_AdjPVal", "4-5_log2FoldChange", "featureAvg")
    
    test_that("Statistics run succeeds without errors", {
        expect_silent(
            normalyzerDE(
                jobName=currName, 
                designPath=currDesignPath, 
                dataPath=currDataPath, 
                outputDir=tempOut, 
                logTrans=TRUE, 
                comparisons=c("4-5"), 
                leastRepCount=1, 
                sampleCol="sample", 
                condCol="group", 
                type="limma",
                batchCol="batch",
                quiet=TRUE
            )
        )
    })
    
    test_that("Statistics (batch) results are identical to previous", {
        
        compare_output_directories(
            "Stats: Batch",
            samples,
            currOutDir,
            currRefDir,
            stat_cols=statCols
        )
    })
}

context("Statistics runs: Technical replicate reduce")
if (statisticsRunTechRepRed) {
    
    currName <- "techrep"
    currDesignPath <- fourColDesignPath
    currDataPath <- fourColDataPath
    currRefDir <- statResultsTechRep
    
    designDf <- read.csv(currDesignPath, sep="\t")
    currOutDir <- paste0(tempOut, "/", currName)
    statCols <- c("4-5_PValue", "4-5_AdjPVal", "4-5_log2FoldChange", "featureAvg")
    
    designDf <- read.csv(designPath, sep="\t")
    samples <- c(
        "s_125amol_1.s_125amol_2",
        "s_125amol_3.s_125amol_4",
        "s_12500amol_1.s_12500amol_2",
        "s_12500amol_3.s_12500amol_4",
        "s_25000amol_1.s_25000amol_2",
        "s_25000amol_3.s_25000amol_4"
    )
    
    test_that("Statistics run succeeds without errors", {
        expect_silent(
            normalyzerDE(
                jobName="techrep",
                designPath=currDesignPath,
                dataPath=currDataPath,
                outputDir=tempOut,
                comparisons=c("4-5"),
                logTrans=TRUE,
                quiet=TRUE,
                techRepCol="techrep",
                sampleCol="sample",
                condCol="group",
                leastRepCount=2
            )   
        )
    })
    
    test_that("Statistics results are identical to previous", {
        compare_output_directories(
            "Stats: Techrep run",
            samples,
            currOutDir,
            currRefDir,
            stat_cols=statCols
        )
    })
    
}

context("Statistics runs: Log transformed")
if (statisticsLogTransform) {
    
    currName <- "normal"
    currDesignPath <- designPath
    currDataPath <- dataPathlog2
    currRefDir <- statResultsDirNormal
    
    designDf <- read.csv(designPath, sep="\t")
    samples <- as.character(designDf$sample)
    currOutDir <- paste0(tempOut, "/", currName)
    statCols <- c("4-5_PValue", "4-5_AdjPVal", "4-5_log2FoldChange", "featureAvg")
    
    test_that("Statistics run succeeds without errors", {
        expect_silent(
            normalyzerDE(
                jobName=currName, 
                designPath=currDesignPath, 
                dataPath=currDataPath, 
                outputDir=tempOut, 
                logTrans=FALSE, 
                comparisons=c("4-5"), 
                leastRepCount=1,
                sampleCol="sample", 
                condCol="group", 
                type="limma",
                quiet=TRUE
            )
        )
    })
    
    test_that("Statistics (log2) results are identical to previous", {
        
        compare_output_directories(
            "Stats: Log2 transformed input",
            samples,
            currOutDir,
            currRefDir,
            stat_cols=statCols
        )
    })
}

context("Statistics runs: Multiple comparisons")
if (statisticsMultipleComparisons) {
    
    currName <- "multiple_comps"
    currDesignPath <- fourColDesignPath
    currDataPath <- fourColDataPath
    currRefDir <- statResultsDirMultipleComps
    
    designDf <- read.csv(currDesignPath, sep="\t")
    samples <- as.character(designDf$sample)
    currOutDir <- paste0(tempOut, "/", currName)
    statCols <- c(
        "4-5_PValue", "4-6_PValue", "5-6_PValue", 
        "4-5_AdjPVal", "4-6_AdjPVal", "5-6_AdjPVal", 
        "4-5_log2FoldChange", "4-6_log2FoldChange", "5-6_log2FoldChange", 
        "featureAvg")
    
    test_that("Statistics run succeeds without errors", {
        expect_silent(
            normalyzerDE(
                jobName=currName, 
                designPath=currDesignPath, 
                dataPath=currDataPath, 
                outputDir=tempOut, 
                logTrans=FALSE, 
                comparisons=c("4-5", "4-6", "5-6"), 
                leastRepCount=1, 
                sampleCol="sample", 
                condCol="group", 
                type="limma",
                quiet=TRUE,
                sigThres=0.2,
                sigThresType="p",
                log2FoldThres=0.5
            )
        )
    })
    
    test_that("Statistics (multiple comparisons) results are identical to previous", {
        
        compare_output_directories(
            "Stats: Multiple comparisons",
            samples,
            currOutDir,
            currRefDir,
            stat_cols=statCols
        )
    })
}

context("SummarizedExperiments run")

if (summarizedExperimentsRun) {
    test_that("Normalization run from SummarizedExperiment succeeds without errors", {
        
        rdf <- read.csv(dataPath, sep="\t")
        ddf <- read.csv(designPath, sep="\t")
        ddf$sample <- as.character(ddf$sample)
        sdf <- rdf[, ddf$sample]
        adf <- rdf[, !(colnames(rdf) %in% ddf$sample)]
        se <- SummarizedExperiment::SummarizedExperiment(
            assays=list(raw=as.matrix(sdf)),
            colData=ddf,
            rowData=adf
        )
        
        expect_silent(
            normalyzer(
                jobName="unit_test_run_norm",
                experimentObj=se,
                outputDir=tempOut, 
                quiet=TRUE
            )
        )
    })
    
    test_that("Statistics run from SummarizedExperiment succeeds without errors", {
        
        rdf <- read.csv(dataPath, sep="\t")
        ddf <- read.csv(designPath, sep="\t")
        ddf$sample <- as.character(ddf$sample)
        sdf <- rdf[, ddf$sample]
        adf <- rdf[, !(colnames(rdf) %in% ddf$sample)]
        se <- SummarizedExperiment::SummarizedExperiment(
            assays=list(raw=as.matrix(sdf)),
            colData=ddf,
            rowData=adf
        )
        
        expect_silent(
            normalyzerDE(
                jobName="unit_test_run_stat",
                experimentObj=se,
                outputDir=tempOut,
                comparisons=c("4-5"),
                logTrans=TRUE,
                quiet=TRUE
            )   
        )
    })
}

