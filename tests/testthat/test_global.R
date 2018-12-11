designPath <- system.file(package="NormalyzerDE", "extdata", "tiny_design.tsv")
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

tempOut <- tempdir()

referenceNormResultsDir <- system.file(package="NormalyzerDE", "extdata", "unit_test_run_norm_reference")
referenceStatResultsDir <- system.file(package="NormalyzerDE", "extdata", "unit_test_run_stat_reference")

# Note: During testing the expect_silent uses sink() to capture all output
# For debugging using the browser() statement while in sink - run
# closeAllConnections()

forceAll <- TRUE

### Normal runs ###
normalRun <- FALSE || forceAll
singleRepRun <- FALSE || forceAll
singleCondRun <- FALSE || forceAll
nonNAEmptyRun <- FALSE || forceAll
logTransformedRun <- FALSE || forceAll
maxQuantPepRun <- FALSE || forceAll
maxQuantProtRun <- FALSE || forceAll
proteiosRun <- FALSE || forceAll
singleAnnotColRun <- FALSE || forceAll
noAnnotRun <- FALSE || forceAll

### Stats report runs ###
statisticsRunNormal <- TRUE
statisticsRunSingleAnnot <- FALSE
statisticsRunNoAnnot <- FALSE
statisticsRunCustomThreshold <- TRUE

### SummarizedExperiments runs ###
summarizedExperimentsRun <- FALSE || forceAll

are_matrices_identical <- function(label, samples, path1, path2, custom_annot=NULL) {
    
    df1 <- read.csv(path1, sep="\t", stringsAsFactors=FALSE)
    df2 <- read.csv(path2, sep="\t", stringsAsFactors=FALSE)
    
    if (is.null(custom_annot)) {
        expect_true(
            all(dim(df1) == dim(df2)), 
            paste("Expected identical dimensions for", label, 
                  "\nFound:", paste(dim(df1), collapse=", "), "and", paste(dim(df2), collapse=", "))
        )
        
        expect_true(
            all(colnames(df1) == colnames(df2)),
            paste("Expected identical column names for", label,
                  "\nFound:\n", paste(colnames(df1), collapse=", "), "\nexpected:\n", 
                  paste(colnames(df2), collapse=", "))
        )
    }
    else {
        expect_true(
            nrow(df1) == nrow(df2), 
            paste("Expected same row count for", label, 
                  "\nFound:", nrow(df1), "and", nrow(df2))
        )
    }

    df1colSums <- round(colSums(df1[, samples], na.rm=TRUE), 4)
    df2colSums <- round(colSums(df2[, samples], na.rm=TRUE), 4)

    expect_true(
        all(df1colSums == df2colSums), 
        paste("Expected identical column sums for", label, 
              "\nFound:\n", paste(df1colSums, collapse="\t"), "\n", paste(df2colSums, collapse="\t"))
    )
    
    df1annots <- df1[, !(colnames(df1) %in% samples)]
    
    if (is.null(custom_annot)) {
        df2annots <- df2[, !(colnames(df2) %in% samples)]
    }
    else {
        df2annots <- custom_annot
    }
    
    expect_true(
        all(df1annots == df2annots), 
        paste("Expected annotations for", label, 
              "\nFound:", paste(head(df1annots), collapse=", "), "and", paste(df2annots, collapse=", "))
    )        
}

compare_output_directories <- function(label, samples, dir1, dir2, expected_count=NULL, ignores=c(), custom_annot=NULL) {
    
    currFiles <- list.files(dir1 , pattern="*.txt", full.names=TRUE)
    refFiles <- list.files(dir2, pattern="*.txt", full.names=TRUE)

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
    
    designDf <- read.csv(designPath, sep="\t")
    samples <- as.character(designDf$sample)
    
    targetCurrFiles <- currFiles[!(names(currFiles) %in% ignores)]
    
    for (name in names(targetCurrFiles)) {
        are_matrices_identical(name, samples, currFiles[[name]], refFiles[[name]], custom_annot=custom_annot)
    }
}

context("Normalization runs: Normal run")

if (normalRun) {
    test_that("Normalization run succeeds without errors", {
        
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
    test_that("Statistics run succeeds without errors", {
        
        expect_silent(
            normalyzerDE(
                jobName="unit_test_run_stat",
                dataPath=dataPath,
                designPath=designPath,
                outputDir=tempOut,
                comparisons=c("4-5"),
                logTrans=TRUE,
                quiet=TRUE
            )   
        )
    })
    
    test_that("Statistics results are identical to previous", {
        
        currOutDir <- paste0(tempOut, "/unit_test_run_stat")
        currStatFiles <- list.files(currOutDir , pattern="*.tsv", full.names=TRUE)
        refStatFiles <- list.files(referenceStatResultsDir, pattern="*.tsv", full.names=TRUE)
        expect_true(length(currStatFiles) == length(refStatFiles), "Number of normalized matrices should be same as reference")
        
        currStatMd5 <- tools::md5sum(currStatFiles)
        refStatMd5 <- tools::md5sum(refStatFiles)
        
        expect_true(all(currStatMd5 == refStatMd5), "MD5-sums for stat matrices should be equal")
    })
}

context("Statistics runs: Single annotation column")

if (statisticsRunSingleAnnot) {
    test_that("Statistics, single annotation column, no errors", {
        
        expect_silent(
            normalyzerDE(
                jobName="unit_test_run_stat",
                dataPath=dataPathSingleAnnot,
                designPath=designPath,
                outputDir=tempOut,
                comparisons=c("4-5"),
                logTrans=TRUE,
                quiet=TRUE
            )   
        )
    })
    
    test_that("Statistics, single annotation column, identical to previous", {
        
        currOutDir <- paste0(tempOut, "/unit_test_run_stat")
        currStatFiles <- list.files(currOutDir , pattern="*.tsv", full.names=TRUE)
        refStatFiles <- list.files(referenceStatResultsDir, pattern="*.tsv", full.names=TRUE)
        expect_true(length(currStatFiles) == length(refStatFiles), "Number of normalized matrices should be same as reference")
        
        currStatMd5 <- tools::md5sum(currStatFiles)
        refStatMd5 <- tools::md5sum(refStatFiles)
        
        expect_true(all(currStatMd5 == refStatMd5), "MD5-sums for stat matrices should be equal")
    })
}

context("Statistics runs: No annotation column")
if (statisticsRunNoAnnot) {
    test_that("Statistics run succeeds without errors", {
        
        expect_silent(
            normalyzerDE(
                jobName="unit_test_run_stat",
                dataPath=dataPath,
                designPath=designPath,
                outputDir=tempOut,
                comparisons=c("4-5"),
                logTrans=TRUE,
                quiet=TRUE
            )   
        )
    })
    
    test_that("Statistics results are identical to previous", {
        
        currOutDir <- paste0(tempOut, "/unit_test_run_stat")
        currStatFiles <- list.files(currOutDir , pattern="*.tsv", full.names=TRUE)
        refStatFiles <- list.files(referenceStatResultsDir, pattern="*.tsv", full.names=TRUE)
        expect_true(length(currStatFiles) == length(refStatFiles), "Number of normalized matrices should be same as reference")
        
        currStatMd5 <- tools::md5sum(currStatFiles)
        refStatMd5 <- tools::md5sum(refStatFiles)
        
        expect_true(all(currStatMd5 == refStatMd5), "MD5-sums for stat matrices should be equal")
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

