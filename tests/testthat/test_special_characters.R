context("Special character support")

test_that("Sample names with spaces are preserved", {
    
    test_data <- matrix(
        c(100, 200, 300,
          150, 250, 350,
          200, 300, 400,
          250, 350, 450),
        nrow = 4, byrow = TRUE
    )
    colnames(test_data) <- c("Sample 1", "Sample 2", "Sample 3")
    
    design <- data.frame(
        sample = c("Sample 1", "Sample 2", "Sample 3"),
        group = c("A", "A", "B")
    )
    
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(raw = test_data),
        colData = design,
        rowData = data.frame(peptide = paste0("Pep", 1:4))
    )
    S4Vectors::metadata(se) <- list(sample = "sample", group = "group")
    
    nds <- getVerifiedNormalyzerObject(
        "special_char_test",
        se,
        threshold = 1,
        requireReplicates = FALSE,
        quiet = TRUE
    )
    
    expect_true(all(c("Sample 1", "Sample 2", "Sample 3") %in% colnames(filterrawdata(nds))))
})

test_that("Unicode characters in sample names are preserved", {
    
    test_data <- matrix(
        c(100, 200, 150),
        nrow = 1
    )
    colnames(test_data) <- c("Sample α", "Sample β", "Sample γ")
    
    design <- data.frame(
        sample = c("Sample α", "Sample β", "Sample γ"),
        group = c("A", "A", "B")
    )
    
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(raw = test_data),
        colData = design,
        rowData = data.frame(feature = "F1")
    )
    S4Vectors::metadata(se) <- list(sample = "sample", group = "group")
    
    nds <- getVerifiedNormalyzerObject(
        "unicode_test",
        se,
        threshold = 1,
        requireReplicates = FALSE,
        quiet = TRUE
    )
    
    expect_true(all(c("Sample α", "Sample β", "Sample γ") %in% colnames(filterrawdata(nds))))
})

test_that("setupRawDataObject preserves special character sample names from file", {
    
    temp_data <- tempfile(fileext = ".tsv")
    temp_design <- tempfile(fileext = ".tsv")
    on.exit({
        unlink(temp_data)
        unlink(temp_design)
    })
    
    data_df <- data.frame(
        Peptide = paste0("P", 1:5),
        `Sample A-1` = runif(5, 100, 1000),
        `Sample A-2` = runif(5, 100, 1000),
        `Sample B.1` = runif(5, 100, 1000),
        check.names = FALSE
    )
    write.table(data_df, temp_data, sep = "\t", row.names = FALSE, quote = FALSE)
    
    design_df <- data.frame(
        sample = c("Sample A-1", "Sample A-2", "Sample B.1"),
        group = c("A", "A", "B")
    )
    write.table(design_df, temp_design, sep = "\t", row.names = FALSE, quote = FALSE)
    
    se <- setupRawDataObject(temp_data, temp_design, inputFormat = "default")
    data_colnames <- colnames(SummarizedExperiment::assay(se))
    
    expect_true(all(c("Sample A-1", "Sample A-2", "Sample B.1") %in% data_colnames))
})

test_that("Design matrix with custom column names preserves spaces", {
    
    temp_design <- tempfile(fileext = ".tsv")
    on.exit(unlink(temp_design))
    
    design_df <- data.frame(
        `Sample ID` = c("A-1", "A-2", "B.1"),
        `Group Name` = c("Control", "Control", "Treatment"),
        check.names = FALSE
    )
    write.table(design_df, temp_design, sep = "\t", row.names = FALSE, quote = FALSE)
    
    result <- loadDesign(temp_design, sampleCol = "Sample ID", groupCol = "Group Name")
    
    expect_true(all(c("Sample ID", "Group Name") %in% colnames(result)))
})

test_that("Special characters survive normalization + output", {
    
    test_data <- matrix(runif(50, min = 100, max = 1000), nrow = 10)
    colnames(test_data) <- c("Sample A-1", "Sample A-2", "Sample B.1", "Sample B.2", "Sample C_1")
    
    design <- data.frame(
        sample = c("Sample A-1", "Sample A-2", "Sample B.1", "Sample B.2", "Sample C_1"),
        group = c("GroupA", "GroupA", "GroupB", "GroupB", "GroupC")
    )
    
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(raw = test_data),
        colData = design,
        rowData = data.frame(peptide = paste0("Pep", 1:10))
    )
    S4Vectors::metadata(se) <- list(sample = "sample", group = "group")
    
    nds <- getVerifiedNormalyzerObject(
        "full_pipeline_test",
        se,
        threshold = 1,
        requireReplicates = FALSE,
        quiet = TRUE
    )
    nr <- suppressWarnings(normMethods(nds, quiet = TRUE))
    nr_eval <- analyzeNormalizations(nr)
    
    output_dir <- tempfile()
    dir.create(output_dir)
    on.exit(unlink(output_dir, recursive = TRUE), add = TRUE)
    
    writeNormalizedDatasets(nr_eval, output_dir)
    
    output_files <- list.files(output_dir, pattern = "-normalized.txt$", full.names = TRUE)
    output_data <- read.table(output_files[1], header = TRUE, sep = "\t",
                              check.names = FALSE, comment.char = "")
    
    expect_true(all(c("Sample A-1", "Sample B.1", "Sample C_1") %in% colnames(output_data)))
})
