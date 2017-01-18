

# Input: Normalized matrices

# We want to target particular 7 vs 7 columns
# Also, we want to slice out potato-named rows into a distinctly separate dataframe

# Let's test this first in an R terminal

source("NormalyzerDataset.R")
source("NormalizationEvaluationResults.R")
source("NormalyzerResults.R")

source("utils.R")
source("inputVerification.R")
source("analyzeResults.R")
source("normMethods.R")
source("normfinder-pipeline.R")

full_report <- "FeatureReport_MQP_20140828-fornormalyzer.txt"
subset_report <- "MQP.subset_500.tsv"
# subset_report <- "BA_MQP.subset_500.tsv"

POT_PAT <- "^sol"
HUMAN_PAT <- "^hum"
SAMPLE_PAT <- "dilA_[23]"
custom_replicate_groups <- c(2,2,2,2,2,2,2,3,3,3,3,3,3,3)

run_review_data_test <- function(subset=FALSE, output_path=NULL, verbose=FALSE) {
    
    start_time <- Sys.time()
    
    if (subset) review_data_path <- paste("../tests/data", subset_report, sep="/")
    else review_data_path <- paste("../tests/data", full_report, sep="/")
    
    print(paste("Processing dataset:", review_data_path))
    
    print("Get normalyzer object")
    job_name <- "review_evaluation"
    norm_obj <- getVerifiedNormalyzerObject(review_data_path, job_name)

    print(paste("Normalizing for window:", rt_window))
    # file_name <- paste0("normalyzer_eval_rt_size_", rt_window)
    print("Before run")
    default_rt_window <- 2
    rt_windows <- c(0.2, 0.5, 1, 2, 3, 5, 8, 13, 17, 21, 27, 34)
    perform_normalyzer_run(norm_obj, job_name, default_rt_window, output_path=output_path, rt_windows=rt_windows, verbose=verbose)

    
    # nr <- analyzeNormalizations(nr, jobName)
    end_time <- Sys.time()
    print(difftime(end_time, start_time))
    print(paste("Start:", start_time, "End:", end_time))
}

perform_normalyzer_run <- function(norm_obj, job_name, default_rt_size, output_path=NULL, rt_windows=c(), verbose=FALSE) {
    
    nr <- normMethods(norm_obj, job_name, normalizeRetentionTime=FALSE)
    
    used_methods_names <- getUsedMethodNames(nr)
    normalization_matrices <- getNormalizationMatrices(nr)
    
    header <- c("method", "tot_na_reduced", "RT settings", 
                "Total potato", "Number sig. p", "Frac. p", "Number sig. FDR", "Frac. FDR",  
                "Total background", "Number sig. p", "Frac. p", "Number sig. FDR", "Frac. FDR")
    
    output <- matrix(ncol=length(header), nrow=0)
    colnames(output) <- header
    
    fdr_thres <- 0.1
    p_thres <- 0.05
    
    for (i in 1:length(used_methods_names)) {
        print(paste("Method: ", used_methods_names[[i]]))
        run_results <- get_analysis_vector(nr, normalization_matrices[[i]], used_methods_names[[i]], POT_PAT, HUMAN_PAT, p_thres=p_thres, fdr_thres=fdr_thres)
        output <- rbind(output, run_results)
    }
    
    normalyzer_filterraw <- nr@nds@filterrawdata
    retention_times <- nr@nds@retentionTimes
    
    print("rt_windows")
    print(rt_windows)
    
    for (rt_window in rt_windows) {
        
        print(paste("Method: RTs, window size:", rt_window))
        
        median_rt <- getSmoothedRTNormalizedMatrix(normalyzer_filterraw, retention_times, medianNormalization, rt_window)
        median_res <- get_analysis_vector(nr, median_rt, "RT-median", POT_PAT, HUMAN_PAT, p_thres=p_thres, fdr_thres=fdr_thres, rt_settings=rt_window)
        output <- rbind(output, median_res)
        
        mean_rt <- getSmoothedRTNormalizedMatrix(normalyzer_filterraw, retention_times, meanNormalization, rt_window)
        mean_res <- get_analysis_vector(nr, mean_rt, "RT-mean", POT_PAT, HUMAN_PAT, p_thres=p_thres, fdr_thres=fdr_thres, rt_settings=rt_window)
        output <- rbind(output, mean_res)
        
        loess_rt <- getSmoothedRTNormalizedMatrix(normalyzer_filterraw, retention_times, performCyclicLoessNormalization, rt_window)
        loess_res <- get_analysis_vector(nr, loess_rt, "RT-loess", POT_PAT, HUMAN_PAT, p_thres=p_thres, fdr_thres=fdr_thres, rt_settings=rt_window)
        output <- rbind(output, loess_res)
    }
    
    if (verbose) {
        print(output)
    }
        
    if (!is.null(output_path)) {
        write.csv(output, file=output_path, quote=FALSE)
    }
    
}

get_analysis_vector <- function(nr, method_data, name, potato_pattern, human_pattern, p_thres=0.05, fdr_thres=0.1, rt_settings="None") {
    
    combined_pattern <- paste0("(", potato_pattern, "|", human_pattern, ")")
    row_name_col <- 3

    prepared_df <- get_prepared_normalyzer_sheet(nr, method_data, SAMPLE_PAT, row_name_col)
    na_filter_df <- na_filter(prepared_df)
    anova_pval <- get_anova_pvals(na_filter_df, custom_replicate_groups)
    anova_fdr <- stats::p.adjust(anova_pval, method="BH")
    
    tot_passing <- nrow(na_filter_df)
    
    potato_pval <- anova_pval[which(grepl(potato_pattern, names(anova_pval)))]
    tot_potato <- length(potato_pval)
    potato_n_p <- length(potato_pval[which(potato_pval <= p_thres)])
    potato_p_frac <- round(100 * potato_n_p / tot_potato, 1)
    
    potato_fdr <- anova_fdr[which(grepl(potato_pattern, names(anova_fdr)))]
    potato_n_fdr <- length(potato_fdr[which(potato_fdr <= fdr_thres)])
    potato_fdr_frac <- round(100 * potato_n_fdr / tot_potato, 1)
    
    background_pval <- anova_pval[which(!grepl(combined_pattern, names(anova_pval)))]
    tot_background <- length(background_pval)
    background_n_p <- length(background_pval[which(background_pval <= p_thres)])
    background_p_frac <- round(100 * background_n_p / tot_background, 1)
    
    background_fdr <- anova_fdr[which(!grepl(combined_pattern, names(anova_fdr)))]
    background_n_fdr <- length(potato_fdr[which(background_fdr <= fdr_thres)])
    background_fdr_frac <- round(100 * background_n_fdr / tot_background, 1)

    analysis_vector <- c(name, tot_passing, rt_settings,
                         tot_potato, potato_n_p, potato_p_frac, potato_n_fdr, potato_fdr_frac, 
                         tot_background, background_n_p, background_p_frac, background_n_fdr, background_fdr_frac)

    analysis_vector
}



get_prepared_normalyzer_sheet <- function(nr, df, col_pattern, row_name_col) {
    
    raw_data <- nr@nds@rawData
    
    header_row <- 2
    names <- raw_data[-1:-2, row_name_col] 
    target_cols <- which(grepl(col_pattern, raw_data[header_row,]))
    
    
    parsed_df <- df[, target_cols]
    rownames(parsed_df) <- names
    colnames(parsed_df) <- custom_replicate_groups
    parsed_df
}



get_normalyzer_df_subset <- function(nr, df, col_pattern, row_pattern, row_name_col, 
                                     inverse_name_pattern=FALSE) {
    
    raw_data <- nr@nds@rawData

    header_row <- 2
    if (inverse_name_pattern) {
        target_rows <- which(!grepl(row_pattern, raw_data[-1:-2, row_name_col]))
    }
    else {
        target_rows <- which(grepl(row_pattern, raw_data[-1:-2, row_name_col]))
    }

    target_cols <- which(grepl(col_pattern, raw_data[header_row,]))
    df[target_rows, target_cols]
}


na_filter <- function(df) {
    na_per_line <- rowSums(is.na(df))
    df_na_removed <- df[na_per_line < ncol(df) / 2, ]
    df_na_removed
}


get_anova_pvals <- function(df, replicate_groups) {
    
    anova_func <- function(sampleIndex) {
        summary(stats::aov(unlist(sampleIndex)~replicate_groups))[[1]][[5]][1]
    }
    
    anovaPVal <- apply(df, 1, anova_func)
    anovaPVal
}

