

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
source("higherOrderNormMethods.R")
source("normfinder-pipeline.R")

source("EvaluationUtils/evaluationPlotting.R")
source("EvaluationUtils/RunEntry.R")

suppressPackageStartupMessages("ggplot2")
library("ggplot2")


full_report <- "FeatureReport_MQP_20140828-fornormalyzer.txt"
subset_report <- "MQP.subset_500.tsv"
# subset_report <- "BA_MQP.subset_500.tsv"

POT_PAT <- "^sol"
HUMAN_PAT <- "^hum"
SAMPLE_PAT <- "dilA_[23]"
COMBINED_PATTERN <- paste0("(", POT_PAT, "|", HUMAN_PAT, ")")
custom_replicate_groups <- c(2,2,2,2,2,2,2,3,3,3,3,3,3,3)


run_review_data_test <- function(output_base, 
                                 subset=FALSE,
                                 do_normal_run=TRUE, 
                                 do_rt_run=TRUE, 
                                 verbose=FALSE,
                                 measure_type="all", 
                                 sig_thres=0.1, 
                                 do_fdr=TRUE,
                                 frame_shifts=3) {
    
    output_path <- paste(output_base, "csv", sep=".")
    plot_path <- paste(output_base, "png", sep=".")
    print(paste("Evaluation run with sig_thres:", sig_thres, "and FDR setting:", do_fdr))
    plot_ylabs <- get_y_label(measure_type)
    score_funcs <- get_norm_methods(measure_type)
    print(paste("Processing using measure:", measure_type))
    
    start_time <- Sys.time()
    
    if (subset) review_data_path <- paste("../tests/data", subset_report, sep="/")
    else review_data_path <- paste("../tests/data", full_report, sep="/")
    
    print(paste("Processing dataset:", review_data_path))
    job_name <- "review_evaluation"
    norm_obj <- getVerifiedNormalyzerObject(review_data_path, job_name)

    if (subset) rt_windows <- c(1, 2, 3)

    else        rt_windows <- c(0.2, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 7, 8, 9, 10, 13, 16, 20)
    nr <- normMethods(norm_obj, job_name, normalizeRetentionTime=FALSE, runNormfinder=FALSE)
    normal_run_entries <- generate_normal_run_results(nr, adjust_fdr=do_fdr, sig_thres=sig_thres)
    rt_run_entries <- generate_rt_run_results(nr, rt_windows, adjust_fdr=do_fdr, sig_thres=sig_thres, frame_shifts=frame_shifts)
    
    print(paste("Normal entries:", length(normal_run_entries)))
    print(paste("RT entries:", length(rt_run_entries)))
    
    if (verbose) {
        print("--- Normal entries ---")
        for (entry in normal_run_entries) {
            print(get_entry_string(entry))
        }
        print("--- RT entries ---")
        for (entry in rt_run_entries) {
            print(get_entry_string(entry))
        }
    }
    
    print("Visualizing results...")
    print(length(score_funcs))
    
    plots <- list()
    for (i in 1:length(score_funcs)) {
        
        score_func <- score_funcs[[i]]
        y_lab <- plot_ylabs[i]
        
        # sub_plot_path <- paste(plot_path, tolower(y_lab), "png", sep=".")
        
        include_legend <- TRUE
        if (i != length(score_funcs)) {
            include_legend <- FALSE
        }
        
        plt <- visualize_results_new(normal_run_entries,
                              rt_run_entries,
                              score_func,
                              title=y_lab,
                              xlab="RT window size (minutes)",
                              ylab=y_lab,
                              verbose=verbose)
        plots[[i]] <- plt
    }

    plot_path <- paste(plot_path, "png", sep=".")
    png(plot_path)
    # multiplot(plotlist=plots, cols=length(plots))
    grid_arrange_shared_legend(plots=plots, ncol=length(plots))
    dev.off()

    end_time <- Sys.time()
    print(difftime(end_time, start_time))
    print(paste("Start:", start_time, "End:", end_time))
}

get_y_label <- function(measure_type) {
    
    if (measure_type == "fscore") {
        plot_ylabs <- "F-score"
    }
    else if (measure_type == "precision") {
        plot_ylabs <- "Precision"
    }
    else if (measure_type == "recall") {
        plot_ylabs <- "Recall"
    }
    else if (measure_type == "all") {
        plot_ylabs <- c("F-score", "Precision", "Recall")
    }
    else {
        stop(paste("Unknown measure_type:", measure_type))
    }
    
    plot_ylabs    
}

get_norm_methods <- function(measure_type) {
    
    if (measure_type == "fscore") {
        score_funcs <- get_f_score
    }
    else if (measure_type == "precision") {
        score_funcs <- get_precision
    }
    else if (measure_type == "recall") {
        score_funcs <- get_recall
    }
    else if (measure_type == "all") {
        score_funcs <- list(get_f_score, get_precision, get_recall)
    }
    else {
        stop(paste("Unknown measure_type:", measure_type))
    }
    
    score_funcs
}

generate_normal_run_results <- function(nr, sig_thres=0.1, adjust_fdr=TRUE) {
    
    used_methods_names <- getUsedMethodNames(nr)
    normalization_matrices <- getNormalizationMatrices(nr)
    
    # header <- c("method", "tot_na_reduced", "rt_settings", "potato_tot", "potato_sig", "back_tot", "back_sig")
    # results_matrix <- data.frame(matrix(ncol=length(header), nrow=0))
    run_entries <- list()
    
    for (i in 1:length(used_methods_names)) {
        print(paste("Method: ", used_methods_names[[i]]))
        run_result <- get_run_entry(nr, normalization_matrices[[i]], used_methods_names[[i]], 
                                           POT_PAT, HUMAN_PAT, sig_thres=sig_thres, adjust_fdr=adjust_fdr)
        
        run_entries[[i]] <- run_result
    }
    run_entries
}

generate_rt_run_results <- function(nr, rt_windows, sig_thres=0.1, adjust_fdr=TRUE, frame_shifts=3) {
    
    used_methods_names <- getUsedMethodNames(nr)
    normalization_matrices <- getNormalizationMatrices(nr)
    
    normalyzer_filterraw <- nr@nds@filterrawdata
    retention_times <- nr@nds@retentionTimes

    run_entries <- list()
    for (rt_window in rt_windows) {

        print(paste("Method: RTs, window size:", rt_window))

        median_rt <- getSmoothedRTNormalizedMatrix(normalyzer_filterraw, retention_times, medianNormalization, rt_window, frame_shifts=frame_shifts)
        median_entry <- get_run_entry(nr, median_rt, "RT-median", POT_PAT, HUMAN_PAT, sig_thres=sig_thres, adjust_fdr=adjust_fdr, rt_settings=rt_window)

        mean_rt <- getSmoothedRTNormalizedMatrix(normalyzer_filterraw, retention_times, meanNormalization, rt_window, frame_shifts=frame_shifts)
        mean_entry <- get_run_entry(nr, mean_rt, "RT-mean", POT_PAT, HUMAN_PAT, sig_thres=sig_thres, adjust_fdr=adjust_fdr, rt_settings=rt_window)

        loess_rt <- getSmoothedRTNormalizedMatrix(normalyzer_filterraw, retention_times, performCyclicLoessNormalization, rt_window, frame_shifts=frame_shifts)
        loess_entry <- get_run_entry(nr, loess_rt, "RT-loess", POT_PAT, HUMAN_PAT, sig_thres=sig_thres, adjust_fdr=adjust_fdr, rt_settings=rt_window)

        run_entries_base_index <- length(run_entries) + 1
        run_entries[[run_entries_base_index]] <- median_entry
        run_entries[[run_entries_base_index + 1]] <- mean_entry
        run_entries[[run_entries_base_index + 2]] <- loess_entry
    }
    
    run_entries
}

# Generate entry object containing information on number of total and DE found in target pattern and background
get_run_entry <- function(nr, method_data, name, potato_pattern, human_pattern, sig_thres, adjust_fdr, rt_settings=NULL) {
    
    combined_pattern <- paste0("(", potato_pattern, "|", human_pattern, ")")
    row_name_col <- 3
    
    col_pattern <- SAMPLE_PAT

    prepared_df <- get_prepared_normalyzer_sheet(nr, method_data, col_pattern, row_name_col)
    na_filter_df <- na_filter(prepared_df)
    anova_pval <- get_anova_pvals(na_filter_df, custom_replicate_groups)
    
    if (adjust_fdr) {
        anova_fdr <- stats::p.adjust(anova_pval, method="BH")
        target_anova <- anova_fdr
    }
    else {
        target_anova <- anova_pval
    }
    
    potato_sig <- target_anova[which(grepl(potato_pattern, names(target_anova)))]
    nbr_potato_tot <- length(potato_sig)
    nbr_potato_sig <- length(nbr_potato_tot[which(potato_sig <= sig_thres)])
    
    back_sig <- target_anova[which(!grepl(combined_pattern, names(target_anova)))]
    nbr_back_tot <- length(back_sig)
    nbr_back_sig <- length(nbr_back_tot[which(back_sig <= sig_thres)])

    tot_passing_na_check <- nrow(na_filter_df)
    potato_entry <- EntryRow(norm_method=name, 
                             tot_rows=tot_passing_na_check,
                             target_tot=nbr_potato_tot,
                             target_sign=nbr_potato_sig,
                             background_tot=nbr_back_tot,
                             background_sign=nbr_back_sig,
                             rt_settings=rt_settings)
    potato_entry
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





