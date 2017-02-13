

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
# source("EvaluationUtils/RunSetting.R")

suppressPackageStartupMessages("ggplot2")
library("ggplot2")


path_base <- "../tests/data"
full_report <- "FeatureReport_MQP_20140828-fornormalyzer.txt"
subset_report <- "MQP.subset_500.tsv"
# subset_report <- "BA_MQP.subset_500.tsv"

POT_PAT <- "^sol"
HUMAN_PAT <- "^hum"
SAMPLE_PAT_BASE <- "dilA_"
SAMPLE_PAT_END <- "]_\\d\\d"
# SAMPLE_PAT_BASE <- "dilA_[23]"
COMBINED_PATTERN <- paste0("(", POT_PAT, "|", HUMAN_PAT, ")")
# custom_replicate_groups <- c(2,2,2,2,2,2,2,3,3,3,3,3,3,3)


run_review_data_test <- function(super_dirname,
                                 subset=FALSE,
                                 do_normal_run=TRUE,
                                 do_rt_run=TRUE,
                                 verbose=FALSE,
                                 measure_type="all",
                                 sig_thres=0.1,
                                 do_fdr=TRUE,
                                 frame_shifts=3,
                                 lowest_window_size=1,
                                 window_merge_method="median",
                                 custom_rt_windows=NULL,
                                 target_replicates=c(1,3),
                                 stat_test="anova") {

    rs <- RunSetting(sig_thres=sig_thres,
                     do_fdr=do_fdr,
                     rt_windows=custom_rt_windows,
                     window_shifts=frame_shifts,
                     lowest_window_size=lowest_window_size,
                     window_merge_method=window_merge_method,
                     quiet=!verbose,
                     subset=subset,
                     super_dirname=super_dirname,
                     sample_comp=target_replicates,
                     stat_test=stat_test
                     )
    
    output_base <- get_run_setting_base(rs)
    output_path <- paste(output_base, "csv", sep=".")
    plot_path <- paste(output_base, "png", sep=".")
    print(paste("Evaluation run with sig_thres:", rs$sig_thres, "and FDR setting:", rs$do_fdr))
    plot_ylabs <- get_y_label("all")
    score_funcs <- get_norm_methods("all")

    start_time <- Sys.time()
    
    if (subset) review_data_path <- paste(path_base, subset_report, sep="/")
    else review_data_path <- paste(path_base, full_report, sep="/")
    
    print(paste("Processing dataset:", review_data_path))
    job_name <- "review_evaluation"
    norm_obj <- getVerifiedNormalyzerObject(review_data_path, job_name)

    if (!is.null(custom_rt_windows))    rt_windows <- custom_rt_windows
    else if (subset)                    rt_windows <- c(1, 2, 3)
    else                                rt_windows <- c(0.2, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 7, 8, 9, 10, 13, 16, 20)
    
    nr <- normMethods(norm_obj, job_name, normalizeRetentionTime=FALSE, runNormfinder=FALSE)
    
    if (do_normal_run) {
        normal_run_entries <- generate_normal_run_results(nr, rs)
    }
        
    if (do_rt_run) {
        rt_run_entries <- generate_rt_run_results(nr, rs)
    }
    
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
        
        include_legend <- TRUE
        if (i != length(score_funcs)) {
            include_legend <- FALSE
        }
        
        plt <- visualize_results_new(normal_run_entries,
                              rt_run_entries,
                              score_func,
                              title=y_lab,
                              xlab="RT window size (min)",
                              ylab=y_lab,
                              verbose=verbose)
        plots[[i]] <- plt
        
        norm_sub_path <- paste(plot_path, tolower(y_lab), "norm", "csv", sep=".")
        rt_sub_path <- paste(plot_path, tolower(y_lab), "rt", "csv", sep=".")
        write_entries_to_file(norm_sub_path, normal_run_entries)
        write_entries_to_file(rt_sub_path, rt_run_entries)
    }

    # plot_path <- paste(plot_path, "png", sep=".")
    png(plot_path)
    # multiplot(plotlist=plots, cols=length(plots))
    grid_arrange_shared_legend(plots=plots, ncol=length(plots))
    dev.off()

    sign_out_path <- paste(plot_path, "sign_matrix", "tsv", sep=".")
    print(paste("Writing significance entries to:", sign_out_path))
    write_significance_matrix(sign_out_path, c(normal_run_entries, rt_run_entries))
    
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

generate_normal_run_results <- function(nr, rs) {

    used_methods_names <- getUsedMethodNames(nr)
    normalization_matrices <- getNormalizationMatrices(nr)
    
    run_entries <- list()
    for (i in 1:length(used_methods_names)) {
        print(paste("Method: ", used_methods_names[[i]]))
        run_result <- get_run_entry(nr, normalization_matrices[[i]], used_methods_names[[i]], 
                                    POT_PAT, HUMAN_PAT, rs)

        run_entries[[i]] <- run_result
    }
    run_entries
}

generate_rt_run_results <- function(nr, rs) {

    used_methods_names <- getUsedMethodNames(nr)
    normalization_matrices <- getNormalizationMatrices(nr)
    
    normalyzer_filterraw <- nr@nds@filterrawdata
    retention_times <- nr@nds@retentionTimes

    run_entries <- list()
    for (rt_window in rs$rt_windows) {

        print(paste("Method: RTs, window size:", rt_window))

        frame_shifts <- rs$window_shifts
        win_size_min <- rs$lowest_window_size
        merge_method <- rs$window_merge_method
        
        median_rt <- getSmoothedRTNormalizedMatrix(normalyzer_filterraw, 
                                                   retention_times, 
                                                   medianNormalization, 
                                                   rt_window, 
                                                   frame_shifts=frame_shifts,
                                                   win_size_min=win_size_min,
                                                   merge_method=merge_method)
        median_entry <- get_run_entry(nr, median_rt, "RT-median", POT_PAT, HUMAN_PAT, rs, current_rt_setting=rt_window)

        mean_rt <- getSmoothedRTNormalizedMatrix(normalyzer_filterraw,
                                                 retention_times,
                                                 meanNormalization,
                                                 rt_window,
                                                 frame_shifts=frame_shifts,
                                                 win_size_min=win_size_min,
                                                 merge_method=merge_method)
        mean_entry <- get_run_entry(nr, mean_rt, "RT-mean", POT_PAT, HUMAN_PAT, rs, current_rt_setting=rt_window)

        loess_rt <- getSmoothedRTNormalizedMatrix(normalyzer_filterraw,
                                                  retention_times,
                                                  performCyclicLoessNormalization,
                                                  rt_window,
                                                  frame_shifts=frame_shifts,
                                                  win_size_min=win_size_min,
                                                  merge_method=merge_method)
        loess_entry <- get_run_entry(nr, loess_rt, "RT-loess", POT_PAT, HUMAN_PAT, rs, current_rt_setting=rt_window)

        run_entries_base_index <- length(run_entries) + 1
        run_entries[[run_entries_base_index]] <- median_entry
        run_entries[[run_entries_base_index + 1]] <- mean_entry
        run_entries[[run_entries_base_index + 2]] <- loess_entry
    }
    
    run_entries
}

# Generate entry object containing information on number of total and DE found in target pattern and background
get_run_entry <- function(nr, method_data, name, potato_pattern, human_pattern, rs, current_rt_setting=NULL) {

    sig_thres <- rs$sig_thres
    adjust_fdr <- rs$do_fdr
    replicate_nbrs <- rs$sample_comp
    stat_test <- rs$stat_test

    combined_pattern <- paste0("(", potato_pattern, "|", human_pattern, ")")
    row_name_col <- 3
    
    col_pattern <- paste0(SAMPLE_PAT_BASE, "[", paste(replicate_nbrs, sep="", collapse=""), SAMPLE_PAT_END)

    all_replicates <- nr@nds@sampleReplicateGroups
    custom_replicate_groups <- all_replicates[which(all_replicates %in% replicate_nbrs)]
    
    prepared_df <- get_prepared_normalyzer_sheet(nr, method_data, col_pattern, row_name_col, replicate_nbrs)
    na_filter_df <- na_filter(prepared_df)
    anova_pval <- get_anova_pvals(na_filter_df, custom_replicate_groups, stat_test=stat_test)
    
    if (adjust_fdr) {
        anova_fdr <- stats::p.adjust(anova_pval, method="BH")
        target_sig_val <- anova_fdr
    }
    else {
        target_sig_val <- anova_pval
    }

    # str(target_sig_val)
    # print(head(target_sig_val))
    # stop("")
        
    potato_sig <- target_sig_val[which(grepl(potato_pattern, names(target_sig_val)))]
    nbr_potato_tot <- length(potato_sig)
    nbr_potato_sig <- length(nbr_potato_tot[which(potato_sig <= sig_thres)])
    
    back_sig <- target_sig_val[which(!grepl(combined_pattern, names(target_sig_val)))]
    nbr_back_tot <- length(back_sig)
    nbr_back_sig <- length(nbr_back_tot[which(back_sig <= sig_thres)])

    tot_passing_na_check <- nrow(na_filter_df)
    potato_entry <- EntryRow(norm_method=name, 
                             tot_rows=tot_passing_na_check,
                             target_tot=nbr_potato_tot,
                             target_sign=nbr_potato_sig,
                             background_tot=nbr_back_tot,
                             background_sign=nbr_back_sig,
                             rt_settings=current_rt_setting,
                             sig_df=target_sig_val)
    potato_entry
}


get_prepared_normalyzer_sheet <- function(nr, df, col_pattern, row_name_col, replicate_nbrs) {

    raw_data <- nr@nds@rawData
    header_row <- raw_data[2,]
    names <- raw_data[-1:-2, row_name_col] 

    target_cols <- which(grepl(col_pattern, header_row))
    parsed_df <- df[, target_cols]

    all_replicates <- nr@nds@sampleReplicateGroups
    target_replicates <- all_replicates[which(all_replicates %in% replicate_nbrs)]

    rownames(parsed_df) <- names
    colnames(parsed_df) <- target_replicates
    
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


na_filter <- function(df, max_na=2) {
    na_per_line <- rowSums(is.na(df))
    df_na_removed <- df[na_per_line <= max_na, ]
    # df_na_removed <- df[na_per_line < ncol(df) / 2, ]
    df_na_removed
}

get_anova_pvals <- function(df, replicate_groups, stat_test="anova") {

    # Is this needed?    
    replicate_groups <- as.factor(replicate_groups)
    
    if (stat_test == "anova") {
        anova_func <- function(sampleIndex) {
            t.test(unlist(sampleIndex)~replicate_groups, var.equal=TRUE)$p.value
        }
    }
    else if (stat_test == "welch") {
        anova_func <- function(sampleIndex) {
            t.test(unlist(sampleIndex)~replicate_groups)$p.value
        }
    }
    else {
        stop(paste("Unknown test type: ", stat_test))
    }
    
    anovaPVal <- apply(df, 1, anova_func)
    anovaPVal
}

write_entries_to_file <- function(out_path, run_entries) {
    
    first_entry_v <- get_entry_vector(run_entries[[1]])
    out_strs <- matrix(ncol=length(first_entry_v), nrow=0)
    first_entry_head <- get_entry_header(run_entries[[1]])
    colnames(out_strs) <- get_entry_header(run_entries[[1]])
    
    for (i in 1:length(run_entries)) {
        
        entry <- run_entries[[i]]
        out_strs <- rbind(out_strs, get_entry_vector(entry))
    }
    
    write.csv(out_strs, file=out_path, quote=FALSE, row.names=FALSE)
}

write_significance_matrix <- function(sign_out_path, entries) {
    
    nbr_features <- entries[[1]]$tot_rows
    nbr_methods <- length(entries)
    
    sign_matrix <- data.frame(matrix(nrow=nbr_features, ncol=0))
    entry_names <- NULL
    
    for (i in 1:nbr_methods) {
        
        entry <- entries[[i]]
        
        if (is.null(entry_names)) {
            entry_names <- names(entry$sig_df)
        }
        
        entry_col <- entry$sig_df

        if (!is.null(entry$rt_settings)) {
            name <- paste(entry$norm_method, entry$rt_settings, sep="_")
        }
        else {
            name <- entry$norm_method
        }
        
        sign_matrix[name] <- entry_col
    }

    feature_names <- entry_names
    unique_rownames <- paste(feature_names, seq(1, entries[[1]]$tot_rows), sep="_")
    rownames(sign_matrix) <- unique_rownames
      
    write.table(sign_matrix, file=sign_out_path, sep="\t", quote=FALSE, col.names=NA)
}



















