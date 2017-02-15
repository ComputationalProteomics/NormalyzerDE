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
source("EvaluationUtils/RunSetting.R")

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
COMBINED_PATTERN <- paste0("(", POT_PAT, "|", HUMAN_PAT, ")")

NAME_FILTER <- "[,]"


run_review_data_test <- function(super_dirname,
                                 subset=FALSE,
                                 do_normal_run=TRUE,
                                 do_rt_run=TRUE,
                                 verbose=FALSE,
                                 measure_type="all",
                                 sig_thres=0.1,
                                 do_fdr=TRUE,
                                 frame_shifts=3,
                                 lowest_window_size=50,
                                 window_merge_method="median",
                                 rt_windows=NULL,
                                 target_replicates=c(1,3),
                                 stat_test="welch") {

    rs <- RunSetting(sig_thres=sig_thres,
                     do_fdr=do_fdr,
                     rt_windows=rt_windows,
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

    nr <- normMethods(norm_obj, job_name, normalizeRetentionTime=FALSE, runNormfinder=FALSE)
    
    if (do_normal_run) {
        normal_run_entries <- generate_normal_run_results(nr, rs)
    }
        
    if (do_rt_run) {
        rt_run_entries <- generate_rt_run_results(nr, rs)
    }
    
    # Postfiltering to match and remove lines only present in some matrices
    
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
        
        norm_sub_path <- paste(output_base, tolower(y_lab), "norm", "csv", sep=".")
        rt_sub_path <- paste(output_base, tolower(y_lab), "rt", "csv", sep=".")
        write_entries_to_file(norm_sub_path, normal_run_entries)
        write_entries_to_file(rt_sub_path, rt_run_entries)
    }

    # plot_path <- paste(plot_path, "png", sep=".")
    png(plot_path)
    # multiplot(plotlist=plots, cols=length(plots))
    grid_arrange_shared_legend(plots=plots, ncol=length(plots))
    dev.off()

    sign_out_path <- paste(output_base, "sign_matrix", "tsv", sep=".")
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

        # print(paste("Method: RTs, window size:", rt_window))

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
get_run_entry <- function(nr, method_data, name, potato_pattern, human_pattern, rs, current_rt_setting=NULL, omit_na=TRUE) {

    sig_thres <- rs$sig_thres
    adjust_fdr <- rs$do_fdr
    replicate_nbrs <- rs$sample_comp
    stat_test <- rs$stat_test

    # print(paste("sig:", sig_thres, "fdr:", adjust_fdr, "repl", replicate_nbrs, "stat test", stat_test))
    
    combined_pattern <- paste0("(", potato_pattern, "|", human_pattern, ")")
    row_name_col <- 3
    
    col_pattern <- paste0(SAMPLE_PAT_BASE, "[", paste(replicate_nbrs, sep="", collapse=""), SAMPLE_PAT_END)

    all_replicates <- nr@nds@sampleReplicateGroups
    custom_replicate_groups <- all_replicates[which(all_replicates %in% replicate_nbrs)]

    prepared_df <- get_prepared_normalyzer_sheet(nr, method_data, col_pattern, row_name_col, replicate_nbrs)
    
    #### Skip this, and let test return NA instead if invalid?
    # print(head(prepared_df))
    # print(head(rownames(prepared_df)))
    
    filter_df <- prepared_df[which(!grepl(NAME_FILTER, rownames(prepared_df))),]
    
    # na_filter_df <- na_filter(prepared_df)
    # print(na_filter_df)
    # print(paste("Length before:", length(prepared_df[,1])))
    # print(paste("Length after:", length(na_filter_df[,1])))
    # stop("")
    # str(na_filter_df)
    
    anova_pval <- get_anova_pvals(filter_df, custom_replicate_groups, stat_test=stat_test)
    
    if (adjust_fdr) {
        anova_fdr <- stats::p.adjust(anova_pval, method="BH")
        target_sig_val <- anova_fdr
    }
    else {
        target_sig_val <- anova_pval
    }

    potato_sig <- target_sig_val[which(grepl(potato_pattern, names(target_sig_val)))]
    back_sig <- target_sig_val[which(!grepl(combined_pattern, names(target_sig_val)))]
    
    if (omit_na) {
        nbr_potato_tot <- length(na.omit(potato_sig))
        nbr_potato_sig <- length(na.omit(potato_sig[which(potato_sig <= sig_thres)]))
        nbr_back_tot <- length(na.omit(back_sig))
        nbr_back_sig <- length(na.omit(back_sig[which(back_sig <= sig_thres)]))
    }
    else {
        nbr_potato_tot <- length(potato_sig)
        nbr_potato_sig <- length(potato_sig[which(potato_sig <= sig_thres)])
        nbr_back_tot <- length(back_sig)
        nbr_back_sig <- length(back_sig[which(back_sig <= sig_thres)])
    }

    total_rows <- length(na.omit(potato_sig)) + length(na.omit(back_sig))
    
    potato_entry <- EntryRow(norm_method=name, 
                             tot_rows=total_rows,
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


# Filter out lines with high NA abundance - based on log matrix for comparability
# na_filter <- function(df, max_na=2) {
#     
#     # log2_data <- nr@data2log2
#     
#     # na_per_line <- rowSums(is.na(log2_data))
#     na_per_line <- rowSums(is.na(df))
#     df_na_removed <- df[na_per_line <= max_na, ]
#     # df_na_removed <- df[na_per_line < ncol(df) / 2, ]
#     df_na_removed
# }

get_anova_pvals <- function(df, replicate_groups, stat_test="anova") {

    # Is this needed?    
    replicate_groups <- as.factor(replicate_groups)
    
    stat_test_func <- function(sampleIndex) {
        
        rep_counts <- table(names(na.omit(sampleIndex)))
        
        if (length(rep_counts) == 2 && min(rep_counts > 1)) {
            if (stat_test == "anova") {
                t.test(unlist(sampleIndex)~replicate_groups, var.equal=TRUE)$p.value
            }
            else if (stat_test == "welch") {
                t.test(unlist(sampleIndex)~replicate_groups)$p.value
            }
            else {
                stop(paste("Unknown test type:", stat_test))
            }
        }
        else {
            NA
        }
    }
    
    anovaPVal <- apply(df, 1, stat_test_func)
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
    
    # nbr_features <- entries[[1]]$target_tot + entries[[1]]$background_tot
    nbr_features <- length(entries[[1]]$sig_df)
    nbr_methods <- length(entries)
    
    # print(paste("nbr_features", nbr_features))
    # stop("")
    
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
        
        print(paste(name, length(entry_col)))
        sign_matrix[name] <- entry_col
    }

    feature_names <- entry_names
    unique_rownames <- paste(feature_names, seq(1, nbr_features), sep="_")
    
    n_t <- table(unique_rownames)
    print(n_t[n_t > 1])
    print(length(n_t))
    
    rownames(sign_matrix) <- unique_rownames
      
    write.table(sign_matrix, file=sign_out_path, sep="\t", quote=FALSE, col.names=NA)
}



















