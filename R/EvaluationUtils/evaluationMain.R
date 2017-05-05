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
source("EvaluationUtils/evaluationUtils.R")
source("EvaluationUtils/outputUtils.R")
source("EvaluationUtils/evaluationStatistics.R")

suppressPackageStartupMessages("ggplot2")
library("ggplot2")
library("limma")


path_base <- "../tests/data"

target_dataset <- "pdx001819"

if (target_dataset == "max_quant_proteios") {
    full_report <- "FeatureReport_MQP_20140828-fornormalyzer.txt"
    subset_report <- "MQP.subset_500.tsv"
    NAME_COLUMNS <- c(3,2)
    SAMPLE_PAT_BASE <- "dilA_"
    SAMPLE_PAT_END <- "]_\\d\\d"
    
    POT_PAT <- "^sol"
    HUMAN_PAT <- "^hum"
    COMBINED_PATTERN <- paste0("(", POT_PAT, "|", HUMAN_PAT, ")")
    
    
} else if (target_dataset == "progenesis") {
    full_report <- "SpikeIn_Results_Progenesis_20130322_REFORMAT.NORMALYZER.txt"
    subset_report <- "progenesis.subset_500.tsv"
    NAME_COLUMNS <- c(2,1)
    SAMPLE_PAT_BASE <- "dilA_"
    SAMPLE_PAT_END <- "]_\\d\\d"
    
    POT_PAT <- "^sol"
    HUMAN_PAT <- "^hum"
    COMBINED_PATTERN <- paste0("(", POT_PAT, "|", HUMAN_PAT, ")")
    
} else if (target_dataset == "dinosaur_test") {
    full_report <- "FeatureReport_dino_nulltona.tsv"
    subset_report <- "FeatureReport_dino_500.tsv"
    NAME_COLUMNS <- c(4,3)
    SAMPLE_PAT_BASE <- "dil"
    SAMPLE_PAT_END <- "]"
    
    POT_PAT <- "^sol"
    HUMAN_PAT <- "^hum"
    COMBINED_PATTERN <- paste0("(", POT_PAT, "|", HUMAN_PAT, ")")
    
} else if (target_dataset == "dinosaur_massmatch") {
    full_report <- "eval_dinosaurs/FeatureReport_dino_massmatch.tsv"
    NAME_COLUMNS <- c(4,3)
    SAMPLE_PAT_BASE <- "130124_dil\\w_"
    SAMPLE_PAT_END <- "]_\\d\\d.mzML"
    
    POT_PAT <- "^sol"
    HUMAN_PAT <- "^hum"
    COMBINED_PATTERN <- paste0("(", POT_PAT, "|", HUMAN_PAT, ")")
    
} else if (target_dataset == "dinosaur_theoretical") {
    full_report <- "eval_dinosaurs/FeatureReport_dino_theoretical.tsv"
    NAME_COLUMNS <- c(4,3)
    SAMPLE_PAT_BASE <- "130124_dil\\w_"
    SAMPLE_PAT_END <- "]_\\d\\d.mzML"
    
} else if (target_dataset == "max_quant_proteios_background") {
    full_report <- "FeatureReport_MQP_20140828-fornormalyzer_without_potato.txt"
    NAME_COLUMNS <- c(3,2)
    SAMPLE_PAT_BASE <- "dil[ABCD]_"
    SAMPLE_PAT_END <- "]_\\d\\d"
    
    POT_PAT <- "^sol"
    HUMAN_PAT <- "^hum"
    COMBINED_PATTERN <- paste0("(", POT_PAT, "|", HUMAN_PAT, ")")
  
} else if (target_dataset == "pdx001819") {
  
    full_report <- "ups_data/pdx001819_report2_features_only_normalyzer_formatted_na.tsv"
    subset_report <- "ups_data/pdx001819_report2_features_only_normalyzer_formatted_na.500.tsv"
    NAME_COLUMNS <- c(4,3)
    SAMPLE_PAT_BASE <- "s_"
    SAMPLE_PAT_END <- "]amol"
    
    POT_PAT <- "^\\w+ups"
    HUMAN_PAT <- "^hum"
    COMBINED_PATTERN <- "^\\w+ups"
      
} else {
    stop(paste("Unknown dataset", target_dataset))
}
# subset_report <- "BA_MQP.subset_500.tsv"


NAME_FILTER <- "[, ]"
# ANNOTATION_FILTER <- "[,| ]"


run_review_data_test <- function(super_dirname,
                                 run_dirname,
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
                                 var_filter_frac=NULL,
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
                     var_filter_frac=var_filter_frac,
                     sample_comp=target_replicates,
                     stat_test=stat_test
                     )

    rs$run_dir <- run_dirname
    
    print(rs$run_dir)
    
    output_base <- get_run_setting_base(rs)
    output_path <- paste(output_base, "csv", sep=".")
    plot_path <- paste(output_base, "pdf", sep=".")
    print(paste("Evaluation run with sig_thres:", rs$sig_thres, "and FDR setting:", rs$do_fdr))
    plot_ylabs <- get_y_label("all")
    score_funcs <- get_norm_methods("all")

    start_time <- Sys.time()
    
    if (subset) review_data_path <- paste(path_base, subset_report, sep="/")
    else review_data_path <- paste(path_base, full_report, sep="/")
    
    print(paste("------ START: PROCESSING DATASET IN PATH", review_data_path))
    print(paste("Running replicates: ", paste(target_replicates, collapse=",")))
    
    # print(paste("Processing dataset:", review_data_path))
    job_name <- "review_evaluation"
    norm_obj <- getVerifiedNormalyzerObject(review_data_path, job_name)

    # browser()
    
    nr <- normMethods(norm_obj, job_name, normalizeRetentionTime=FALSE, runNormfinder=FALSE)
    
    if (do_normal_run) {
        normal_run_entries <- generate_normal_run_results(nr, rs, row_name_cols=NAME_COLUMNS)
    }
        
    if (do_rt_run) {
        rt_run_entries <- generate_rt_run_results(nr, rs, name_col=3, row_name_cols=NAME_COLUMNS)
    }
    
    # Postfiltering to match and remove lines only present in some matrices
    
    print(paste("Normal entries:", length(normal_run_entries)))
    if (do_rt_run) {
        print(paste("RT entries:", length(rt_run_entries)))
    }
    
    if (verbose) {
        print("--- Normal entries ---")
        for (entry in normal_run_entries) {
            print(get_entry_string(entry))
        }
        
        if (do_rt_run) {
            print("--- RT entries ---")
            for (entry in rt_run_entries) {
                print(get_entry_string(entry))
            }
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
        
        if (!do_rt_run) {
            rt_run_entries <- c()
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
        
        if (do_rt_run) {
            write_entries_to_file(rt_sub_path, rt_run_entries)
        }
    }
    
    write_data_matrices(output_base, normal_run_entries, rt_run_entries, nr)

    pdf(plot_path)
    grid_arrange_shared_legend(plots=plots, ncol=length(plots), plot_info=gsub(".*\\/", "", output_base))
    dev.off()

    sign_out_path <- paste(output_base, "sign_matrix", "tsv", sep=".")
    print(paste("Writing significance entries to:", sign_out_path))
    
    
    # write_significance_matrix(sign_out_path, c(normal_run_entries, rt_run_entries))
    
    end_time <- Sys.time()
    print(difftime(end_time, start_time))
    print(paste("Start:", start_time, "End:", end_time))
}


generate_normal_run_results <- function(nr, rs, row_name_cols) {

    used_methods_names <- getUsedMethodNames(nr)
    normalization_matrices <- getNormalizationMatrices(nr)
    
    run_entries <- list()
    for (i in 1:length(used_methods_names)) {
        print(paste("Method: ", used_methods_names[[i]]))
        
        run_result <- get_run_entry(nr, normalization_matrices[[i]], used_methods_names[[i]], 
                                    POT_PAT, rs, row_name_cols=row_name_cols)
        
        run_entries[[i]] <- run_result
    }
    run_entries
}

generate_rt_run_results <- function(nr, rs, row_name_cols, name_col=NULL) {

    used_methods_names <- getUsedMethodNames(nr)
    normalization_matrices <- getNormalizationMatrices(nr)
    
    normalyzer_filterraw <- nr@nds@filterrawdata
    retention_times <- nr@nds@retentionTimes

    if (!is.null(name_col)) {
        rownames(normalyzer_filterraw) <- nr@nds@rawData[-1:-2, name_col]
    }
    
    row_names <- get_rownames(nr, NAME_COLUMNS)
    col_names <- get_colnames(nr)
    
    run_entries <- list()
    for (rt_window in rs$rt_windows) {

        frame_shifts <- rs$window_shifts
        win_size_min <- rs$lowest_window_size
        merge_method <- rs$window_merge_method
        
        raw_copy <- nr@data2log2

        rownames(raw_copy) <- row_names
        colnames(raw_copy) <- col_names[2:length(col_names)]
        
        print(paste("Current rt setting:", rt_window))
        
        print("median rt")
        median_rt <- getSmoothedRTNormalizedMatrix(normalyzer_filterraw, 
                                                   retention_times, 
                                                   medianNormalization, 
                                                   rt_window, 
                                                   frame_shifts=frame_shifts,
                                                   win_size_min=win_size_min,
                                                   merge_method=merge_method,
                                                   debug_matrix_header=col_names,
                                                   debug_matrix_rownames=row_names,
                                                   debug_id=debug_id,
                                                   debug_samples=debug_samples)
        
        median_entry <- get_run_entry(nr, median_rt, "RT-median", POT_PAT, rs, current_rt_setting=rt_window, row_name_cols=row_name_cols)

        print("mean rt")
        mean_rt <- getSmoothedRTNormalizedMatrix(normalyzer_filterraw,
                                                 retention_times,
                                                 meanNormalization,
                                                 rt_window,
                                                 frame_shifts=frame_shifts,
                                                 win_size_min=win_size_min,
                                                 merge_method=merge_method,
                                                 debug_matrix_header=col_names,
                                                 debug_matrix_rownames=row_names,
                                                 debug_id=debug_id,
                                                 debug_samples=debug_samples)
        mean_entry <- get_run_entry(nr, mean_rt, "RT-mean", POT_PAT, rs, current_rt_setting=rt_window, row_name_cols=row_name_cols)

        print("loess rt")
        loess_rt <- getSmoothedRTNormalizedMatrix(normalyzer_filterraw,
                                                  retention_times,
                                                  performCyclicLoessNormalization,
                                                  rt_window,
                                                  frame_shifts=frame_shifts,
                                                  win_size_min=win_size_min,
                                                  merge_method=merge_method,
                                                  debug_matrix_header=col_names,
                                                  debug_matrix_rownames=row_names,
                                                  debug_id=debug_id,
                                                  debug_samples=debug_samples)
        loess_entry <- get_run_entry(nr, loess_rt, "RT-loess", POT_PAT, rs, current_rt_setting=rt_window, row_name_cols=row_name_cols)

        run_entries_base_index <- length(run_entries) + 1
        run_entries[[run_entries_base_index]] <- median_entry
        run_entries[[run_entries_base_index + 1]] <- mean_entry
        run_entries[[run_entries_base_index + 2]] <- loess_entry
    }
    
    run_entries
}

# Generate entry object containing information on number of total and DE found in target pattern and background
get_run_entry <- function(nr, method_data, name, target_pattern, rs, row_name_cols, current_rt_setting=NULL, omit_na=TRUE) {

    sig_thres <- rs$sig_thres
    adjust_fdr <- rs$do_fdr
    replicate_nbrs <- rs$sample_comp
    stat_test <- rs$stat_test
    var_filter_frac <- rs$var_filter_frac

    # combined_pattern <- paste0("(", target_pattern, "|", human_pattern, ")")
    # row_name_col <- c(2,3)
    # col_pattern <- paste0(SAMPLE_PAT_BASE, "[", paste(replicate_nbrs, sep="", collapse=""), SAMPLE_PAT_END)
    
    all_replicates <- nr@nds@sampleReplicateGroups
    custom_replicate_groups <- all_replicates[which(all_replicates %in% replicate_nbrs)]

    # browser()
    
    prepared_df <- get_prepared_normalyzer_sheet(nr, method_data, row_name_cols, replicate_nbrs, var_filter_frac=var_filter_frac)
    filter_df <- prepared_df[which(!grepl(NAME_FILTER, rownames(prepared_df))),]

    if (stat_test == "welch" || stat_test == "student") {
        target_sig_val <- get_anova_pvals(filter_df, custom_replicate_groups, ttest_type=stat_test, adjust_fdr=T)
        print(paste("Running regular test: ", stat_test))
    }
    else if (stat_test == "limma") {
        target_sig_val <- get_limma_pvals(filter_df, custom_replicate_groups, adjust_fdr=T)
        print("Running Limma test")
    }
    else {
        stop(paste("Unknown stat_test:", stat_test))
    }
    
    potato_sig <- target_sig_val[which(grepl(target_pattern, names(target_sig_val)))]
    back_sig <- target_sig_val[which(!grepl(target_pattern, names(target_sig_val)))]
    
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
    
    potato_entry <- EntryRow(nr=nr,
                             method_data=filter_df,
                             norm_method=name, 
                             tot_rows=total_rows,
                             target_tot=nbr_potato_tot,
                             target_sign=nbr_potato_sig,
                             background_tot=nbr_back_tot,
                             background_sign=nbr_back_sig,
                             rt_settings=current_rt_setting,    
                             sig_df=target_sig_val)
    potato_entry
}


get_prepared_normalyzer_sheet <- function(nr, df, row_name_cols, replicate_nbrs, var_filter_frac) {

    raw_data <- nr@nds@rawData
    replicate_numbers <- raw_data[1,which(as.integer(raw_data[1,]) > 0)]
    replicate_pattern <- paste0("^(", replicate_nbrs[1], "|", replicate_nbrs[2], ")$")
    
    names_cols <- raw_data[-1:-2, row_name_cols]
    row_names <- paste(names_cols[,1], names_cols[,2], sep="|")
    
    target_cols <- which(grepl(replicate_pattern, replicate_numbers))
    data_df <- df[, target_cols]

    all_replicates <- nr@nds@sampleReplicateGroups
    target_replicates <- all_replicates[which(all_replicates %in% replicate_nbrs)]
    
    rownames(data_df) <- row_names
    colnames(data_df) <- target_replicates
    
    filter_contrast <- getRowNAFilterContrast(data_df, target_replicates, var_filter_frac)
    filtered_df <- data_df[filter_contrast,]
    
    filtered_df
}


















