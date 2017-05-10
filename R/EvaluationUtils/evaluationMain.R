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

# target_dataset <- "pdx001819"
target_dataset <- "max_quant_proteios"

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
                     stat_test=stat_test,
                     run_dir=run_dirname
                     )

    output_base <- get_run_setting_base(rs)
    output_path <- paste(output_base, "csv", sep=".")
    print(paste("Evaluation run with sig_thres:", rs$sig_thres, "and FDR setting:", rs$do_fdr))

    start_time <- Sys.time()
    
    if (subset) review_data_path <- paste(path_base, subset_report, sep="/")
    else review_data_path <- paste(path_base, full_report, sep="/")
    
    print(paste("------ START: PROCESSING DATASET IN PATH", review_data_path))
    print(paste("Running replicates: ", paste(target_replicates, collapse=",")))
    
    job_name <- "review_evaluation"
    norm_obj <- getVerifiedNormalyzerObject(review_data_path, job_name)

    nr <- normMethods(norm_obj, job_name, normalizeRetentionTime=FALSE, runNormfinder=FALSE)
    
    if (do_normal_run) {
        normal_run_entries <- generate_normal_run_results(nr, rs, row_name_cols=NAME_COLUMNS)
    }

    if (do_rt_run) {
        rt_run_entries <- generate_rt_run_results(nr, rs, name_col=3, row_name_cols=NAME_COLUMNS)
    }
    
    print(paste("Normal entries:", length(normal_run_entries)))
    if (do_rt_run) {
        print(paste("RT entries:", length(rt_run_entries)))
    }
    
    roc_list <- list()
    roc_list[["roc_plot_normal"]] <- generate_roc_plot(normal_run_entries, rt_run_entries, only_normal=T)
    roc_list[["roc_plot_normal_04cut"]] <- generate_roc_plot(normal_run_entries, rt_run_entries, only_normal=T, y_cutoff=0.4)
    roc_list[["roc_plot_normal_06cut"]] <- generate_roc_plot(normal_run_entries, rt_run_entries, only_normal=T, y_cutoff=0.6)
    roc_list[["roc_plot_normal_08cut"]] <- generate_roc_plot(normal_run_entries, rt_run_entries, only_normal=T, y_cutoff=0.8)
    roc_list[["roc_plot_all"]] <- generate_roc_plot(normal_run_entries, rt_run_entries, only_normal=F)
    score_plots <- generate_score_plots(normal_run_entries, rt_run_entries, output_base)
    # generate_results(normal_run_entries, rt_run_entries, output_base, score_plots, roc_plot_normal, roc_plot_all)
    generate_results(normal_run_entries, rt_run_entries, output_base, score_plots, roc_list)
    
    end_time <- Sys.time()
    print(difftime(end_time, start_time))
    print(paste("Start:", start_time, "End:", end_time))
}


generate_roc_plot <- function(normal_run_entries, rt_run_entries, use_fdr=T, only_normal=T, debug=T, sig_coord_thres=0.1, y_cutoff=0) {
    
    comb_df <- data.frame(pvals=c(), type=c(), method=c())
    sig_coord_df <- data.frame(x=c(), y=c(), method=c())
    
    if (!only_normal) {
        all_entries <- c(normal_run_entries, rt_run_entries)
    }
    else {
        all_entries <- normal_run_entries
    }
    
    plot_ylabs <- get_y_label("all")
    
    for (method_i in 1:length(all_entries)) {
        
        target_obj <- all_entries[[method_i]]
        
        if (use_fdr) {
            target_sig <- all_entries[[method_i]][["diff_sig_fdr"]]
            background_sig <- all_entries[[method_i]][["back_sig_fdr"]]
        }
        else {
            target_sig <- all_entries[[method_i]][["diff_sig_p"]]
            background_sig <- all_entries[[method_i]][["back_sig_p"]]
        }

        if (debug) {
            write.csv(target_sig, file="/home/jakob/Desktop/debug/diff.csv")
            write.csv(background_sig, file="/home/jakob/Desktop/debug/back.csv")
        }

        found_pos <- 0
        found_neg <- 0
        last_x_coord<- 0
        last_y_coord <- 0
        x_vals <- c(0)
        y_vals <- c(0)
        
        sig_coord <- NULL
        
        tot_sig <- length(target_sig)
        tot_back <- length(background_sig)
                
        target_df <- cbind(pvals=target_sig, type=1, sample=method_i)
        background_df <- cbind(pvals=background_sig, type=0, sample=method_i)
        
        unsorted_df <- rbind(target_df, background_df)
        rownames(unsorted_df) <- seq(1, nrow(unsorted_df))
        sorted_df <- unsorted_df[order(unsorted_df[, "pvals"]),]

        for (i in 1:nrow(sorted_df)) {

            row <- sorted_df[i,]
            if (row["type"] == 0) {
                found_neg <- found_neg + 1
                last_x_coord <- found_neg / tot_back
            }
            else {
                found_pos <- found_pos + 1
                last_y_coord <- found_pos / tot_sig
            }
            x_vals <- c(x_vals, last_x_coord)
            y_vals <- c(y_vals, last_y_coord)
            
            if (is.null(sig_coord) && row["pvals"] > 0.1) {
                sig_coord <- data.frame(x=last_x_coord, y=last_y_coord, method=target_obj$norm_method)
                sig_coord_df <- rbind(sig_coord_df, sig_coord)
            }
        }
        
        coord_df <- data.frame(FPR=x_vals, TPR=y_vals, sample=target_obj$norm_method)
        comb_df <- rbind(comb_df, coord_df)
    }

    plt <- ggplot() + 
        geom_line(data=comb_df, aes(FPR, TPR, color=sample)) + 
        ggtitle("ROC curve") +
        geom_point(data=sig_coord_df, aes(x=x, y=y, color=method)) +
        ylim(y_cutoff, NA)
    
    if (debug) {
        ggsave("/home/jakob/Desktop/debug/test.png", plot=plt)
        write.csv(comb_df, file="/home/jakob/Desktop/debug/coord.csv")
    }
    
    plt
}


generate_score_plots <- function(normal_run_entries, rt_run_entries, output_base, do_rt_run=TRUE) {
    
    score_funcs <- get_norm_methods("all")
    plot_ylabs <- get_y_label("all")

    print("Visualizing results...")
    print(length(score_funcs))
    
    score_plots <- list()
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
                                     verbose=FALSE)
        score_plots[[i]] <- plt
        
        norm_sub_path <- paste(output_base, tolower(y_lab), "norm", "csv", sep=".")
        rt_sub_path <- paste(output_base, tolower(y_lab), "rt", "csv", sep=".")
        write_entries_to_file(norm_sub_path, normal_run_entries)
        
        if (do_rt_run) {
            write_entries_to_file(rt_sub_path, rt_run_entries)
        }
    }
    
    score_plots
}


generate_normal_run_results <- function(nr, rs, row_name_cols) {

    used_methods_names <- getUsedMethodNames(nr)
    normalization_matrices <- getNormalizationMatrices(nr)
    
    run_entries <- list()
    for (i in 1:length(used_methods_names)) {
        
        print(paste("Method: ", used_methods_names[[i]]))
        
        run_result <- get_run_entry(nr, normalization_matrices[[i]], used_methods_names[[i]], 
                                    POT_PAT, rs, row_name_cols=row_name_cols)
        
        # browser()
        
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

    prepared_df <- get_prepared_normalyzer_sheet(nr, method_data, row_name_cols, replicate_nbrs, var_filter_frac=var_filter_frac)
    filter_df <- prepared_df[which(!grepl(NAME_FILTER, rownames(prepared_df))),]

    if (stat_test == "welch" || stat_test == "student") {
        target_sig_val_p <- get_anova_pvals(filter_df, custom_replicate_groups, ttest_type=stat_test, adjust_fdr=F)
        target_sig_val_fdr <- get_anova_pvals(filter_df, custom_replicate_groups, ttest_type=stat_test, adjust_fdr=T)
        print(paste("Running regular test: ", stat_test))
    }
    else if (stat_test == "limma") {
        target_sig_val_p <- get_limma_pvals(filter_df, custom_replicate_groups, adjust_fdr=F)
        target_sig_val_fdr <- get_limma_pvals(filter_df, custom_replicate_groups, adjust_fdr=T)
        print("Running Limma test")
    }
    else {
        stop(paste("Unknown stat_test:", stat_test))
    }
    
    diff_sig_p <- target_sig_val_p[which(grepl(target_pattern, names(target_sig_val_p)))]
    back_sig_p <- target_sig_val_p[which(!grepl(target_pattern, names(target_sig_val_p)))]
    
    diff_sig_fdr <- target_sig_val_fdr[which(grepl(target_pattern, names(target_sig_val_fdr)))]
    back_sig_fdr <- target_sig_val_fdr[which(!grepl(target_pattern, names(target_sig_val_fdr)))]
    
    if (omit_na) {
        nbr_potato_tot <- length(na.omit(diff_sig_fdr))
        nbr_potato_sig <- length(na.omit(diff_sig_fdr[which(diff_sig_fdr <= sig_thres)]))
        nbr_back_tot <- length(na.omit(back_sig_fdr))
        nbr_back_sig <- length(na.omit(back_sig_fdr[which(back_sig_fdr <= sig_thres)]))
    }
    else {
        nbr_potato_tot <- length(diff_sig_fdr)
        nbr_potato_sig <- length(diff_sig_fdr[which(diff_sig_fdr <= sig_thres)])
        nbr_back_tot <- length(back_sig_fdr)
        nbr_back_sig <- length(back_sig_fdr[which(back_sig_fdr <= sig_thres)])
    }

    total_rows <- length(na.omit(back_sig_fdr)) + length(na.omit(back_sig_fdr))
    
    entry <- EntryRow(nr=nr,
                      method_data=filter_df,
                      norm_method=name, 
                      tot_rows=total_rows,
                      target_tot=nbr_potato_tot,
                      target_sign=nbr_potato_sig,
                      background_tot=nbr_back_tot,
                      background_sign=nbr_back_sig,
                      rt_settings=current_rt_setting,    
                      sig_df=target_sig_val_fdr,
                      diff_sig_p=na.omit(diff_sig_p),
                      back_sig_p=na.omit(back_sig_p),
                      diff_sig_fdr=na.omit(diff_sig_fdr),
                      back_sig_fdr=na.omit(back_sig_fdr)
                     )
    
    entry
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


















