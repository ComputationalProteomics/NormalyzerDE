

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




header <- c("method", "tot_na_reduced", "rt_settings", "potato_tot", "potato_sig", "back_tot", "back_sig")

EntryRow <- function(norm_method, tot_rows, target_tot, target_sign, 
                     background_tot, background_sign, rt_settings=NULL) {
    
    me <- list(
        norm_method = norm_method,
        tot_rows = tot_rows,
        target_tot = target_tot,
        target_sign = target_sign,
        background_tot = background_tot,
        background_sign = background_sign,
        rt_settings = rt_settings
    )
    
    class(me) <- append(class(me), "EntryRow")
    return(me)
}

get_entry_vector <- function(e) {
    
    base <- c(e$norm_method,
      e$tot_rows,
      e$target_tot,
      e$target_sign,
      e$background_tot,
      e$background_sign)
    
    if (!is.null(e$rt_settings)) {
        base <- c(base, e$rt_settings)
    }
    
    base
}

get_entry_header <- function(e) {
    
    header <- c("method", "tot_na_reduced", "rt_settings", "potato_tot", "potato_sig", "back_tot", "back_sig")
    if (!is.null(e$rt_settings)) {
        header <- c(header, e$rt_settings)
    }
    header
}

run_review_data_test <- function(subset=FALSE, output_path=NULL, verbose=FALSE, plot_path=NULL, do_normal_run=TRUE, do_rt_run=TRUE) {
    
    start_time <- Sys.time()
    
    if (subset) review_data_path <- paste("../tests/data", subset_report, sep="/")
    else review_data_path <- paste("../tests/data", full_report, sep="/")
    
    print(paste("Processing dataset:", review_data_path))
    
    job_name <- "review_evaluation"
    norm_obj <- getVerifiedNormalyzerObject(review_data_path, job_name)

    default_rt_window <- 2
    
    if (subset) {
        rt_windows <- c(1, 2, 3)
    }
    else {
        rt_windows <- c(0.2, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 7, 8, 9, 10, 13, 16, 20)
    }
    
    nr <- normMethods(norm_obj, job_name, normalizeRetentionTime=FALSE, runNormfinder=FALSE)
    
    normal_run_entries <- generate_normal_run_results(nr, adjust_fdr=FALSE)
    rt_run_entries <- generate_rt_run_results(nr, rt_windows, adjust_fdr=FALSE)
    
    print(paste("Normal entries:", length(normal_run_entries)))
    print(paste("RT entries:", length(rt_run_entries)))
    
    # stop("Testing - for now")
    
    new_visualize_results(plot_path, normal_run_entries, rt_run_entries)

    stop("Testing - for now")
        
    result_dfs <- slice_results_to_dfs(result_df)

    if (is.null(plot_path)) {
        print("No plot path specified - No plotting")
    }
    else {
        visualize_results(result_dfs, plot_path, "rt_settings", "frac_p", 
                          title="RT normalizations with varying time window, average of overlapping windows",
                          xlab="RT window size (minutes)",
                          ylab="Amount of features detected as differentially expressed (%)")
    }
    
    end_time <- Sys.time()
    print(difftime(end_time, start_time))
    print(paste("Start:", start_time, "End:", end_time))
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

generate_rt_run_results <- function(nr, rt_windows, sig_thres=0.1, adjust_fdr=TRUE) {
    
    used_methods_names <- getUsedMethodNames(nr)
    normalization_matrices <- getNormalizationMatrices(nr)
    
    normalyzer_filterraw <- nr@nds@filterrawdata
    retention_times <- nr@nds@retentionTimes

    run_entries <- list()
    
    print("Performing RT runs")

    # run_entries_base_index <- 1
    for (rt_window in rt_windows) {

        print(paste("Method: RTs, window size:", rt_window))

        median_rt <- getSmoothedRTNormalizedMatrix(normalyzer_filterraw, retention_times, medianNormalization, rt_window)
        median_entry <- get_run_entry(nr, median_rt, "RT-median", POT_PAT, HUMAN_PAT, sig_thres=sig_thres, adjust_fdr=adjust_fdr, rt_settings=rt_window)

        mean_rt <- getSmoothedRTNormalizedMatrix(normalyzer_filterraw, retention_times, meanNormalization, rt_window)
        mean_entry <- get_run_entry(nr, mean_rt, "RT-mean", POT_PAT, HUMAN_PAT, sig_thres=sig_thres, adjust_fdr=adjust_fdr, rt_settings=rt_window)

        loess_rt <- getSmoothedRTNormalizedMatrix(normalyzer_filterraw, retention_times, performCyclicLoessNormalization, rt_window)
        loess_entry <- get_run_entry(nr, loess_rt, "RT-loess", POT_PAT, HUMAN_PAT, sig_thres=sig_thres, adjust_fdr=adjust_fdr, rt_settings=rt_window)

        run_entries_base_index <- length(run_entries) + 1
          
        run_entries[[run_entries_base_index]] <- median_entry
        run_entries[[run_entries_base_index + 1]] <- mean_entry
        run_entries[[run_entries_base_index + 2]] <- loess_entry
        
        # run_entries <- c(run_entries, median_entry, mean_entry, loess_entry)
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


# # Returns vector on the form: short_header <- c("method", "tot_na_reduced", "rt_settings", "pot_or_back",
# # "tot", "nbr_sig_p", "frac_p", "nbr_sig_fdr", "frac_fdr")
# get_significance_vector <- function(analysis_name, tot_na_reduced, anova_vals, sub_pattern, sig_thres, pot_or_back, inverse_match=FALSE, rt_settings=NULL) {
# 
#     if (!inverse_match) {
#         sub_sig <- anova_vals[which(grepl(sub_pattern, names(anova_vals)))]
#     }
#     else {
#         sub_sig <- anova_vals[which(!grepl(sub_pattern, names(anova_vals)))]
#     }
#     
#     tot_sub <- length(sub_sig)
#     sub_n_sig <- length(sub_sig[which(sub_sig <= sig_thres)])
#     sub_sig_frac <- round(100 * sub_n_sig / tot_sub, 1)
# 
#     if (!is.null(rt_settings)) {
#         result_vector <- c(analysis_name, tot_na_reduced, rt_settings, tot_sub, sub_n_sig, sub_sig_frac)
#     }
#     else {
#         result_vector <- c(analysis_name, tot_na_reduced, tot_sub, sub_n_sig, sub_sig_frac)
#     }
# }


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


new_visualize_results <- function(output_path, normal_objects, rt_objects, title="[title]", xlab="[xlab]", ylab="[ylab]") {
    
    # Visualize normal entry with ggplot
    
    dark_colors <- c("darkgreen", "darkblue", "darkred")
    bright_colors <- c("green", "blue", "red")
    
    all_colors <- c("green", "darkgreen", "blue", "darkblue", "red", "darkred")
    techniques <- c("Loess, background", "Loess, potato", "Mean, background", "Mean, potato", "Median, background", "Median, potato")
    
    plt <- ggplot() #+ scale_colour_manual(labels=techniques, values=all_colors)
    plt <- plt + ggtitle(title)
    plt <- plt + xlab(xlab)
    plt <- plt + ylab(ylab)
    
    full_df <- data.frame(matrix(nrow=0, ncol=3))
    colnames(full_df) <- c("x", "y_level", "norm_method")
    
    # RT plotting
    
    for (i in 1:length(rt_objects)) {
        
        rt_obj <- rt_objects[[i]]
        
        print(paste("Target obj:", rt_obj))
        
        y_level <- rt_obj$target_sign / rt_obj$target_tot
        
        x <- rt_obj$rt_settings
        print(x)
        
        # norm_method <- paste(rt_obj$norm_method, rt_obj$rt_settings, sep="_")
        
        print(paste("Norm method: ", rt_obj$norm_method))
        
        df <- data.frame(x, y_level, rt_obj$norm_method)
        colnames(df) <- c("x", "y_level", "norm_method")
        
        full_df <- rbind(full_df, df)
        
    }
    
    # Normal plotting
    
    min_rt_setting <- min(sapply(rt_run_entries, function(x) {x$rt_settings}))
    max_rt_setting <- max(sapply(rt_run_entries, function(x) {x$rt_settings}))
    x <- c(min_rt_setting, max_rt_setting)
    
    for (i in 1:length(normal_objects)) {
        
        obj <- normal_objects[[i]]
        y_level <- obj$target_sign / obj$target_tot
        df <- data.frame(x, y_level, obj$norm_method)
        colnames(df) <- c("x", "y_level", "norm_method")
        full_df <- rbind(full_df, df)
    }

    plt <- plt + geom_line(data=full_df, aes_string(x="x", y="y_level", colour="norm_method"))
    plt <- plt + geom_point(data=full_df, aes_string(x="x", y="y_level", colour="norm_method"))
    
    plt
}

# visualize_results <- function(result_dfs, output_path, x_col, y_col, 
#                               title="Default title", xlab="Default xlab", ylab="Default ylab") {
# 
#     dark_colors <- c("darkgreen", "darkblue", "darkred")
#     bright_colors <- c("green", "blue", "red")
#     
#     all_colors <- c("green", "darkgreen", "blue", "darkblue", "red", "darkred")
#     techniques <- c("Loess, background", "Loess, potato", "Mean, background", "Mean, potato", "Median, background", "Median, potato")
#     
#     plt <- ggplot() + scale_colour_manual(labels=techniques, values=all_colors)
#     plt <- plt + ggtitle(title)
#     plt <- plt + xlab(xlab)
#     plt <- plt + ylab(ylab)
#     
#     for (result_df in result_dfs) {
#         
#         plt <- plt + geom_line(data=result_df, aes_string(x=x_col, y=y_col, colour="method_and_settings"))
#         plt <- plt + geom_point(data=result_df, aes_string(x=x_col, y=y_col, colour="method_and_settings"))
#     }
#     
#     ggsave(output_path)
# }
# 
# slice_results_to_dfs <- function(result_m) {
#     
#     result_dfs <- list()
# 
#     convert_to_numericals <- c("tot_na_reduced", "rt_settings", "pot_or_back",
#                                "tot", "nbr_sig_p", "frac_p", "nbr_sig_fdr", "frac_fdr")
#     
#     method_settings_col <- paste0(result_m$method, result_m$pot_or_back)
#     result_m <- cbind(result_m, method_and_settings=method_settings_col)
# 
#     unique_techniques <- unique(result_m$method_and_settings)
#     
#     list_index <- 1
#     for (tech in unique_techniques) {
#         
#         tmp_df <- result_m[which(result_m$method_and_settings == tech),]
#         tmp_df[["rt_settings"]] <- as.numeric(tmp_df[["rt_settings"]])
#         
#         for (field in convert_to_numericals) {
#             tmp_df[[field]] <- as.numeric(tmp_df[[field]])
#         }
#         
#         result_dfs[[list_index]] <- tmp_df
#         list_index <- list_index + 1
#     }
#     
#     result_dfs
# }






