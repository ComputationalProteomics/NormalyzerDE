# require("EvaluationUtils/evaluationMain.R")
# require("EvaluationUtils/RunSetting.R")
# require("EvaluationUtils/evaluationUtils.R")
# require("EvaluationUtils/hardcoded_screen_values.R")
# require("utils.R")

hardcoded_screen_values <- function(do_full_run, super_dirname, subset=T, debug=F) {
    
    if (do_full_run) {
        sig_thres <- c(0.1, 0.05)
        do_fdr <- c(TRUE)
        # rt_windows <- c(1,2)
        rt_windows <- c(seq(1, 20, 2))
        frame_shifts <- c(1, 3, 5)
        lowest_window_size <- c(100, 200, 500)
        window_merge_method <- c("median")
        max_cores <- 7
        target_replicates <- list(c(2,3))
        stat_test <- c("welch")
        var_filter_fracs <- c(0)
        
        quiet <- FALSE
    }
    else {
        sig_thres <- c(0.1)
        do_fdr <- c(TRUE)
        # rt_windows <- c(seq(0.5, 5, 0.5), seq(6, 15, 1))
        rt_windows <- c(1, 3, 5, 7, 10)
        frame_shifts <- c(1)
        lowest_window_size <- c(100)
        window_merge_method <- c("median")
        max_cores <- 7
        target_replicates <- list(c(5,2), c(5,3), c(5,4), c(5,6), c(5,7), c(5,8), c(5,9))
        stat_test <- c("welch")
        quiet <- FALSE
        var_filter_fracs <- c(0.5)
        # var_filter_fracs <- c(seq(0,1,0.3))
        
        verbose <- TRUE
    }
    
    screen_values(sig_thresholds = sig_thres,
                  do_fdr = do_fdr,
                  rt_windows = rt_windows,
                  window_shifts = frame_shifts,
                  lowest_window_size = lowest_window_size,
                  window_merge_method = window_merge_method,
                  max_cores = max_cores,
                  quiet_processing = quiet,
                  subset = subset,
                  super_dirname = super_dirname,
                  sample_comparisons = target_replicates,
                  stat_test = stat_test,
                  var_filter_fracs = var_filter_fracs,
                  debug = debug)
}