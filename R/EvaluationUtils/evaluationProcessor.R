

source("EvaluationUtils/evaluationUtils.R")
source("EvaluationUtils/RunSetting.R")
source("utils.R")


do_small_run <- function(out_path="EvaluationUtils/automated_tests") {
    
    run_dir <- generate_output_dir("small")
    run_base <- paste(run_dir, "/", "small_out", sep="")
    
    print(paste("Writing content to:", run_dir))
    run_review_data_test(run_base, 
                         subset=TRUE,
                         measure_type="all",
                         sig_thres=0.1,
                         do_fdr=FALSE)
    print(paste("Content written to:", run_dir))
}

do_large_run <- function() {

    run_dir <- generate_output_dir("large")
    run_base <- paste(run_dir, "/", "large_out", sep="")
    
    print(paste("Writing content to:", run_dir))
    run_review_data_test(run_base, 
                         subset=FALSE,
                         measure_type="all",
                         sig_thres=0.1,
                         do_fdr=FALSE)
    print(paste("Content written to:", run_dir))
}


generate_output_dir <- function(description, out_path, create_dir=TRUE) {
    datestamp <- get_datestamp_string()
    run_dir <- paste(out_path, "/", datestamp, "_", description, sep="")
    if (create_dir) {
        createDirectory(run_dir)
    }
    run_dir
}

do_multiple_runs <- function(run_setting_list, max_cores, testrun_only=FALSE) {
    
    library("parallel")
    
    print(paste("Will process", length(run_setting_list), "entries on", max_cores, "cores"))
    
    # print(length(run_setting_list))
    
    if (testrun_only) {
        do_run_from_runsetting(run_setting_list[[1]])
        stop("")
    }
    
    total_runs <- length(run_setting_list)
    current_run <- 1
    while (current_run <= total_runs) {

        running_procs <- list()
        
        for (core in 1:max_cores) {

            current_runsetting <- run_setting_list[[current_run]]
            print(paste("Starting processing run:", current_run, "/", total_runs, "at", Sys.time()))
            # running_procs[[core]] <- mcparallel(do_run_from_runsetting(current_runsetting))
            running_procs[[core]] <- mcparallel(do_run_from_runsetting(current_runsetting), detached=current_runsetting$quiet)
            
            current_run <- current_run + 1
            if (current_run > total_runs) {
                break
            }
            
            # Sys.sleep(2)
        }

        print("Waiting for running processes to finish...")
        mccollect(running_procs)
    }
}


screen_values <- function(sig_thresholds, do_fdr, rt_windows, window_shifts, 
                          lowest_window_size, window_merge_method, small_run=FALSE, max_cores=1,
                          quiet_processing=TRUE, subset=FALSE, super_dirname=NULL) {
    
    index <- 1
    total_runs <- length(sig_thresholds) * length(do_fdr) * length(window_shifts) * length(lowest_window_size) * length(window_merge_method)

    # ongoing_runs <- list()
    # for (core in 1:max_cores) {
        
    run_settings <- list()
        
    for (sig in sig_thresholds) {
        for (fdr in do_fdr) {
            for (win_shift in window_shifts) {
                for (low_win_size in lowest_window_size) {
                    for (merge_method in window_merge_method) {
                        
                        # print(paste("Screening sig", sig, "fdr", fdr, "win_shift", win_shift, "low_win_size", low_win_size, "merge_method", merge_method))

                        run_settings[[index]] <- RunSetting(sig_thres=sig, 
                                                            do_fdr=fdr, 
                                                            rt_windows=rt_windows, 
                                                            window_shifts=win_shift, 
                                                            lowest_window_size=low_win_size, 
                                                            window_merge_method=merge_method,
                                                            quiet=quiet_processing,
                                                            subset=subset,
                                                            super_dirname=super_dirname)
                        
                        index <- index + 1       
                    }
                }
            }
        }
    }
    
    do_multiple_runs(run_settings, max_cores)
    # print(paste("Total screened: ", index - 1))
}

hardcoded_screen_values <- function(do_full_run, super_dirname) {
    
    if (do_full_run) {
        sig_thres <- c(0.01, 0.05, 0.1, 0.2)
        do_fdrs <- c(TRUE, FALSE)
        rt_windows <- c(seq(0.2, 1, 0.2), seq(2, 5, 0.5), seq(5, 20, 1), seq(22, 30, 2))
        nbr_frame_shifts <- c(1, 3, 5)
        fix_window_bottom <- c(1, 2, 5, 10, 20)
        merge_method <- c("mean", "median")
        max_cores <- 7
        quiet <- FALSE
        subset <- FALSE
    }
    else {
        sig_thres <- c(0.1)
        do_fdrs <- c(FALSE)
        rt_windows <- c(1,2,3)
        nbr_frame_shifts <- c(1, 3)
        fix_window_bottom <- c(1, 5)
        merge_method <- c("mean", "median")
        max_cores <- 4
        quiet <- FALSE
        subset <- TRUE
    }
    
    screen_values(sig_thresholds = sig_thres,
                  do_fdr = do_fdrs,
                  rt_windows = rt_windows,
                  window_shifts = nbr_frame_shifts,
                  lowest_window_size = fix_window_bottom,
                  window_merge_method = merge_method,
                  max_cores = max_cores,
                  quiet_processing = quiet,
                  subset = subset,
                  super_dirname = super_dirname)
}
















