source("EvaluationUtils/evaluationMain.R")
source("EvaluationUtils/RunSetting.R")
source("EvaluationUtils/evaluationUtils.R")
source("EvaluationUtils/hardcoded_screen_values.R")
source("utils.R")


do_small_run <- function(out_path="EvaluationUtils/regular_tests") {
    
    run_dir <- generate_output_dir("small", out_path)
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


do_multiple_runs <- function(run_setting_list, max_cores, testrun_only=FALSE) {
    
    library("parallel")
    
    print(paste("Will process", length(run_setting_list), "entries on", max_cores, "cores"))
    
    if (testrun_only) {
        for (i in length(run_setting_list)) {
            do_run_from_runsetting(run_setting_list[[i]])
        }
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
        }

        print("Waiting for running processes to finish...")
        mccollect(running_procs)
    }
}


screen_values <- function(sig_thresholds, do_fdr, rt_windows, window_shifts, 
                          lowest_window_size, window_merge_method, stat_test=c("anova"), small_run=FALSE, max_cores=1,
                          quiet_processing=TRUE, subset=FALSE, super_dirname=NULL, sample_comparisons=list(c(2,3)),
                          debug=F) {
    
    index <- 1
    total_runs <- length(sig_thresholds) * length(do_fdr) * length(window_shifts) * length(lowest_window_size) * length(window_merge_method)

    run_settings <- list()
        
    for (sig in sig_thresholds) {
        for (fdr in do_fdr) {
            for (win_shift in window_shifts) {
                for (low_win_size in lowest_window_size) {
                    for (merge_method in window_merge_method) {
                        for (sample_comp in sample_comparisons) {
                            for (test in stat_test) {
                                # print(paste("Screening sig", sig, "fdr", fdr, "win_shift", win_shift, "low_win_size", low_win_size, "merge_method", merge_method))
                                
                                run_settings[[index]] <- RunSetting(sig_thres=sig, 
                                                                    do_fdr=fdr, 
                                                                    rt_windows=rt_windows, 
                                                                    window_shifts=win_shift, 
                                                                    lowest_window_size=low_win_size, 
                                                                    window_merge_method=merge_method,
                                                                    quiet=quiet_processing,
                                                                    subset=subset,
                                                                    super_dirname=super_dirname,
                                                                    sample_comp=sample_comp,
                                                                    stat_test=test)

                                index <- index + 1
                            }
                        }
                    }
                }
            }
        }
    }
    
    do_multiple_runs(run_settings, max_cores, testrun_only = debug)
    # print(paste("Total screened: ", index - 1))
}


















