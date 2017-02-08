source("EvaluationUtils/evaluationUtils.R")





RunSetting <- function(sig_thres, do_fdr, rt_windows, window_shifts,
                       lowest_window_size, window_merge_method, quiet=FALSE, subset=FALSE) {
    
    me <- list(
        sig_thres = sig_thres,
        do_fdr = do_fdr,
        rt_windows = rt_windows,
        window_shifts = window_shifts,
        lowest_window_size = lowest_window_size,
        window_merge_method = window_merge_method,
        quiet = quiet,
        subset = subset
    )
    
    class(me) <- append(class(me), "RunSetting")
    return(me)
}

get_run_setting_descr_string <- function(rs) {
    paste(rs$sig_thres, rs$do_fdr, rs$window_shifts, rs$lowest_window_size, rs$window_merge_method, sep="_")
}

do_run_from_runsetting <- function(rs) {
    
    descr_string <- get_run_setting_descr_string(rs)
    out_dir <- generate_output_dir(descr_string, create_dir=TRUE)
    run_base <- paste(out_dir, "/", descr_string, sep="")
    
    # run_review_data_test()
    
    print("WIll do run!")
    
    run_review_data_test(run_base, 
                         subset=rs$subset,
                         measure_type="all",
                         sig_thres=rs$sig_thres,
                         do_fdr=rs$do_fdr,
                         frame_shifts=rs$window_shifts,
                         lowest_window_size=rs$lowest_window_size,
                         window_merge_method=rs$window_merge_method)
}
