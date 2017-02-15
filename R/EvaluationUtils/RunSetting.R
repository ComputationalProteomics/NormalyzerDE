# source("EvaluationUtils/evaluationMain.R")


RunSetting <- function(sig_thres, do_fdr, rt_windows, window_shifts,
                       lowest_window_size, window_merge_method, sample_comp, stat_test, quiet=FALSE, subset=FALSE,
                       super_dirname=NULL) {
    
    me <- list(
        sig_thres = sig_thres,
        do_fdr = do_fdr,
        rt_windows = rt_windows,
        window_shifts = window_shifts,
        lowest_window_size = lowest_window_size,
        window_merge_method = window_merge_method,
        quiet = quiet,
        subset = subset,
        super_dirname = super_dirname,
        sample_comp = sample_comp,
        stat_test = stat_test
    )
    
    class(me) <- append(class(me), "RunSetting")
    return(me)
}

get_run_setting_descr_string <- function(rs) {
    
    if (rs$do_fdr) {
        fdr_str <- "fdr"
    }
    else {
        fdr_str <- "nofdr"
    }
    
    paste(rs$sig_thres, fdr_str, rs$window_shifts, rs$lowest_window_size, rs$window_merge_method, rs$stat_test, paste(rs$sample_comp, sep="", collapse="vs"), sep="_")
}

get_run_setting_base <- function(rs) {
    descr_string <- get_run_setting_descr_string(rs)
    out_dir <- generate_output_dir(descr_string, rs$super_dirname, create_dir=FALSE)
    run_base <- paste(out_dir, "/", descr_string, sep="")
    run_base
}

do_run_from_runsetting <- function(rs) {
    
    descr_string <- get_run_setting_descr_string(rs)
    
    # descr_string <- paste(rs$super_dirname, "/", orig_descr_string, sep="")

    out_dir <- generate_output_dir(descr_string, rs$super_dirname, create_dir=TRUE)
    # run_base <- paste(out_dir, "/", descr_string, sep="")

    print("Will do run!")
    run_review_data_test(rs$super_dirname, 
                         subset=rs$subset,
                         measure_type="all",
                         sig_thres=rs$sig_thres,
                         do_fdr=rs$do_fdr,
                         frame_shifts=rs$window_shifts,
                         lowest_window_size=rs$lowest_window_size,
                         window_merge_method=rs$window_merge_method,
                         rt_windows=rs$rt_windows,
                         target_replicates=rs$sample_comp,
                         stat_test=rs$stat_test)
}


generate_output_dir <- function(description, out_path, create_dir=TRUE) {
    datestamp <- get_datestamp_string()
    run_dir <- paste(out_path, "/", datestamp, "_", description, sep="")
    if (create_dir) {
        createDirectory(run_dir)
    }
    run_dir
}



