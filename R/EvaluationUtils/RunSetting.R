# source("EvaluationUtils/evaluationMain.R")


RunSetting <- function(sig_thres, do_fdr, rt_windows, window_shifts, var_filter_frac,
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
        stat_test = stat_test,
        var_filter_frac = var_filter_frac,
        run_dir = NULL
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
    
    paste(rs$sig_thres, fdr_str, rs$window_shifts, rs$lowest_window_size, rs$window_merge_method, rs$stat_test, 
          rs$var_filter_frac, paste(rs$sample_comp, sep="", collapse="vs"), sep="_")
}

get_run_setting_base <- function(rs) {
    
    descr_string <- get_run_setting_descr_string(rs)
    
    if (is.null(rs$run_dir)) {
        stop("No run_dir assigned!")
    }
    
    run_base <- paste(rs$run_dir, "/", descr_string, sep="")
    run_base
}

do_run_from_runsetting <- function(rs) {
    
    descr_string <- get_run_setting_descr_string(rs)
    
    out_dir <- generate_output_dir(rs)

    if (length(rs$rt_windows) > 0) {
        do_rt_run <- TRUE
    }
    else {
        do_rt_run <- FALSE
    }
    
    print("Will do run!")
    run_review_data_test(rs$super_dirname,
                         rs$run_dir,
                         subset=rs$subset,
                         measure_type="all",
                         sig_thres=rs$sig_thres,
                         do_fdr=rs$do_fdr,
                         frame_shifts=rs$window_shifts,
                         lowest_window_size=rs$lowest_window_size,
                         window_merge_method=rs$window_merge_method,
                         rt_windows=rs$rt_windows,
                         target_replicates=rs$sample_comp,
                         stat_test=rs$stat_test,
                         var_filter_frac=rs$var_filter_frac,
                         do_rt_run=do_rt_run)
}


setup_output_dir_path <- function(rs) {
    
    out_path <- rs$super_dirname
    
    description <- get_run_setting_descr_string(rs)
    datestamp <- get_datestamp_string()
    run_dir <- paste(out_path, "/", datestamp, "_", description, sep="")
    rs$run_dir <- run_dir
    rs
}


generate_output_dir <- function(rs) {
    
    createDirectory(rs$run_dir)
}



