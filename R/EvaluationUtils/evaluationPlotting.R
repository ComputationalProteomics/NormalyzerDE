library("ggplot2")
library("RColorBrewer")


generate_roc_plot <- function(normal_run_entries, rt_run_entries, use_fdr=T, only_normal=F, debug=T, sig_coord_thres=0.1, y_cutoff=0, title="ROC curve") {
    
    comb_df <- data.frame(pvals=c(), type=c(), method=c())
    sig_coord_df <- data.frame(x=c(), y=c(), method=c())
    
    if (!only_normal) {
        all_entries <- c(normal_run_entries, rt_run_entries)
    }
    else {
        all_entries <- normal_run_entries
    }
    
    plot_ylabs <- get_y_label("all")
    aucs <- c()

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
            
            if (is.null(sig_coord) && row["pvals"] > sig_coord_thres) {
                sig_coord <- data.frame(x=last_x_coord, y=last_y_coord, method=target_obj$detailed_name)
                sig_coord_df <- rbind(sig_coord_df, sig_coord)
            }
        }
        
        coord_df <- data.frame(FPR=x_vals, TPR=y_vals, sample=target_obj$detailed_name)
        comb_df <- rbind(comb_df, coord_df)
        
        aucs <- c(aucs, calculate_auc(x_vals, y_vals))
    }
    
    orig_lvs <- levels(comb_df$sample)
    new_lvs <- paste(round(aucs, 3), orig_lvs, sep=", ")
    levels(comb_df$sample) <- new_lvs
    sig_coord_df$method <- new_lvs
        
    colorCount <- length(unique(comb_df$sample))
    getPalette <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))
    
    plt <- ggplot() + 
        geom_line(data=comb_df, aes(FPR, TPR, color=sample)) + 
        ggtitle(title) +
        geom_point(data=sig_coord_df, aes(x=x, y=y, color=method)) +
        ylim(y_cutoff, NA) +
        scale_color_manual(values=getPalette(colorCount))
    
    if (debug) {
        ggsave("/home/jakob/Desktop/debug/test.png", plot=plt)
        write.csv(comb_df, file="/home/jakob/Desktop/debug/coord.csv")
    }
    
    plt
}

calculate_auc <- function(x_coord, y_coord) {

    tot_sum <- 0
    prev_x <- x_coord[1]
    for (i in 2:length(x_coord)) {
        
        x <- x_coord[i]
        y <- y_coord[i]
        x_diff <- x - prev_x
        
        if (x_diff != 0) {
            diff_area <- y * x_diff
            # print(paste("y", y, "x_diff", x_diff, "diff_area", diff_area))        
            tot_sum <- tot_sum + diff_area
        }
        prev_x <- x
    }
    
    # print(tot_sum)
    tot_sum
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


visualize_results_new <- function(normal_objects, rt_objects, target_measure_func,
                                  title="[title]", xlab="[xlab]", ylab="[ylab]",
                                  verbose=FALSE) {
    
    plt <- ggplot()  #+ scale_colour_manual(labels=techniques, values=all_colors)
    plt <- plt + ggtitle(title)
    plt <- plt + xlab(xlab)
    plt <- plt + ylab(ylab)

    full_df <- data.frame(matrix(nrow=0, ncol=3))
    colnames(full_df) <- c("x", "y_level", "norm_method")
    
    if (length(rt_objects) > 0) {
        min_rt_setting <- min(sapply(rt_objects, function(x) { x$rt_settings }))
        max_rt_setting <- max(sapply(rt_objects, function(x) { x$rt_settings }))
        x_lim <- c(min_rt_setting, max_rt_setting)
        rt_df <- get_rt_plotting_df(rt_objects, target_measure_func)
    }
    else {
        x_lim <- c(0,1)
    }
    
    normal_df <- get_normal_plotting_df(normal_objects, x_lim, target_measure_func)
    
    if (length(rt_objects) > 0) {
        full_df <- rbind(full_df, rt_df, normal_df)
    }
    else {
        full_df <- rbind(full_df, normal_df)
    }

    if (verbose) {
        print("--- Normal entries ---")
        print(normal_df)
        
        if (length(rt_objects) > 0) {
            print("--- RT entries ---")
            print(rt_df)
        }
    }
    
    plt <- plt + geom_line(data=full_df, aes_string(x="x", y="y_level", colour="norm_method")) + ylim(0, 1)
    plt <- plt + geom_point(data=full_df, aes_string(x="x", y="y_level", colour="norm_method"))

    plt
}

get_rt_plotting_df <- function(rt_objects, target_measure_func) {
    
    rt_df <- data.frame(matrix(nrow=0, ncol=3))
    for (i in 1:length(rt_objects)) {
        
        rt_obj <- rt_objects[[i]]
        y_level <- target_measure_func(rt_obj)
        x <- rt_obj$rt_settings
        df <- data.frame(x, y_level, rt_obj$norm_method)
        colnames(df) <- c("x", "y_level", "norm_method")
        rt_df <- rbind(rt_df, df)
    }
    
    rt_df
}

get_normal_plotting_df <- function(normal_objects, x_lim, target_measure_func) {
    
    normal_df <- data.frame(matrix(nrow=0, ncol=3))
    for (i in 1:length(normal_objects)) {
        
        obj <- normal_objects[[i]]
        y_level <- target_measure_func(obj)
        df <- data.frame(x_lim, y_level, obj$norm_method)
        colnames(df) <- c("x", "y_level", "norm_method")
        normal_df <- rbind(normal_df, df)
    }
    
    normal_df
}






