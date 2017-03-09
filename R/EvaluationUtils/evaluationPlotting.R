library("ggplot2")

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

