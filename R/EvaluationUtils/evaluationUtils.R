
get_y_label <- function(measure_type) {
    
    if (measure_type == "fscore") {
        plot_ylabs <- "F-score"
    }
    else if (measure_type == "precision") {
        plot_ylabs <- "Precision"
    }
    else if (measure_type == "recall") {
        plot_ylabs <- "Recall"
    }
    else if (measure_type == "all") {
        plot_ylabs <- c("F-score", "Precision", "Recall")
    }
    else {
        stop(paste("Unknown measure_type:", measure_type))
    }
    
    plot_ylabs    
}

get_norm_methods <- function(measure_type) {
    
    if (measure_type == "fscore") {
        score_funcs <- get_f_score
    }
    else if (measure_type == "precision") {
        score_funcs <- get_precision
    }
    else if (measure_type == "recall") {
        score_funcs <- get_recall
    }
    else if (measure_type == "all") {
        score_funcs <- list(get_f_score, get_precision, get_recall)
    }
    else {
        stop(paste("Unknown measure_type:", measure_type))
    }
    
    score_funcs
}

get_rownames <- function(nr, row_name_cols) {
    
    raw_data <- nr@nds@rawData
    names_cols <- raw_data[-1:-2, row_name_cols]
    row_names <- paste(names_cols[,2], names_cols[,1], sep="|")
    row_names
}

get_colnames <- function(nr) {
    raw_data <- nr@nds@rawData
    header_row <- raw_data[2,which(raw_data[1,] > 0)]
    
    # print(paste("Found header row:", header_row))
    header_row
}









