write_entries_to_file <- function(out_path, run_entries) {
    
    first_entry_v <- get_entry_vector(run_entries[[1]])
    out_strs <- matrix(ncol=length(first_entry_v), nrow=0)
    first_entry_head <- get_entry_header(run_entries[[1]])
    
    print(out_strs)
    print(get_entry_header(run_entries[[1]]))
    
    colnames(out_strs) <- get_entry_header(run_entries[[1]])
    
    for (i in 1:length(run_entries)) {
        
        entry <- run_entries[[i]]
        out_strs <- rbind(out_strs, get_entry_vector(entry))
    }
    
    write.csv(out_strs, file=out_path, quote=FALSE, row.names=FALSE)
}

write_significance_matrix <- function(sign_out_path, entries) {
    
    nbr_features <- length(entries[[1]]$sig_df)
    nbr_methods <- length(entries)
    
    sign_matrix <- data.frame(matrix(nrow=nbr_features, ncol=0))
    entry_names <- NULL
    
    browser()
    
    for (i in 1:nbr_methods) {
        
        entry <- entries[[i]]
        if (is.null(entry_names)) {
            entry_names <- names(entry$sig_df)
        }
        
        entry_col <- entry$sig_df
        if (!is.null(entry$rt_settings)) {
            name <- paste(entry$norm_method, entry$rt_settings, sep="_")
        }
        else {
            name <- entry$norm_method
        }
        
        sign_matrix[name] <- entry_col
    }
    
    feature_names <- entry_names
    unique_rownames <- paste(feature_names, seq(1, nbr_features), sep="_")
    
    n_t <- table(unique_rownames)
    print(n_t[n_t > 1])
    print(length(n_t))
    
    rownames(sign_matrix) <- unique_rownames
    
    write.table(sign_matrix, file=sign_out_path, sep="\t", quote=FALSE, col.names=NA)
}

write_data_matrices <- function(output_base, normal_run_entries, rt_run_entries, nr) {
    
    for (run_entry in normal_run_entries) {
        
        out_df <- cbind(sig_val=run_entry$sig_df, run_entry$method_data)
        output_path <- paste0(output_base, "_", run_entry$norm_method, ".data.tsv")
        write.table(out_df, file=output_path, sep="\t", quote=FALSE, col.names=NA)
    }
    
    for (run_entry in rt_run_entries) {
        
        out_df <- cbind(sig_val=run_entry$sig_df, run_entry$method_data)
        output_path <- paste0(output_base, "_", run_entry$norm_method, run_entry$rt_settings, ".data.tsv")
        write.table(out_df, file=output_path, sep="\t", quote=FALSE, col.names=NA)
    }
}

