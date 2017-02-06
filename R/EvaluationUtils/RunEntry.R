header <- c("method", "tot_na_reduced", "rt_settings", "potato_tot", "potato_sig", "back_tot", "back_sig")

EntryRow <- function(norm_method, tot_rows, target_tot, target_sign, 
                     background_tot, background_sign, rt_settings=NULL) {
    
    me <- list(
        norm_method = norm_method,
        tot_rows = tot_rows,
        target_tot = target_tot,
        target_sign = target_sign,
        background_tot = background_tot,
        background_sign = background_sign,
        rt_settings = rt_settings
    )
    
    class(me) <- append(class(me), "EntryRow")
    return(me)
}

get_entry_vector <- function(e) {
    
    base <- c(e$norm_method,
              e$tot_rows,
              e$target_tot,
              e$target_sign,
              e$background_tot,
              e$background_sign)
    
    if (!is.null(e$rt_settings)) {
        base <- c(base, e$rt_settings)
    }
    
    base
}

get_entry_header <- function(e) {
    
    header <- c("method", "tot_na_reduced", "rt_settings", "potato_tot", "potato_sig", "back_tot", "back_sig")
    if (!is.null(e$rt_settings)) {
        header <- c(header, e$rt_settings)
    }
    header
}

get_precision <- function(e) {
    
    target_tot <- e$target_tot
    target_sign <- e$target_sign
    background_tot <- e$background_tot
    background_sign <- e$background_sign
    
    e$target_sign / (e$target_sign + e$background_sign)
}

get_recall <- function(e) {
    
    e$target_sign / e$target_tot
}

get_f_score <- function(e) {

    prec <- get_precision(e)
    recall <- get_recall(e)
    
    2 * (prec * recall) / (prec + recall)
}

get_entry_string <- function(e) {
    target_fields <- c(e$norm_method, e$tot_rows, 
                       e$target_tot, e$target_sign,
                       e$background_tot, e$background_sign,
                       e$rt_settings)
    
    paste(target_fields, sep="\t")
}


