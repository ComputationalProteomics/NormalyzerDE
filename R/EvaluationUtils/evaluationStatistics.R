
get_anova_pvals <- function(df, replicate_groups, ttest_type="welch", adjust_fdr=T) {
    
    replicate_groups <- as.factor(replicate_groups)
    
    stat_test_func <- function(sampleIndex) {
        
        rep_counts <- table(names(na.omit(sampleIndex)))
        
        if (length(rep_counts) == 2 && min(rep_counts > 1)) {
            if (ttest_type == "student") {
                t.test(unlist(sampleIndex)~replicate_groups, var.equal=TRUE)$p.value
            }
            else if (ttest_type == "welch") {
                t.test(unlist(sampleIndex)~replicate_groups)$p.value
            }
            else {
                stop(paste("Unknown test type:", ttest_type))
            }
        }
        else {
            NA
        }
    }
    
    anovaPVals <- apply(df, 1, stat_test_func)
    
    if (!adjust_fdr) {
        anovaPVals
    }
    else {
        stats::p.adjust(anovaPVals, method="BH")
    }
}


get_limma_pvals <- function(df, replicate_groups, adjust_fdr=T) {
    
    design_names <- c("c1", "c2")
    contrasts <- c("c1-c2")
    
    limma_fit <- get_limma_fit(df, factor(replicate_groups), design_names, contrasts)
    limma_table <- get_limma_table(limma_fit, comp_nbr=1)
    
    if (!adjust_fdr) {
        sig_vals <- limma_table[,"P.Value"]
    }
    else {
        sig_vals <- limma_table[,"adj.P.Val"]
    }
    
    names(sig_vals) <- rownames(df)
    sig_vals
}


get_limma_fit <- function(data_df, header, design_names, my_contrasts) {
    
    design <- model.matrix(~0+header)
    
    colnames(design) <- design_names
    matrix <- data.matrix(data_df)
    
    contrast.matrix <- makeContrasts(contrasts=my_contrasts, levels=design)
    
    fit <- lmFit(matrix, design)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    
    fit2
}

get_limma_table <- function(limma_fit, comp_nbr) {
    topTable(limma_fit, coef=comp_nbr, number=Inf, sort.by="none")
}