

# Input: Normalized matrices

# We want to target particular 7 vs 7 columns
# Also, we want to slice out potato-named rows into a distinctly separate dataframe

# Let's test this first in an R terminal

source("NormalyzerDataset.R")
source("NormalizationEvaluationResults.R")
source("NormalyzerResults.R")

source("utils.R")
source("inputVerification.R")
source("analyzeResults.R")
source("normMethods.R")
source("normfinder-pipeline.R")

full_report <- "FeatureReport_BA_MQP_20140828-fornormalyzer_annotated.txt"
subset_report <- "MQP.subset_500.tsv"
# subset_report <- "BA_MQP.subset_500.tsv"

potato_pattern <- "^sol"
sample_pattern <- "dilA_[23]"
custom_replicate_groups <- c(2,2,2,2,2,2,2,3,3,3,3,3,3,3)

run_review_data_test <- function(subset=FALSE, output_path=NULL) {
    
    start_time <- Sys.time()
    
    if (subset) review_data_path <- paste("../tests/data", subset_report, sep="/")
    else review_data_path <- paste("../tests/data", full_report, sep="/")
    
    job_name <- "review_evaluation"

    print(paste("Processing dataset:", review_data_path))
    
    normalizeRetentionTime <- TRUE
    retentionTimeWindow <- 1
    
    print("Get normalyzer object")
    norm_obj <- getVerifiedNormalyzerObject(review_data_path, job_name)
    # job_dir <- setupJobDir(job_name)

    print("Perform analyzations")
    nr <- normMethods(norm_obj, jobName, normalizeRetentionTime=normalizeRetentionTime,
                      retentionTimeWindow=retentionTimeWindow)
    
    # first_test_matrix <- nr@data2log2
    
    used_methods_names <- getUsedMethodNames(nr)
    normalization_matrices <- getNormalizationMatrices(nr)

    header <- c("method", "pot_de", "pot_tot", "pot_frac", "non_pot_de", "non_pot_tot", "non_pot_frac")
    output <- matrix(ncol=7, nrow=0)
    colnames(output) <- header
    
    for (i in 1:length(used_methods_names)) {
        print(paste("Method: ", used_methods_names[[i]]))
        analysis_vector <- get_analysis_vector(nr, normalization_matrices[[i]], used_methods_names[[i]])
        output <- rbind(output, analysis_vector)
    }
    
    print(output)
    
    if (!is.null(output_path)) {
        write.csv(output, file=output_path, quote=FALSE)
    }
    
    # nr <- analyzeNormalizations(nr, jobName)
    end_time <- Sys.time()
    print(difftime(end_time, start_time))
    print(paste("Start:", start_time, "End:", end_time))
}


get_analysis_vector <- function(nr, method_data, name) {
    
    row_name_col <- 3
    potato_sub <- get_normalyzer_df_subset(nr, method_data, sample_pattern, potato_pattern, row_name_col, inverse_name_pattern=FALSE)
    non_potato_sub <- get_normalyzer_df_subset(nr, method_data, sample_pattern, potato_pattern, row_name_col, inverse_name_pattern=TRUE)
    
    potato_sub_na_reduced <- na_filter(potato_sub)
    non_potato_sub_na_reduced <- na_filter(non_potato_sub)
    
    potato_de <- get_number_de_anova(potato_sub_na_reduced, custom_replicate_groups)
    non_potato_de <- get_number_de_anova(non_potato_sub_na_reduced, custom_replicate_groups)
    
    potato_frac <- round(100 * length(potato_de) / nrow(potato_sub_na_reduced), 1)
    non_potato_frac <- round(100 * length(non_potato_de) / nrow(non_potato_sub_na_reduced), 1)
    
    analysis_vector <- c(name, length(potato_de), nrow(potato_sub_na_reduced), potato_frac, length(non_potato_de), nrow(non_potato_sub_na_reduced), non_potato_frac)
    analysis_vector
    # print(paste0("No potato: ", length(non_potato_de), "/", nrow(non_potato_sub_na_reduced), " - ", non_potato_frac))
}


get_normalyzer_df_subset <- function(nr, df, col_pattern, row_pattern, row_name_col, 
                                     inverse_name_pattern=FALSE) {
    
    raw_data <- nr@nds@rawData

    header_row <- 2
    if (inverse_name_pattern) {
        target_rows <- which(!grepl(row_pattern, raw_data[-1:-2, row_name_col]))
    }
    else {
        target_rows <- which(grepl(row_pattern, raw_data[-1:-2, row_name_col]))
    }

    target_cols <- which(grepl(col_pattern, raw_data[header_row,]))
    df[target_rows, target_cols]
}


na_filter <- function(df) {
    na_per_line <- rowSums(is.na(df))
    df_na_removed <- df[na_per_line < ncol(df) / 2, ]
    df_na_removed
}


get_number_de_anova <- function(df, replicate_groups, do_fdr=FALSE, thres=0.05) {
    
    # Filter NAs
    # na_per_line <- rowSums(is.na(df))
    # df_na_removed <- df[na_per_line < ncol(df) / 2, ]
    
    anova_func <- function(sampleIndex) {
        summary(stats::aov(unlist(sampleIndex)~replicate_groups))[[1]][[5]][1]
    }
    
    anovaPVal <- apply(df, 1, anova_func)
    
    anovaPVal[anovaPVal < 0.05]
    
    
    # anovaFDR <- cbind(anovaFDR, stats::p.adjust(anovaPVal[, methodIndex], method="BH"))

}


# norm_obj@rawData[2,] - Column names
# norm_obj@rawData[,3] - Here we target names

# rd[which(rd[,3] == "sol72"),]
# rd[which(grepl("^sol", rd[,3])),]

# > length(rd[which(grepl("^", rd[,3])),])
# [1] 507093
# > length(rd[which(grepl("^", rd[,3], invert=TRUE)),])
# Error in grepl("^", rd[, 3], invert = TRUE) : 
#     unused argument (invert = TRUE)
# > length(rd[which(!grepl("^", rd[,3])),])
# [1] 0
# > length(rd[which(!grepl("^sol", rd[,3])),])
# [1] 487512
# > length(rd[which(grepl("^sol", rd[,3])),])
# [1] 19581
# > 487512 + 19581
# [1] 507093

# head(rd[c(which(grepl("^sol", rd[,3]))),1][, c(which(grepl("^dilA_[23]", rd[2,])), 3)])



# This is the way!
# > sol_rows <- which(grepl("^sol", rd[,3]))
# > dil_cols <- which(grepl("^dilA_[23]", rd[2,]))
# head(rd[c(2, sol_rows),c(3, dil_cols)])