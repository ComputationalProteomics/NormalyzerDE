

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

full_report <- "FeatureReport_BA_MQP_20140828-fornormalyzer_annotated.txt"
subset_report <- "BA_MQP.subset_500.tsv"

run_review_data_test <- function(subset=FALSE) {
    
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
    potato_sub <- get_normalyzer_df_subset(nr, "^dilA_[23]", "^sol", 3, inverse_name_pattern=FALSE)
    non_potato_sub <- get_normalyzer_df_subset(nr, "^dilA_[23]", "^sol", 3, inverse_name_pattern=TRUE)

    
    
    # nr <- analyzeNormalizations(nr, jobName)
    
    end_time <- Sys.time()
    print(difftime(end_time, start_time))
    print(paste("Start:", start_time, "End:", end_time))
}


get_normalyzer_df_subset <- function(nr, col_pattern, row_pattern, row_name_col, 
                                     inverse_name_pattern=FALSE) {
    
    raw_data <- nr@nds@rawData
    df <- nr@data2log2
    
    header_row <- 2
    
    if (inverse_name_pattern) {
        target_rows <- which(!grepl(row_pattern, raw_data[-1:-2, row_name_col]))
    }
    else {
        target_rows <- which(grepl(row_pattern, raw_data[-1:-2, row_name_col]))
    }

    target_cols <- which(grepl(col_pattern, raw_data[header_row,]))
    
    print(target_rows)
    
    df[target_rows, target_cols]
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