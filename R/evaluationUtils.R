

# Input: Normalized matrices

# We want to target particular 7 vs 7 columns
# Also, we want to slice out potato-named rows into a distinctly separate dataframe

# Let's test this first in an R terminal


my_test <- function(target_matrix) {
    
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
}

