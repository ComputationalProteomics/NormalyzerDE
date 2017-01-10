# This is an R-version of Normfinder including only the simplest 
# output and for the case of more than one group.

# Input data is a space seperated text file.
# The first entry is the number of genes in the experiment. 
# Then follows row by row the pcr values for all the samples, each 
# row is the values for one gene.
# Finally there is a line giving the group labels: these must be 
# integers.
# The input file is not allowed to contain further material.

# Output: In the output genes are numbered consequtively following 
# the rows of the input. The first results are the gene numbers and 
# the quality measures. Next follows the same values but arranged 
# in increasing order of the quality measure. 
# The gene with the smallest value of the quality measure is selected 
# and combined with each of the other genes to produce a quality 
# measure for the combination of two genes. These quality measures 
# are then listed. The entry given for the gene that was selected in 
# the first round is the original quality measure for that gene 
# multiplied by the square root of k/(k-1), where k is the number 
# of genes.

# THIS SCRIPT FOR NORMFINDER IS MODIFIED FOR USE IN NORMALYZER. BY USING 
# NORMALYZER YOU AGREE TO ACCEPT THE ORIGINAL TERMS OF USE FOR NORMFINDER. 
# TERMS AND CONDITIONS OF USE FOR NORMFINDER CAN BE FOUND IN THIS LINK: 
# http://moma.dk/normfinder-software

normfinder <- function(nds) {

    filterrawdata1 <- nds@normfinderFilterRawData
    getEDdata <- nds@inputHeaderValues
    
    gr <- getEDdata[-which(as.numeric(getEDdata) < 1)]
    
    filnavn <- filterrawdata1[, which(as.numeric(filterrawdata1[1,]) > 0)]
    filnavn <- as.matrix(filnavn[-(1:2), ])
    filnavn <- apply(filnavn, 2, function(x) as.numeric(gsub(",", ".", x)))
    class(filnavn) <- "numeric"

    k <- nrow(filnavn)
    gr <- unlist(gr)
    da <- as.matrix(filnavn)
    m <- length(levels(as.factor(gr)))
    medgr <- as.numeric(levels(as.factor(unlist(gr))))
    
    medgen <- c(1:k)                # all genes are included in the analysis
    
    y1 <- log(da)                   # log of pcr values

    y <- y1[medgen, gr == medgr[1]]
    for (i in 2:m) {
        y <- cbind(y, y1[medgen, gr == medgr[i]])
    }                               # data rearranged according to group

    ngr <- rep(0, m)                # number of samples in each group
    for (i in 1:m) {
        ngr[i] <- sum(gr == medgr[i])
    }
    grny <- rep(c(1:m), ngr)        # group labels for data in y
    n <- sum(ngr)

    mei <- apply(y, 1, mean)
    mej <- apply(y, 2, mean)
    me <- mean(mej)
    
    # Estimates of the variances when group label is not considered
    
    a <- rep(0, k)
    for (i in 1:k) {
        a[i] <- sum((y[i, ] - mej - mei[i] + me)^2) / (n - 1)}
    
    b <- sum(a)
    varnogroup <- (a - b / (k * k - k)) / (1 - 2 / k)
    
    # Estimates of variances for the m groups and k genes
    
    meigr <- matrix(rep(0, k*m), k, m)
    for (j in 1:m){
        meigr[,j] <- apply(y[, grny == j], 1, mean)
    }
    
    megr <- rep(0, m)
    for (j in 1:m) {
        megr[j] <- mean(meigr[, j])
    }
    
    g <- y
    for (j in 1:n) {
        g[, j] <- y[, j] - meigr[, grny[j]] - mej[j] + megr[grny[j]]
    }
    
    vargroupall <- matrix(rep(0, m*k), m, k)
    for (j in 1:m) {
        a <- rep(0, k)
        for (i in 1:k) {
            a[i] <- sum((g[i, grny == j])^2) / (ngr[j] - 1)
        }
        b <- sum(a)
        vargroupall[j, ] <- (a - b / (k * k - k)) / (1 - 2 / k)
    }
    
    varmin <- vargroupall#takes time
    for (i in 1:m) {
        z <- y[, grny == i]
        for (j in 1:k){
            varpair <- rep(0, k)
            for (j1 in 1:k) {
                varpair[j1] <- stats::var(z[j, ] - z[j1, ])
            }
            varmin[i,j] <- min(varpair[varpair > 0]) / 4
        }
    }
    
    #intragroupvariation
    vargroupall=ifelse(vargroupall<0,varmin,vargroupall)
    
    # Variances have been estimated and corrected if negative values
    # were encountered
    
    # quality measure for each gene is calculated
    
    dif <- meigr
    m1i <- apply(dif, 1, mean)
    m1j <- apply(dif, 2, mean)
    m1 <- mean(m1i)
    
    #intergroup variation
    for (i in 1:k) {
        for (j in 1:m) {
            dif[i, j] <- dif[i, j] - m1i[i] - m1j[j] + m1
        }
    }
     
    va <- vargroupall
    for (j in 1:m) {
        va[j,] <- va[j,] / ngr[j]}
    
    tau <- -1
    
    tau <- sum(dif * dif) / ((m-1) * (k-1)) - mean(va)

    if (tau < 0) {
        tau <- 0
    }
    
    dnew <- dif * tau / (tau + t(va))
    vanew <- t(va + tau * va / (tau + va))
    
    
    qm=abs(dnew)+sqrt(vanew)
    qmaal=apply(qm,1,mean)
    
    # We now look for the best combination of the gene with the smallest 
    # value of the quality measure and one more gene
    Hkg <- NULL
    G1 <- 1
    G2 <- 1
    var <- 1000
    k2 <- k
    
    for(count in 1:k2) {
        b <- order(qmaal)[count]
        qmaaldob <- rep(0, k2)
        qmaaldob[b] <- qmaal[b] * sqrt(k2 / (k2 - 1))
        
        for (j in c(1:k)[-b]) {
            a <- c(b,j)
            a1 <- dnew[a,]
            a2 <- apply(a1, 2, mean) * sqrt(k / (k - 2))
            b1 <- vanew[a,]
            b2 <- apply(b1, 2, mean) / 2
            qmaaldob[j] <- mean(abs(a2) + sqrt(b2))
        }
        
        if (var > qmaaldob[order(qmaaldob)][1]) {
            var <- qmaaldob[order(qmaaldob)][1]
            G1 <- b
            G2 <- order(qmaaldob)[1]
        }
    }
    
    Hkg <- c(G1, G2)
    
    return(filterrawdata1[(Hkg + 2), ])
}


