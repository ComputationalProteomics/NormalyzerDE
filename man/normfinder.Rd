\name{normfinder}
\alias{normfinder}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
NormFinder
}
\description{
NormFinder
}
\usage{
normfinder(filterrawdata1, getEDdata)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filterrawdata1}{
%%     ~~Describe \code{filterrawdata1} here~~
}
  \item{getEDdata}{
%%     ~~Describe \code{getEDdata} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (filterrawdata1, getEDdata) 
{
    gr <- getEDdata[-which(getEDdata < 1)]
    filnavn <- filterrawdata1[, -(1:(length(getEDdata) - length(gr)))]
    filnavn <- (as.matrix((filnavn[-(1:2), ])))
    class(filnavn) <- "numeric"
    k = nrow(filnavn)
    gr <- unlist(gr)
    da <- as.matrix(filnavn)
    m <- length(levels(as.factor(gr)))
    medgr <- as.numeric(levels(as.factor(unlist(gr))))
    medgen = c(1:k)
    y1 = log(da)
    y = y1[medgen, gr == (medgr[1])]
    for (i in 2:m) {
        y = cbind(y, y1[medgen, gr == (medgr[i])])
    }
    ngr = rep(0, m)
    for (i in 1:m) {
        ngr[i] = sum(gr == (medgr[i]))
    }
    grny = rep(c(1:m), ngr)
    n = sum(ngr)
    mei = apply(y, 1, mean)
    mej = apply(y, 2, mean)
    me = mean(mej)
    a = rep(0, k)
    for (i in 1:k) {
        a[i] = sum((y[i, ] - mej - mei[i] + me)^2)/(n - 1)
    }
    b = sum(a)
    varnogroup = (a - b/(k * k - k))/(1 - 2/k)
    meigr = matrix(rep(0, k * m), k, m)
    for (j in 1:m) {
        meigr[, j] = apply(y[, grny == j], 1, mean)
    }
    megr = rep(0, m)
    for (j in 1:m) {
        megr[j] = mean(meigr[, j])
    }
    g = y
    for (j in 1:n) {
        g[, j] = y[, j] - meigr[, grny[j]] - mej[j] + megr[grny[j]]
    }
    vargroupall = matrix(rep(0, m * k), m, k)
    for (j in 1:m) {
        a = rep(0, k)
        for (i in 1:k) {
            a[i] = sum((g[i, grny == j])^2)/(ngr[j] - 1)
        }
        b = sum(a)
        vargroupall[j, ] = (a - b/(k * k - k))/(1 - 2/k)
    }
    varmin = vargroupall
    for (i in 1:m) {
        z = y[, grny == i]
        for (j in 1:k) {
            varpair = rep(0, k)
            for (j1 in 1:k) {
                varpair[j1] = var(z[j, ] - z[j1, ])
            }
            varmin[i, j] = min(varpair[varpair > 0])/4
        }
    }
    vargroupall = ifelse(vargroupall < 0, varmin, vargroupall)
    dif = meigr
    m1i = apply(dif, 1, mean)
    m1j = apply(dif, 2, mean)
    m1 = mean(m1i)
    for (i in 1:k) {
        for (j in 1:m) {
            dif[i, j] = dif[i, j] - m1i[i] - m1j[j] + m1
        }
    }
    va = vargroupall
    for (j in 1:m) {
        va[j, ] = va[j, ]/ngr[j]
    }
    tau = -1
    tau = sum(dif * dif)/((m - 1) * (k - 1)) - mean(va)
    if (tau < 0) {
        tau = 0
    }
    dnew = dif * tau/(tau + t(va))
    vanew = t(va + tau * va/(tau + va))
    qm = abs(dnew) + sqrt(vanew)
    qmaal = apply(qm, 1, mean)
    Hkg <- NULL
    G1 <- 1
    G2 <- 1
    var <- 1000
    k2 <- k
    for (count in 1:k2) {
        b = order(qmaal)[count]
        qmaaldob = rep(0, k2)
        qmaaldob[b] = qmaal[b] * sqrt(k2/(k2 - 1))
        for (j in c(1:k)[-b]) {
            a = c(b, j)
            a1 = dnew[a, ]
            a2 = apply(a1, 2, mean) * sqrt(k/(k - 2))
            b1 = vanew[a, ]
            b2 = apply(b1, 2, mean)/2
            qmaaldob[j] = mean(abs(a2) + sqrt(b2))
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
