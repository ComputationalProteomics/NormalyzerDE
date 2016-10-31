\name{normalyzer}
\alias{normalyzer}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Normalyzer
}
\description{
Normalyzer: a tool for rapid evaluation of normalization methods for omics datasets
}
\usage{
normalyzer(file, name)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{
  File could be a tab separated text file or a dataframe or a matrix.
}
  \item{name}{
  Name of the project
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
function (file, name) 
{
    source("normMethods.R")
    normalizeddata <- normMethods(file, name)
    print("Finished Normalizaiton")
    source("analyzeAndPlot.R")
    analyzeAndPlot(normalizeddata, name)
    print(paste("Done .. Results are stored in '", getwd(), "/", 
        name, "' folder", sep = ""))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
