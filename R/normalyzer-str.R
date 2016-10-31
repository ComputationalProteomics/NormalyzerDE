


normalyzer<-function(datafile,getjob){
  #require(Rcmdr)
  #require(PerformanceAnalytics)
  #require(vsn)
  #require(preprocessCore)
  #require(limma)
  #require(MASS)
  #require(abind)
  #require(e1071)
  #require(ape)
  #require(raster)
  #require(car)
  #require(gridExtra)
  #require(ggplot2)
  print("Normalizing data....")
  try.result<-try(normalizeddata<-normMethods(datafile,getjob))
  if(inherits(try.result,"try-error")){
  return(try.result)
  }
print("Finished Normalization")
print("Analyzing data....")
  
try.result<-try(analyzeAndPlot(normalizeddata,getjob))
  if(inherits(try.result,"try-error")){
    return(try.result)
  }
print("Done .. Results are stored in the working directory")
}  



