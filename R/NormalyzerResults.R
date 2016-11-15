NormalyzerResults <- setClass("NormalyzerResults",
                              slots=c(
                                  data2log2 = "matrix",
                                  data2GI = "matrix",
                                  data2ctr = "matrix",
                                  data2med = "matrix",
                                  data2mean = "matrix",
                                  data2ctrlog = "matrix",
                                  data2vsn = "matrix",
                                  data2quantile = "matrix",
                                  data2mad = "matrix",
                                  data2loess = "matrix",
                                  data2vsnrep = "matrix",
                                  data2limloess = "matrix",
                                  
                                  globalfittedRLR = "matrix",
                                  fittedLR = "matrix"
                              ),
                              prototype=prototype())

