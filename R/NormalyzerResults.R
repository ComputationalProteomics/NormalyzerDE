NormalyzerResults <- setClass("NormalyzerResults",
                              slots=c(
                                  nds = "NormalyzerDataset",
                                  
                                  methodnames = "character",
                                  normfinderMaxThreshold="numeric",
                                  globalNormalizationMinThreshold="numeric",
                                  
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
                              prototype=prototype(nds=NULL, normfinderMaxThreshold=1000, globalNormalizationMinThreshold=100))



setGeneric(name="performBasicNormalizations", function(nr) standardGeneric("performBasicNormalizations"))
setMethod("performBasicNormalizations", "NormalyzerResults",
          function(nr) {
              
              nds <- nr@nds

              # Only basic setup at this point?
              
              nr@data2log2 <- log2(nds@filterrawdata)
              nr@data2GI <- matrix(nrow=nrow(nds@filterrawdata), ncol=ncol(nds@filterrawdata), byrow=T)
              nr@data2ctr <- matrix(nrow=nrow(nds@filterrawdata), ncol=ncol(nds@filterrawdata), byrow=T)
              nr@data2med <- matrix(nrow=nrow(nds@filterrawdata), ncol=ncol(nds@filterrawdata), byrow=T) 
              nr@data2mean <- matrix(nrow=nrow(nds@filterrawdata), ncol=ncol(nds@filterrawdata), byrow=T)
              
              # print(nr)
              
              nr
          })