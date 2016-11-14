

NormalyzerObject <- setClass("NormalyzerObject",
                             slots = c(rawData = "matrix",
                                       filterrawdata = "matrix",
                                       normfinderFilterRawData = "matrix",
                                       
                                       inputHeaderValues = "vector",
                                       sampleReplicateGroups = "vector",
                                       
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
                                       fittedLR = "matrix"))

## TestShow
## Template experimental method
setGeneric(name="testshow", def=function(object){})
setMethod(f="testshow", signature="NormalyzerObject",
    definition=function(object) {
        cat("This is a class-method testprint!")
    }
)






# This can be used by running:
# > testshow(obj)

setClass("BMI", representation(weight="numeric", size="numeric"))

setMethod("show", "BMI",
    function(object) {
        cat("BMI=", object@weight / (object@size ^ 2), " \n ")
    }          
)