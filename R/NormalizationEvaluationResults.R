NormalizationEvaluationResults <- setClass("NormalizationEvaluationResults",
                                            slots=c(
                                                avgcvmem = "matrix",
                                                avgmadmem = "matrix",
                                                avgvarmem = "matrix",
                                                avgcvmempdiff = "numeric",
                                                avgmadmempdiff = "numeric",
                                                avgvarmempdiff = "numeric",
                                                nonsiganfdrlist = "numeric",
                                                nonsiganfdrlistcvpdiff = "numeric",
                                                anfdr = "matrix",
                                                kwfdr = "matrix"
                                            ),
                                            prototype=prototype())
                                            
