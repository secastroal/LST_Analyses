# Function to personalize some analysis options in the Mplus syntax ----

write.Mplus.options <- function(usevariables, cluster = NULL, analysis_type = c("GENERAL", "TWOLEVEL"),
                                  estimator = c("ML", "MLR", "BAYES"), iterations, h1iterations = NULL, 
                                  processors = NULL){
  usevariables <- paste("USEVAR", "=", paste(usevariables, collapse = "\n"), ";", sep = " ")
  if(!is.null(cluster)){
    cluster <- paste("CLUSTER", "=", cluster, ";", sep = " ")
  }
  analysis_type <- paste("TYPE", "=", analysis_type, ";", sep = " ")
  estimator <- paste("ESTIMATOR", "=", estimator, ";", sep = " ")
  if(estimator != "BAYES"){
    iterations <- paste("ITERATIONS", "=", iterations, ";", sep = " ")
    h1iterations <- paste("H1ITERATIONS", "=", h1iterations, ";", sep = " ")
    processors <- NULL
  }else{
    iterations <- paste("BITERATIONS", "=", paste0("(", iterations, ")"), ";", sep = " " )
    h1iterations <- NULL
    processors <- paste("PROCESSORS", "=", processors, sep = " ")
  }
  
  analysis_syntax <- paste(usevariables,
                           cluster,
                           "\nANALYSIS:",
                           analysis_type,
                           estimator,
                           iterations, 
                           h1iterations,
                           processors,
                           sep = "\n")
  cat(analysis_syntax)
  return(analysis_syntax)
}

