# Function to verify whether the model estimation in Mplus was terminated normally without any warning or error message.
# The input is an mplus.model object returned by the function MplusAutomation::readModels. 

check.mplus <- function(mplus.model, file.path){
  fit <- mplus.model
  H1check <- grep(pattern = "INCREASE THE NUMBER OF H1 ITERATIONS.", readLines(con = file.path))
  if(length(fit$parameters) == 0){
    return("Non-convergence")
  }else{
    if(!("est" %in% names(fit$parameters$unstandardized) && ("se" %in% names(fit$parameters$unstandardized) | "posterior_sd" %in% names(fit$parameters$unstandardized)))){
      return("Non-convergence")
    }else{
      if(length(fit$errors)!=0){
        if("USE THE FBITERATIONS OPTION TO INCREASE THE NUMBER OF ITERATIONS BY A FACTOR" %in% unlist(fit$errors) &&
           length(fit$errors) == 1){
          if(length(fit$warnings)!=0){
            if(paste0(length(fit$warnings)," WARNING(S) FOUND IN THE INPUT INSTRUCTIONS") %in% unlist(fit$warnings)){
              if(length(H1check>0)){
                return("Errors/Warnings")
              }else{
                return("Ok")
              }
            }else{
              return("Errors/Warnings")
            } 
          }else{
            if(length(H1check>0)){
              return("Errors/Warnings")
            }else{
              return("Ok")
            }
          }
        }else{
             return("Errors/Warnings")
           }
      }else{
        if(length(fit$warnings)!=0){
          if(paste0(length(fit$warnings)," WARNING(S) FOUND IN THE INPUT INSTRUCTIONS") %in% unlist(fit$warnings)){
            if(length(H1check>0)){
              return("Errors/Warnings")
            }else{
              return("Ok")
            }
          }else{
            return("Errors/Warnings")
          } 
        }else{
          if(length(H1check>0)){
            return("Errors/Warnings")
          }else{
            return("Ok")
          }
        }
      }
    }
  }
  
}


#files with Errors/Warnings to test

#fit.ok <- readModels(paste0(getwd(),"/Mplus_files_Results/OldAnalyses/cuts_long_n100_i4_nt100_old.out"))
#fit.warning.1 <- readModels(paste0(getwd(),"/Mplus_files/Test_analyses/cuts10.out")) 
#fit.warning.2 <- readModels(paste0(getwd(),"/Mplus_files/Test_analyses/cuts15.out")) 
#fit.warning.3 <- readModels(paste0(getwd(),"/Mplus_files/Test_analyses/tso4.out"))
#fit.warning.4 <- readModels(paste0(getwd(),"/Mplus_files/Test_analyses/tso60.out"))
#fit.warning.5 <- readModels(paste0(getwd(),"/ML_Mplus_files/mlcuts_m5.out")) 
#fit.warning.6 <- readModels(paste0(getwd(),"/Mplus_files_Results/OldAnalyses/cuts_trunc_n100_i4_nt5.out")) 
#fit.errors.1 <- readModels(paste0(getwd(),"/Mplus_files/Test_analyses/cuts25.out")) 
#fit.errors.2 <- readModels(paste0(getwd(),"/Mplus_files/Test_analyses/tso5.out"))
#fit.errors.3 <- readModels(paste0(getwd(),"/Mplus_files/Test_analyses/msst60.out"))
#fit.errors.4 <- readModels(paste0(getwd(),"/Mplus_files_Results/OldAnalyses/tso_long_n100_i4_nt10.out"))
#fit.fatal.1 <- readModels(paste0(getwd(),"/Mplus_files/Test_analyses/tso6.out"))
#fit.fatal.2 <- readModels(paste0(getwd(),"/Mplus_files_Results/OldAnalyses/cuts_wide_n100_i4_nt100_old.out"))
#fit.fatal.3 <- readModels(paste0(getwd(),"/Mplus_files_Results/OldAnalyses/cuts_wide_n100_i4_nt60_old.out"))
#fit.bayes.1 <- readModels(paste0(getwd(), "/Mplus_files_Results/TryBayes/tso_long_n100_i4_nt60_na0thin100.out"))
#fit.mlfatal <- readModels(paste0(getwd(), "/Mplus_files/Test_analyses/msst_long_bayes_10_34.out"))

#check.mplus(fit.ok, paste0(getwd(),"/Mplus_files_Results/OldAnalyses/cuts_long_n100_i4_nt100_old.out")) #ok
#check.mplus(fit.warning.1, paste0(getwd(),"/Mplus_files/Test_analyses/cuts10.out")) # W/E
#check.mplus(fit.warning.2, paste0(getwd(),"/Mplus_files/Test_analyses/cuts15.out")) # W/E
#check.mplus(fit.warning.3, paste0(getwd(),"/Mplus_files/Test_analyses/tso4.out")) # W/E
#check.mplus(fit.warning.4, paste0(getwd(),"/Mplus_files/Test_analyses/tso60.out")) #Non Convergence
#check.mplus(fit.warning.5, paste0(getwd(),"/ML_Mplus_files/mlcuts_m5.out")) #ok
#check.mplus(fit.warning.6, paste0(getwd(),"/Mplus_files_Results/OldAnalyses/cuts_trunc_n100_i4_nt5.out")) #ok
#check.mplus(fit.errors.1, paste0(getwd(),"/Mplus_files/Test_analyses/cuts25.out")) # W/E
#check.mplus(fit.errors.2, paste0(getwd(),"/Mplus_files/Test_analyses/tso5.out")) #Non Convergence
#check.mplus(fit.errors.3, paste0(getwd(),"/Mplus_files/Test_analyses/msst60.out")) # W/E
#check.mplus(fit.errors.4, paste0(getwd(),"/Mplus_files_Results/OldAnalyses/tso_long_n100_i4_nt10.out")) #Ok
#check.mplus(fit.fatal.1, paste0(getwd(),"/Mplus_files/Test_analyses/tso6.out")) #Non Convergence
#check.mplus(fit.fatal.2, paste0(getwd(),"/Mplus_files_Results/OldAnalyses/cuts_wide_n100_i4_nt100_old.out")) #Non-Convergence
#check.mplus(fit.fatal.3, paste0(getwd(),"/Mplus_files_Results/OldAnalyses/cuts_wide_n100_i4_nt60_old.out")) #Non-convergence
#check.mplus(fit.bayes.1, paste0(getwd(), "/Mplus_files_Results/TryBayes/tso_long_n100_i4_nt60_na0thin100.out")) #Ok
#check.mplus(fit.mlfatal, paste0(getwd(), "/Mplus_files/Test_analyses/msst_long_bayes_10_34.out")) #Non-Convergence
