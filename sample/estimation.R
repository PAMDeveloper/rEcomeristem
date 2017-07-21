#Test Samara R

# Gestion des dates
lct <- Sys.getlocale("LC_TIME"); Sys.setlocale("LC_TIME", "C")

library(rgenoud)
library(stringr)
library(recomeristem)

# setwd("D:/Modeles/Samara R")

meteo = recomeristem::getMeteo_from_files("D:/PAMStudio_dev/data/ecomeristem/ng")
paramInit = recomeristem::getParameters_from_files("D:/PAMStudio_dev/data/ecomeristem/ng")
vObs = recomeristem::rcpp_get_vObs("D:/PAMStudio_dev/data/ecomeristem/ng/vobs_G1_C_BFF2015.txt")
VarOfInterest <- c("pht","lig")

resulttmp = recomeristem::rcpp_run_from_dataframe(paramInit, meteo)
obsRed = recomeristem::rcpp_reduceVobs(vObs, resulttmp)
resRed = recomeristem::rcpp_reduceResults(resulttmp, vObs)

# Param?tres ? calibrer	Observations	V1	V2	V3
# LeafLengthMax (mm)	Longueur de la plus grande feuille	1000	1000	900
# CoefLeafWLRatio	Largeur/Longueur de la +G feuille	0.08	0.08	0.12
# CoefLeafDeath	Vari?t? STAY GREEN	-	-	0.002

ParamOfInterest <- c("coeff_lifespan","Epsib","density_IN", "Ict", "MGR_init", "plasto_init", "SLAp", "slope_length_IN", "slope_LL_BL_at_PI")
NbParam <- length(ParamOfInterest)
MinValue <- c(450, 5, 0.05, 1.8, 7, 30, 40, 0.2, 0.05)
MaxValue <- c(850, 6.5, 0.15, 2.5, 14, 50, 120, 0.8, 0.5)
Bounds <- matrix(c(MinValue,MaxValue),ncol=2)
Obs <- length(obsRed)

paramInit[which(paramInit=="nbleaf_pi"),2] = 18
paramInit[which(paramInit=="nb_leaf_stem_elong"),2] = 9
paramInit[which(paramInit=="phenostage_to_end_filling"),2] = 30
paramInit[which(paramInit=="phenostage_to_maturity"),2] = 31

param <- paramInit

Optim_Samara_funct <- function(p){
    error_v <- c()
    paramInitTrans[ParamOfInterest] <- p
    param <- data.frame(Name=ParamList,Values=unlist(paramInitTrans[1,]))
    Res_Samara <- rsamara::rcpp_run_from_dataframe(param,meteo)[[1]]
    print(Res_Samara)
    colnames(Res_Samara) <- str_to_lower(colnames(Res_Samara))
    Res_Samara$time <- as.Date(as.character(Res_Samara$time),format="%Y-%b-%d  %H:%M:%S")
    for (v in VarOfInterest){
      IdDate <- pmatch(as.character(ListObs[[v]][[2]]),as.character(Res_Samara$time))
      Simul_V_Date <- Res_Samara[IdDate,v]
      error_v <- c(error_v,ListObs[[v]][[1]]-Simul_V_Date)
    }
    return(sum(error_v^2))
  }

 test <- c(900,0.12,0.002)
#
#
 Optim_Samara_funct(test)




resOptim <- genoud(Optim_Samara_funct, NbParam, max=FALSE, pop.size=1000, max.generations=100, wait.generations=10,
       hard.generation.limit=TRUE, starting.values=NULL, MemoryMatrix=TRUE,
       Domains=Bounds)


