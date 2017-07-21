#Test Samara R

# Gestion des dates
lct <- Sys.getlocale("LC_TIME"); Sys.setlocale("LC_TIME", "C")

library(rgenoud)
library(recomeristem)

# setwd("D:/Modeles/Samara R")

meteo = recomeristem::getMeteo_from_files("D:/PAMStudio_dev/data/ecomeristem/ng")
paramInit = recomeristem::getParameters_from_files("D:/PAMStudio_dev/data/ecomeristem/ng")
vObs = recomeristem::rcpp_get_vObs("D:/PAMStudio_dev/data/ecomeristem/ng/vobs_G1_C_BFF2015.txt")
VarOfInterest <- c("pht","lig")
ParamList <- as.character(paramInit$Name)
paramInitTrans <- as.data.frame(matrix(ncol=length(ParamList)))
paramInitTrans[1,] <- paramInit$Values

resulttmp = recomeristem::rcpp_run_from_dataframe(paramInit, meteo)
obsRed = recomeristem::rcpp_reduceVobs(vObs, resulttmp)
resRed = recomeristem::rcpp_reduceResults(resulttmp, vObs)

ParamOfInterest <- c("coeff_lifespan","Epsib","density_IN", "Ict", "MGR_init", "plasto_init", "SLAp", "slope_length_IN", "slope_LL_BL_at_PI")
NbParam <- length(ParamOfInterest)
MinValue <- c(450, 5, 0.05, 1.8, 7, 30, 40, 0.2, 0.05)
MaxValue <- c(850, 6.5, 0.15, 2.5, 14, 50, 120, 0.8, 0.5)

#ParamOfInterest <- c("Ict")
#NbParam <- length(ParamOfInterest)
#MinValue <- c(3.0)
#MaxValue <- c(5.5)
Bounds <- matrix(c(MinValue,MaxValue),ncol=2)
Obs <- length(obsRed)

paramInit[which(paramInit=="nbleaf_pi"),2] = 18
paramInit[which(paramInit=="nb_leaf_stem_elong"),2] = 9
paramInit[which(paramInit=="phenostage_to_end_filling"),2] = 30
paramInit[which(paramInit=="phenostage_to_maturity"),2] = 31

param <- paramInit

Optim_Ecomeristem_funct <- function(p){
    error <- c()
    Res_ecomeristem <- recomeristem::rcpp_run_from_dataframe(param,meteo)
    res = recomeristem::rcpp_reduceResults(Res_ecomeristem, vObs)
    diff = obsRed - res
    return(sum(diff^2, na.rm=T))
}

resOptim <- genoud(Optim_Ecomeristem_funct, NbParam, max=FALSE, pop.size=200, max.generations=50000, wait.generations=100,
       hard.generation.limit=TRUE, starting.values=NULL, MemoryMatrix=TRUE,
       Domains=Bounds)


