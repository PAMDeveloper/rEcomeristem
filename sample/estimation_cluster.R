#Script estimating parameters for EcoMeristem model
#Authors : Florian Larue, Gregory Beurier, Lauriane Rouan, Delphine Luquet
#-- (PAM, AGAP, BIOS, CIRAD)
###Set informations for parameter estimation###
FPath <- ""
VName <- "vobs_moy.txt"
VECName <- "vobs_et.txt"
#Temoin
ParamOfInterest <- c("Epsib", "Ict", "MGR_init", "plasto_init", "phyllo_init", "ligulo_init", "coef_MGR_PI", "slope_length_IN", "density_IN2", "coef_plasto_PI", "coef_phyllo_PI", "coef_ligulo_PI", "slope_LL_BL_at_PI")
MinValue <- c(3, 0.5, 6, 15, 20, 20, -0.5, 0.5, 0.08, 1, 1, 1, 0.0)
MaxValue <- c(8, 2.5, 14, 45, 45, 45, 0.5, 1, 0.3, 3.0, 3.0, 3.0, 0.4)
coefIncrease <- 10
MaxIter <- 3000
Penalty <- 10 #Penalty for simulation outside of SD (RMSE * Penalty)
SolTol <- 0.01 #will be multiplied by the number of observed variables
ACluster <- TRUE  #parallel for machines with at least 4 cores
coresNb <- 48
###End informations for parameter estimation###

#Install and load packages
list.of.packages <- c("parallel","DEoptim", "truncnorm","data.table")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, library, character.only=TRUE)
library("recomeristem")


#Get data for simulation
setwd(FPath)
meteo <- recomeristem::getMeteo_from_files(FPath)
paramInit <- recomeristem::getParameters_from_files(FPath)
vObs <- recomeristem::rcpp_get_vObs(paste(FPath,paste("/",VName, sep=""), sep=""))
obsET <- recomeristem::rcpp_get_vObs(paste(FPath,paste("/",VECName, sep=""), sep=""))
ParamList <- as.character(paramInit$Name)
paramInitTrans <- as.data.frame(matrix(ncol=length(ParamList)))
colnames(paramInitTrans) <- ParamList
paramInitTrans[1,] <- paramInit$Values
NbParam <- length(ParamOfInterest)
Bounds <- matrix(c(MinValue,MaxValue),ncol=2)

#First model run and matrix reduction
resulttmp <- recomeristem::rcpp_run_from_dataframe(paramInit, meteo)
obsRed <- recomeristem::rcpp_reduceVobs(vObs, resulttmp)
obsETRed <- recomeristem::rcpp_reduceVobs(obsET, resulttmp)
resRed <- recomeristem::rcpp_reduceResults(resulttmp, vObs)
VarList <- names(obsRed)
obsLength <- nrow(obsRed)
coeff <- c()
obsCoef <- rep(1,ncol(obsRed))
SolTol <- SolTol * length(VarList)
res <- list()

#Functions
optimEcomeristem <- function(p){
  if("phyllo_init" %in% ParamOfInterest && "plasto_init" %in% ParamOfInterest && "ligulo_init" %in% ParamOfInterest) {
    if(p[match("phyllo_init",ParamOfInterest)] < p[match("plasto_init",ParamOfInterest)]) {
      return(99999)
    } else if(p[match("ligulo_init",ParamOfInterest)] < p[match("phyllo_init",ParamOfInterest)]) {
      return(99999)
    } else if(p[match("phyllo_init",ParamOfInterest)]*p[match("coef_phyllo_PI",ParamOfInterest)] < p[match("plasto_init",ParamOfInterest)]*p[match("coef_plasto_PI",ParamOfInterest)]) {
      return(99999)
    } else if(p[match("ligulo_init",ParamOfInterest)]*p[match("coef_ligulo_PI",ParamOfInterest)] < p[match("phyllo_init",ParamOfInterest)]*p[match("coef_phyllo_PI",ParamOfInterest)]) {
      return(99999)
    }
  }
  paramInitTrans[ParamOfInterest] <- p
  parameters <- data.frame(Name=ParamList, Values=unlist(paramInitTrans[1,]))
  Res_ecomeristem <- recomeristem::rcpp_run_from_dataframe(parameters,meteo)
  res <- recomeristem::rcpp_reduceResults(Res_ecomeristem, vObs)
  diff1 <- abs(1-(data.table::between(res,obsRed-obsETRed,obsRed+obsETRed)))
  diff2 <- diff1 * Penalty
  diff3 <-  replace(diff2, diff2 == 0,1)
  diff <- ((((obsRed - res)/obsRed)^2)*diff3)*coeff
  return(sum(sqrt((colSums(diff, na.rm=T))/(colSums(!is.na(diff)))),na.rm=T))
}
optimisation <- function(Optimizer, MaxIter, SolTol, NbParam, Bounds, NbEnv, SDate, EDate) {
  coeff <<- coefCompute()
  if(ACluster) {
    resOptim <- DEoptim(optimEcomeristem, lower=Bounds[,1], upper=Bounds[,2], DEoptim.control(VTR=SolTol,itermax=MaxIter, strategy=2, cluster=cl, packages=c("recomeristem"), parVar=c("meteo","vObs","obsRed", "paramInitTrans", "ParamOfInterest", "ParamList", "obsETRed","obsCoef","Penalty","coeff")))
  } else {
    resOptim <- DEoptim(optimEcomeristem, lower=Bounds[,1], upper=Bounds[,2], DEoptim.control(VTR=SolTol,itermax=MaxIter, strategy=2))
  }
  res$par <- resOptim$optim$bestmem
  return(list(res,resOptim))
}
coefCompute <- function() {
  tmp <- c()
  for(i in 0:(obsLength-1)) {
    tmp <- c(tmp, (obsCoef+(i*(obsCoef/(100/coefIncrease)))))
  }
  return(data.frame(matrix(unlist(tmp),nrow=obsLength,byrow=T)))
}
savePar <- function(name = Sys.time()) {
  resPar <- matrix(as.vector(res$par), ncol=length(ParamOfInterest))
  write.table(resPar, file="par.csv", sep=",", append=T, dec=".",col.names = F,row.names = name)
}

#Optimisation run
set.seed(224)
if(ACluster) {
  nbCores <- coresNb
  cl <- makeCluster(nbCores, outfile="clusterlog.txt")
  clusterEvalQ(cl, library(recomeristem,rgenoud))
  clusterExport(cl, varlist = c("meteo","vObs","obsRed", "paramInitTrans", "ParamOfInterest", "ParamList", "obsETRed", "obsCoef","Penalty","coeff"))
} #parallel
resOptim <- optimisation(Optimizer, MaxIter, SolTol, NbParam, Bounds, NbEnv, SDate, EDate)
res <- resOptim[[1]]
savePar()
stopCluster(cl)
