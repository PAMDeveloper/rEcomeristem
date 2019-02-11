#Script estimating parameters for EcoMeristem model
#Authors : Florian Larue, Gregory Beurier, Lauriane Rouan, Delphine Luquet
#-- (PAM, AGAP, BIOS, CIRAD)

###GET ARGUMENTS (for cluster)
# args = commandArgs(trailingOnly=TRUE)
# if (length(args)==0) {
#   stop("At least one argument must be supplied (path).n", call.=FALSE)
# }

###START INSTALL AND PACKAGES LOAD###
#.libPaths('/home/beurier/src/R/library')
library(recomeristem)
library(Rcpp)
library(parallel)
library(DEoptim)
library(data.table)

###SET INFORMATION FOR ESTIMATION###
FPath <- args[[1]]
nbCores <- 30
WorkDir <- "/homedir/beurier/work/"
DataDir <- paste(WorkDir,"phenoarch/rep",FPath,'/WW',sep="")
setwd(DataDir)
vName <- "vobs_moy"
paramOfInterest <- c("Epsib","Ict","MGR_init","plasto_init","phyllo_init","ligulo_init",
                     "coef_phyllo_PI","coef_ligulo_PI","leaf_length_to_IN_length","SLAp",
                     "coef_MGR_PI","slope_LL_BL_at_PI")
minValue <- c(6, 0.5, 6, 25, 25, 25, 1.0, 1.0, 0.1, 15, -0.5, 0.0)
maxValue <- c(20, 2.5, 14, 45, 45, 45, 3.0, 3.0, 0.2, 35, 0.5, 0.25)
coefIncrease <- 10
maxIter <- 2
relTol <- 0.01 #estimation stops if unable to reduce RMSE by (reltol * rmse) after steptol steps
stepTol <- 5000
solTol <- 0.0 #will be multiplied by the number of observed variables
clusterA <- TRUE  #parallel for machines with at least 4 cores

###FUNCTIONS###
optimEcomeristem <- function(p) {
  if(isInit == FALSE) {
    recomeristem::init_simu(param, meteo, obs, "env1")
    isInit <<- TRUE
  }
  #transform param over real interval
  p <- (p-bounds[,1])/(bounds[,2]-bounds[,1]) * (maxValue-minValue) + minValue

  if("phyllo_init" %in% paramOfInterest && "plasto_init" %in% paramOfInterest && "ligulo_init" %in% paramOfInterest) {
    if(p[match("phyllo_init",paramOfInterest)] < p[match("plasto_init",paramOfInterest)]) {
      return(99999)
    } else if(p[match("ligulo_init",paramOfInterest)] < p[match("phyllo_init",paramOfInterest)]) {
      return(99999)
    }
  }
  res <- recomeristem::launch_simu("env1", paramOfInterest, p)
  diff <- ((((obs - res)/obs)^2))*coeff
  diff <- sum(sqrt((colSums(diff, na.rm=T))/(colSums(!is.na(diff)))),na.rm=T)
  return(diff)
}
optimisation <- function(maxIter, solTol, relTol, stepTol, bounds) {
  if(clusterA) {
    resOptim <- DEoptim(optimEcomeristem, lower=bounds[,1], upper=bounds[,2], DEoptim.control(reltol=relTol,steptol=stepTol,VTR=solTol,itermax=maxIter,strategy=2,cluster=cl,packages=c("recomeristem"),parVar=c("meteo","obs", "paramOfInterest","coeff","isInit","param","bounds","minValue","maxValue")))
  } else {
    resOptim <- DEoptim(optimEcomeristem, lower=bounds[,1], upper=bounds[,2], DEoptim.control(reltol=relTol,steptol=stepTol,VTR=solTol,itermax=maxIter,strategy=2))
  }
  result$optimizer <- "Diffential Evolution Optimization"
  result$par <- resOptim$optim$bestmem
  result$value <- resOptim$optim$bestval
  result$iter <- resOptim$optim$iter
  return(list(result,resOptim))
}
coefCompute <- function(length,obsCoef) {
  tmp <- c()
  for(i in 0:(length-1)) {
    tmp <- c(tmp, (obsCoef+(i*(obsCoef/(100/coefIncrease)))))
  }
  return(data.frame(matrix(unlist(tmp),nrow=length,byrow=T)))
}
savePar <- function(pot) {
  resPar <<- matrix(as.vector(c(pot,result$value, result$par)), ncol=length(paramOfInterest)+2)
  write.table(resPar, file=paste("par_rep",FPath,".csv",sep=""), sep=",", append=T, dec=".",col.names=F,row.names = F)
}
res0 <- matrix(as.vector(c(NA,NA,minValue)), ncol=length(paramOfInterest)+2)
res0 <- rbind(res0,as.vector(c(NA,NA,maxValue)))
write.table(res0, file=paste("par_rep",FPath,".csv",sep=""), sep=",", append=F, dec=".",col.names=c("Plant","RMSE",paramOfInterest),row.names = F)
#For every plant in one rep
for(i in 32:33) {
  ###INIT SIMULATIONS###
  pot <- (FPath-1)*420+i
  traitment <- (((pot-1)%/%30)+1)%%2
  bordure <- min(pot%%30,pot%%31)
  if(traitment == 0 && bordure != 0 && pot != 1) {
    obsCSV <- read.table(paste(DataDir,"/",vName,"_P",pot,".txt",sep=""), header=T)
    meteo <- recomeristem::getMeteo_from_files(DataDir)
    param <- recomeristem::getParameters_from_files(DataDir)
    obs <- recomeristem::get_clean_obs(paste(DataDir,"/",vName,"_P",pot,".txt",sep=""))
    param[which(param$Name=="coef_lin_IN_diam"),"Values"] <- obsCSV[nrow(obsCSV),"stemdiam"]
    isInit <- FALSE
    minV <- rep(0,length(minValue))
    maxV <- rep(1,length(maxValue))
    bounds <- matrix(c(minV,maxV),ncol=2)
    obsCoef <- rep(1,ncol(obs))
    nbParam <- length(paramOfInterest)
    result <- list()
    VarList <- names(obs)

    #Optimisation run
    coeff <<- coefCompute(nrow(obs),obsCoef)
    set.seed(1337)
    if(clusterA) {
      nbCores <- nbCores
      cl <- makeCluster(nbCores)
      clusterEvalQ(cl, library(recomeristem, lib.loc='/homedir/beurier/src/R/library'))

    } #parallel
    system.time(resOptim <- optimisation(maxIter, solTol, relTol, stepTol, bounds))
    stopCluster(cl)
    result <- resOptim[[1]]
    result$par <- (result$par-bounds[,1])/(bounds[,2]-bounds[,1]) * (maxValue-minValue) + minValue
    savePar(pot)
  }
}
