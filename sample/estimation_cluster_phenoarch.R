#Script estimating parameters for EcoMeristem model
#Authors : Florian Larue, Gregory Beurier, Lauriane Rouan, Delphine Luquet
#-- (PAM, AGAP, BIOS, CIRAD)

###GET ARGUMENTS (for cluster)
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (path).n", call.=FALSE)
}

###START INSTALL AND PACKAGES LOAD###
.libPaths('/home/beurier/src/R/library')
library(recomeristem)
library(Rcpp)
library(parallel)
library(DEoptim)
library(data.table)

###SET INFORMATION FOR ESTIMATION###
FPath <- as.numeric(args[[1]])
nbCores <- 30
WorkDir <- "/homedir/beurier/work/"
DataDir <- paste0(WorkDir,"phenoarch/rep",FPath,'/WW')
MeteoDir <- paste0(WorkDir,"phenoarch/meteo_local")
setwd(DataDir)
vName <- "vobs_moy"
paramOfInterest <- c("Epsib","Ict","MGR_init","plasto_init","phyllo_init","ligulo_init",
                     "leaf_length_to_IN_length","SLAp","gdw")
minValue <- c(6, 0.5, 6, 30.0, 30.0, 30.0, 0.1, 35, 0.03)
maxValue <- c(20, 2.5, 14, 40.0, 40.0, 45.0, 0.2, 60, 1.0)
coefIncrease <- 0
maxIter <- 30
nbplants <- 35
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
  diff <- sum((colSums(diff, na.rm=T))/(colSums(!is.na(diff))),na.rm=T)
  return(diff)
}
optimisation <- function(maxIter, solTol, relTol, stepTol, bounds) {
  if(clusterA) {
    resOptim <- DEoptim(optimEcomeristem, lower=bounds[,1], upper=bounds[,2], DEoptim.control(itermax=maxIter,strategy=2,cluster=cl,packages=c("recomeristem"),parVar=c("meteo","obs", "paramOfInterest","coeff","isInit","param","bounds","minValue","maxValue")))
  } else {
    resOptim <- DEoptim(optimEcomeristem, lower=bounds[,1], upper=bounds[,2], DEoptim.control(itermax=maxIter,strategy=2))
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

stemdiamprob <- c()
#For every plant in one rep
for(i in 1:nbplants) {
  ###INIT SIMULATIONS###
  pot <- (FPath-1)*420+i
  traitment <- (((pot-1)%/%30)+1)%%2
  bordure <- min(pot%%30,(pot-1)%%30)
  if(traitment == 0 && bordure != 0 && pot != 1) {
    obsCSV <- read.table(paste0(DataDir,"/",vName,"_P",pot,".txt"), header=T)

    meteoCSV <- read.table(paste0(MeteoDir,"/plant_",i,".csv"),sep=";", header=T)
    meteoCSV <- meteoCSV[!duplicated(meteoCSV$date),c("Temperatured","PARd","ETPd","irrigd","P")]
    meteoCSV <- meteoCSV[-nrow(meteoCSV),]
    row.names(meteoCSV) <- NULL
    colnames(meteoCSV) <- c("Temperature","Par","Etp","Irrigation","P")
    meteo <- meteoCSV

    param <- recomeristem::getParameters_from_files(DataDir)
    obs <- recomeristem::get_clean_obs(paste(DataDir,"/",vName,"_P",pot,".txt",sep=""))
    if(!is.na(obsCSV[nrow(obsCSV),"stemdiam"])) {
      param[which(param$Name=="coef_lin_IN_diam"),"Values"] <- obsCSV[nrow(obsCSV),"stemdiam"]
    } else {
      param[which(param$Name=="coef_lin_IN_diam"),"Values"] <- 1.5
      stemdiamprob <- c(stemdiamprob,i)
    }
    param[which(param$Name=="RU1"),"Values"] <- obsCSV[nrow(obsCSV),"RU1"]
    param[which(param$Name=="pf"),"Values"] <- obsCSV[nrow(obsCSV),"PF"]
    param[which(param$Name=="swc_init"),"Values"] <- obsCSV[1,"swc"]

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
