#Script estimating parameters for EcoMeristem model
#Authors : Florian Larue, Gregory Beurier, Lauriane Rouan, Delphine Luquet
#-- (PAM, AGAP, BIOS, CIRAD)

###SET INFORMATION FOR ESTIMATION###
path <- "D:/Workspace/estimworkspace/Sorghum/2017/phenoarch/G23"
vName <- "vobs_moy.txt"
paramOfInterest <- c("Epsib","Ict","MGR_init","plasto_init","phyllo_init","ligulo_init","density_IN2","coef_plasto_PI",
                     "coef_phyllo_PI","coef_ligulo_PI","leaf_length_to_IN_length")
minValue <- c(6, 0.5, 6, 25, 25, 25, 0.01, 1.0, 1.0, 1.0, 0.1)
maxValue <- c(20, 2.5, 14, 45, 45, 45, 0.3, 3.0, 3.0, 3.0, 0.2)
coefIncrease <- 2
maxIter <- 10000
solTol <- 0.01 #will be multiplied by the number of observed variables
relTol <- 0.001 #estimation stops if unable to reduce RMSE by (reltol * rmse) after steptol steps
stepTol <- maxIter #see above
clusterA <- TRUE  #parallel for machines with at least 4 cores
pots <- c(154,478,938,1672)

###START INSTALL AND PACKAGES LOAD###
list.of.packages <- c("Rcpp","recomeristem","parallel","DEoptim")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
invisible(lapply(list.of.packages, library, character.only=TRUE))

###INIT SIMULATIONS###
setwd(paste(path,"/rep2/",sep=""))
meteo1 <- recomeristem::getMeteo_from_files(paste(path,"/rep1",sep=""))
meteo2 <- recomeristem::getMeteo_from_files(paste(path,"/rep2",sep=""))
meteo3 <- recomeristem::getMeteo_from_files(paste(path,"/rep3",sep=""))
meteo4 <- recomeristem::getMeteo_from_files(paste(path,"/rep4",sep=""))

param1 <- recomeristem::getParameters_from_files(paste(path,"/rep1",sep=""))
param2 <- recomeristem::getParameters_from_files(paste(path,"/rep2",sep=""))
param3 <- recomeristem::getParameters_from_files(paste(path,"/rep3",sep=""))
param4 <- recomeristem::getParameters_from_files(paste(path,"/rep4",sep=""))

obs1 <- recomeristem::get_clean_obs(paste(path,"/rep1/",vName,sep=""))
obs2 <- recomeristem::get_clean_obs(paste(path,"/rep2/",vName,sep=""))
obs3 <- recomeristem::get_clean_obs(paste(path,"/rep3/",vName,sep=""))
obs4 <- recomeristem::get_clean_obs(paste(path,"/rep4/",vName,sep=""))

obsCoef1 <- rep(1,ncol(obs1))
obsCoef2 <- rep(1,ncol(obs2))
obsCoef3 <- rep(1,ncol(obs3))
obsCoef4 <- rep(1,ncol(obs4))

###INIT PARAMFILE FOR SAVING###
paramInit <- recomeristem::getParameters_from_files(paste(path,"/rep3",sep=""))
ParamList <- as.character(paramInit$Name)
paramInitTrans <- as.data.frame(matrix(ncol=length(ParamList)))
colnames(paramInitTrans) <- ParamList
paramInitTrans[1,] <- paramInit$Values

###INIT ESTIMATION###
isInit <- FALSE
bounds <- matrix(c(minValue,maxValue),ncol=2)
nbParam <- length(paramOfInterest)
result <- list()
VarList <- names(obs1)

###FUNCTIONS###
optimEcomeristem <- function(p) {
  if(isInit == FALSE) {
    recomeristem::init_simu(param1, meteo1, obs1, "env1")
    recomeristem::init_simu(param2, meteo2, obs2, "env2")
    recomeristem::init_simu(param3, meteo3, obs3, "env3")
    recomeristem::init_simu(param4, meteo4, obs4, "env4")
    isInit <<- TRUE
  }
  if("phyllo_init" %in% paramOfInterest && "plasto_init" %in% paramOfInterest && "ligulo_init" %in% paramOfInterest) {
    if(p[match("phyllo_init",paramOfInterest)] < p[match("plasto_init",paramOfInterest)]) {
      return(99999)
    } else if(p[match("ligulo_init",paramOfInterest)] < p[match("phyllo_init",paramOfInterest)]) {
      return(99999)
    } else if(p[match("phyllo_init",paramOfInterest)]*p[match("coef_phyllo_PI",paramOfInterest)] < p[match("plasto_init",paramOfInterest)]*p[match("coef_plasto_PI",paramOfInterest)]) {
      return(99999)
    } else if(p[match("ligulo_init",paramOfInterest)]*p[match("coef_ligulo_PI",paramOfInterest)] < p[match("phyllo_init",paramOfInterest)]*p[match("coef_phyllo_PI",paramOfInterest)]) {
      return(99999)
    }
  }
  #res1 <- recomeristem::launch_simu("env1", paramOfInterest, p)
  res2 <- recomeristem::launch_simu("env2", paramOfInterest, p)
  #res3 <- recomeristem::launch_simu("env3", paramOfInterest, p)
  #res4 <- recomeristem::launch_simu("env4", paramOfInterest, p)

  #diff1 <- ((((obs1 - res1)/obs1)^2))*coeff1
  #diff1 <- sum(sqrt((colSums(diff1, na.rm=T))/(colSums(!is.na(diff1)))),na.rm=T)

  diff2 <- ((((obs2 - res2)/obs2)^2))*coeff2
  diff2 <- sum(sqrt((colSums(diff2, na.rm=T))/(colSums(!is.na(diff2)))),na.rm=T)

  #diff3 <- ((((obs3 - res3)/obs3)^2))*coeff3
  #diff3 <- sum(sqrt((colSums(diff3, na.rm=T))/(colSums(!is.na(diff3)))),na.rm=T)

  #diff4 <- ((((obs4 - res4)/obs4)^2))*coeff4
  #diff4 <- sum(sqrt((colSums(diff4, na.rm=T))/(colSums(!is.na(diff4)))),na.rm=T)

  return(diff2)
}
optimisation <- function(Optimizer, maxIter, solTol, bounds) {
  if(clusterA && detectCores() >= 4) {
    resOptim <- DEoptim(optimEcomeristem, lower=bounds[,1], upper=bounds[,2], DEoptim.control(reltol=relTol,steptol=stepTol,VTR=solTol,itermax=maxIter,strategy=2,cluster=cl,packages=c("recomeristem"),parVar=c("meteo1","meteo2","meteo3","meteo4","obs1","obs2","obs3","obs4","param1","param2","param3","param4","coeff1","coeff2","coeff3","coeff4","paramOfInterest","isInit")))
  } else {
    resOptim <- DEoptim(optimEcomeristem, lower=bounds[,1], upper=bounds[,2], DEoptim.control(reltol=relTol,steptol=stepTol,VTR=solTol,itermax=maxIter,strategy=2))
  }
  result$optimizer <- "Diffential Evolution Optimization"
  result$par <- resOptim$optim$bestmem
  result$value <- resOptim$optim$bestval
  result$iter <- resOptim$optim$iter
  return(list(result,resOptim))
}
resPlot <- function(obs,meteo,param) {
  VarList <- names(obs)
  p <- result$par
  recomeristem::init_simu(param, meteo, obs, "resP")
  res <<- recomeristem::launch_simu("resP", paramOfInterest, p)
  obsT <- data.frame(day = seq(1,max(obs$day,na.rm=T)))
  for(v in VarList) {
    if(v != "day") {
      obsT[,v] <- 1
    }
  }
  recomeristem::init_simu(param, meteo, obsT, "total_res")
  resT <<- recomeristem::launch_simu("total_res", paramOfInterest, p)

  plotF <- function(x) {
    plot(resT[[x]], type="l", pch=3, xlab="DAS", ylab=x,xaxs="i", yaxs="i", ylim=c(0,max(max(resT[[x]],na.rm=T),max(obs[[x]],na.rm=T),na.rm=T)*1.1), xlim=c(0,max(obs$day,na.rm=T)+5))
    points(obs$day, obs[[x]], type="p", col="red")
    diff <- ((obs[[x]] - res[[x]])/obs[[x]])^2
    return(sqrt((sum(diff, na.rm=T))/(sum(!is.na(diff)))))
  }
  VarListP <- VarList[VarList != "day"]
  sapply(VarListP, plotF)
}
resAPlot <- function() {
  p <- result$par
  recomeristem::init_simu(param1, meteo1, obs1, "res1")
  recomeristem::init_simu(param2, meteo2, obs2, "res2")
  recomeristem::init_simu(param3, meteo3, obs3, "res3")
  recomeristem::init_simu(param4, meteo4, obs4, "res4")
  resP1 <<- recomeristem::launch_simu("res1", paramOfInterest, p)
  resP2 <<- recomeristem::launch_simu("res2", paramOfInterest, p)
  resP3 <<- recomeristem::launch_simu("res3", paramOfInterest, p)
  resP4 <<- recomeristem::launch_simu("res4", paramOfInterest, p)
  obsT1 <- data.frame(day = seq(1,max(obs1$day,na.rm=T)))
  obsT2 <- data.frame(day = seq(1,max(obs2$day,na.rm=T)))
  obsT3 <- data.frame(day = seq(1,max(obs3$day,na.rm=T)))
  obsT4 <- data.frame(day = seq(1,max(obs4$day,na.rm=T)))
  for(v in VarList) {
    if(v != "day") {
      obsT1[,v] <- 1
      obsT2[,v] <- 1
      obsT3[,v] <- 1
      obsT4[,v] <- 1
    }
  }
  recomeristem::init_simu(param1, meteo1, obsT1, "total_res1")
  resT1 <<- recomeristem::launch_simu("total_res1", paramOfInterest, p)
  recomeristem::init_simu(param2, meteo2, obsT2, "total_res2")
  resT2 <<- recomeristem::launch_simu("total_res2", paramOfInterest, p)
  recomeristem::init_simu(param3, meteo3, obsT3, "total_res3")
  resT3 <<- recomeristem::launch_simu("total_res3", paramOfInterest, p)
  recomeristem::init_simu(param4, meteo4, obsT4, "total_res4")
  resT4 <<- recomeristem::launch_simu("total_res4", paramOfInterest, p)
  plotF <- function(x) {
    ymax <- max(resT1[[x]],resT2[[x]],resT3[[x]],resT4[[x]],obs1[[x]],obs2[[x]],obs3[[x]],obs4[[x]],na.rm=T)*1.1
    xmax <- max(obs1$day,obs2$day,obs3$day,obs4$day,na.rm=T)+5
    plot(resT1[[x]], type="l", pch=3, xlab="DAS", ylab=x,xaxs="i", yaxs="i", ylim=c(0,ymax), xlim=c(0,xmax), col="black")
    lines(resT2[[x]], col="red")
    lines(resT3[[x]], col="blue")
    lines(resT4[[x]], col="green")
    points(obs1$day, obs1[[x]], type="p", col="black")
    points(obs2$day, obs2[[x]], type="p", col="red")
    points(obs3$day, obs3[[x]], type="p", col="blue")
    points(obs4$day, obs4[[x]], type="p", col="green")

    diff1 <- ((((obs1[[x]] - resP1[[x]])/obs1[[x]])^2))
    diff1 <- sqrt((sum(diff1, na.rm=T))/(sum(!is.na(diff1))))
    diff2 <- ((((obs2[[x]] - resP2[[x]])/obs2[[x]])^2))
    diff2 <- sqrt((sum(diff2, na.rm=T))/(sum(!is.na(diff2))))
    diff3 <- ((((obs3[[x]] - resP3[[x]])/obs3[[x]])^2))
    diff3 <- sqrt((sum(diff3, na.rm=T))/(sum(!is.na(diff3))))
    diff4 <- ((((obs4[[x]] - resP4[[x]])/obs4[[x]])^2))
    diff4 <- sqrt((sum(diff4, na.rm=T))/(sum(!is.na(diff4))))
    return(c(diff1, diff2, diff3, diff4))
  }
  VarListP <- VarList[VarList != "day"]
  sapply(VarListP, plotF)
}
savePar <- function(name = Sys.Date()) {
  resPar <- matrix(as.vector(c(result$value, result$par)), ncol=length(paramOfInterest)+1)
  write.table(resPar, file=paste("par_",name,".csv"), sep=",", append=F, dec=".",col.names=c("RMSE",paramOfInterest),row.names = F)
}
saveAPlot <- function(name = Sys.Date()) {
  pdf(paste(name,"_T.pdf",sep=""))
  resAPlot()
  dev.off()
}

savePlot <- function(obs,meteo,param) {
  pdf(paste("rep.pdf",sep=""))
  resPlot(obs,meteo,param)
  dev.off()
}
saveParF <- function() {
  bestp <- as.vector(result$par)
  paramInitTrans[paramOfInterest] <- bestp
  parameters <- data.frame(Name=ParamList, Values=unlist(paramInitTrans[1,]), row.names=NULL)
  write.table(parameters, "ECOMERISTEM_parameters.txt", sep="=", dec=".", quote=F, row.names=F, col.names=F)
}
coefCompute <- function(length,obsCoef) {
  tmp <- c()
  for(i in 0:(length-1)) {
    tmp <- c(tmp, (obsCoef+(i*(obsCoef/(100/coefIncrease)))))
  }
  return(data.frame(matrix(unlist(tmp),nrow=length,byrow=T)))
}

###OPTIMISATION RUN###
set.seed(1337)
coeff1 <<- coefCompute(nrow(obs1),obsCoef1)
coeff2 <<- coefCompute(nrow(obs2),obsCoef2)
coeff3 <<- coefCompute(nrow(obs3),obsCoef3)
coeff4 <<- coefCompute(nrow(obs4),obsCoef4)
if(clusterA && detectCores() >= 4) {
  nbCores <- detectCores() - 2
  cl <- makeCluster(nbCores, outfile="clusterlog.txt")
  clusterEvalQ(cl, library(recomeristem,rgenoud))
  clusterExport(cl, varlist = c("meteo1","meteo2","meteo3","meteo4","obs1","obs2","obs3","obs4","param1","param2","param3","param4","coeff1","coeff2","coeff3","coeff4","paramOfInterest","isInit"))
} #parallel
time <- system.time(resOptim <- optimisation(Optimizer, maxIter, solTol, bounds))
result <- resOptim[[1]]
result$time <- time
resOptim <- resOptim[[2]]
print(paste("End of Estimation. Elapsed time :", time[[3]],"s"))
stopCluster(cl)
