#Script estimating parameters for EcoMeristem model
#Authors : Florian Larue, Gregory Beurier, Lauriane Rouan, Delphine Luquet
#-- (PAM, AGAP, BIOS, CIRAD)

###SET INFORMATION FOR ESTIMATION###
path <- "D:/Workspace/estimworkspace/Sorghum/2017/phenoarch/G23"
vName <- "vobs_moy.txt"
paramOfInterest <- c("Epsib","Ict","MGR_init","plasto_init","phyllo_init","ligulo_init","density_IN2",
                     "coef_phyllo_PI","coef_ligulo_PI","leaf_length_to_IN_length")
minValue <- c(6, 0.5, 6, 25, 25, 25, 0.01, 1.0, 1.0, 1.0, 0.1)
maxValue <- c(20, 2.5, 14, 45, 45, 45, 0.3, 3.0, 3.0, 3.0, 0.2)
coefIncrease <- 0
maxIter <- 5000
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
setwd(path)

meteo <- recomeristem::getMeteo_from_files(path)
param <- recomeristem::getParameters_from_files(path)

obs <- recomeristem::get_clean_obs(paste(path,vName,sep="/"))
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
    recomeristem::init_simu(param, meteo, obs, "env")
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

  res <- recomeristem::launch_simu("env", paramOfInterest, p)

  diff1 <- ((((obs1 - res[1:nrow(obs1),])/obs1)^2))*coeff1
  diff1 <- sum(sqrt((colSums(diff1, na.rm=T))/(colSums(!is.na(diff1)))),na.rm=T)

  diff2 <- ((((obs2 - res[1:nrow(obs2),])/obs2)^2))*coeff2
  diff2 <- sum(sqrt((colSums(diff2, na.rm=T))/(colSums(!is.na(diff2)))),na.rm=T)

  diff3 <- ((((obs3 - res[1:nrow(obs3),])/obs3)^2))*coeff3
  diff3 <- sum(sqrt((colSums(diff3, na.rm=T))/(colSums(!is.na(diff3)))),na.rm=T)

  diff4 <- ((((obs4 - res[1:nrow(obs4),])/obs4)^2))*coeff4
  diff4 <- sum(sqrt((colSums(diff4, na.rm=T))/(colSums(!is.na(diff4)))),na.rm=T)

  return((diff1+diff2+diff3+diff4)/4)
}
optimisation <- function(Optimizer, maxIter, solTol, bounds) {
  if(clusterA && detectCores() >= 4) {
    resOptim <- DEoptim(optimEcomeristem, lower=bounds[,1], upper=bounds[,2], DEoptim.control(reltol=relTol,steptol=stepTol,VTR=solTol,itermax=maxIter,strategy=2,cluster=cl,packages=c("recomeristem"),parVar=c("meteo","obs","obs1","obs2","obs3","obs4","param","coeff1","coeff2","coeff3","coeff4","paramOfInterest","isInit")))
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
  recomeristem::init_simu(param, meteo, obs, "res")
  resP <<- recomeristem::launch_simu("res", paramOfInterest, p)
  obsT <- data.frame(day = seq(1,max(obs$day,na.rm=T)))
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
      obsT[,v] <- 1
    }
  }
  recomeristem::init_simu(param, meteo, obsT, "total_res")
  resT <<- recomeristem::launch_simu("total_res", paramOfInterest, p)
  plotF <- function(x) {
    ymax <- max(resT[[x]],obs1[[x]],obs2[[x]],obs3[[x]],obs4[[x]],na.rm=T)*1.1
    xmax <- max(obs1$day,obs2$day,obs3$day,obs4$day,na.rm=T)+5
    plot(resT[[x]], type="l", pch=3, xlab="DAS", ylab=x,xaxs="i", yaxs="i", ylim=c(0,ymax), xlim=c(0,xmax), col="black")
    points(obs1$day, obs1[[x]], type="p", col="black")
    points(obs2$day, obs2[[x]], type="p", col="red")
    points(obs3$day, obs3[[x]], type="p", col="blue")
    points(obs4$day, obs4[[x]], type="p", col="green")

    diff1 <- ((((obs1[[x]] - resP[1:nrow(obs1),x])/obs1[[x]])^2))
    diff1 <- sqrt((sum(diff1, na.rm=T))/(sum(!is.na(diff1))))
    diff2 <- ((((obs2[[x]] - resP[1:nrow(obs2),x])/obs2[[x]])^2))
    diff2 <- sqrt((sum(diff2, na.rm=T))/(sum(!is.na(diff2))))
    diff3 <- ((((obs3[[x]] - resP[1:nrow(obs3),x])/obs3[[x]])^2))
    diff3 <- sqrt((sum(diff3, na.rm=T))/(sum(!is.na(diff3))))
    diff4 <- ((((obs4[[x]] - resP[1:nrow(obs4),x])/obs4[[x]])^2))
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
savePlot <- function(name = Sys.Date()) {
  pdf(paste(name,"_T.pdf",sep=""))
  resAPlot()
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
  clusterExport(cl, varlist = c("meteo","obs","obs1","obs2","obs3","obs4","param","coeff1","coeff2","coeff3","coeff4","paramOfInterest","isInit"))
} #parallel
time <- system.time(resOptim <- optimisation(Optimizer, maxIter, solTol, bounds))
result <- resOptim[[1]]
result$time <- time
resOptim <- resOptim[[2]]
print(paste("End of Estimation. Elapsed time :", time[[3]],"s"))
stopCluster(cl)
