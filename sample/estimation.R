#Script estimating parameters for EcoMeristem model
#Authors : Florian Larue, Gregory Beurier, Lauriane Rouan, Delphine Luquet
#-- (PAM, AGAP, BIOS, CIRAD)
###Set informations for parameter estimation###
PPath <- "D:/Workspace/estimlisa_fvobs/2015"
MPath <- "D:/Workspace/estimlisa_fvobs/2015"
VPath <- "D:/Workspace/estimlisa_fvobs/2015"
VName <- "vobs_G1moyINT_C_BFF2015woleaf.txt"
VECName <- "vobs_G1_C_BFF2015_ET_INTwoleaf.txt"
ParamOfInterest <- c("Epsib", "Ict", "MGR_init", "plasto_init", "phyllo_init", "ligulo_init", "leaf_length_to_IN_length", "coef_MGR_PI", "slope_length_IN", "slope_LL_BL_at_PI", "density_IN1", "density_IN2", "coef_plasto_PI", "coef_phyllo_PI", "coef_ligulo_PI", "SLAp")
MinValue <- c(3, 0.5, 5, 20, 20, 20,  0.01, -0.5, 0.0, 0.0, 0.005, 0.1, 1.0, 1.0, 1.0, 20)
MaxValue <- c(8, 2.5, 15, 50, 50, 50 , 0.5, 0.5, 1, 0.5, 0.1, 0.3, 3.5, 3.5, 3.5, 60)
obsCoef <- c(1,1,1,1,1,1,1)
coefIncrease <- 1
Optimizer <- "D" #(D = DE, G = RGenoud)
RmseM <- "RTEST" #(RS = RSME-sum, REC = RMSE-ET, RC = RMSE-coef, RECC = RMSE-ET-coef)
MaxIter <- 5000
Penalty <- 10 #Penalty for simulation outside of SD (RMSE * Penalty)
SolTol <- 0.01 #will be multiplied by the number of observed variables
ACluster <- TRUE  #parallel for machines with at least 4 cores
###End informations for parameter estimation###

#Install and load packages
list.of.packages <- c("Rcpp","recomeristem","parallel","pso", "rgenoud", "ABCoptim", "DEoptim", "truncnorm","data.table")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, library, character.only=TRUE)

#Dates management
lct <- Sys.getlocale("LC_TIME"); Sys.setlocale("LC_TIME", "C")
setwd(VPath)

#Get data for simulation
meteo <- recomeristem::getMeteo_from_files(MPath)
paramInit <- recomeristem::getParameters_from_files(PPath)
vObs <- recomeristem::rcpp_get_vObs(paste(VPath,paste("/",VName, sep=""), sep=""))
obsET <- recomeristem::rcpp_get_vObs(paste(VPath,paste("/",VECName, sep=""), sep=""))
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
res <- list()
coeff <- c()
SolTol <- SolTol * length(VarList)

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
  lastp <<- p
  Res_ecomeristem <- recomeristem::rcpp_run_from_dataframe(parameters,meteo)
  res <- recomeristem::rcpp_reduceResults(Res_ecomeristem, vObs)
  switch(RmseM,
         "RS" = {
           diff <- ((obsRed - res)/obsRed)^2
           return(sum(sqrt((colSums(diff, na.rm=T))/(colSums(!is.na(diff)))),na.rm=T))
         },
         "REC" = {
           diff1 <- abs(1-(data.table::between(res,obsRed-obsETRed,obsRed+obsETRed)))
           diff2 <- diff1 * Penalty
           diff3 <-  replace(diff2, diff2 == 0,1)
           diff <- (((obsRed - res)/obsRed)^2)*diff3
           return(sum(sqrt((colSums(diff, na.rm=T))/(colSums(!is.na(diff)))),na.rm=T))
         },
         "RC" = {
           diff <- (((obsRed - res)/obsRed)^2)*coeff
           rmse = sqrt((colSums(diff, na.rm=T))/(colSums(!is.na(diff))))
           return(sum(rmse*obsCoef, na.rm=T))
         },
         "RECC" = {
           diff1 <- abs(1-(data.table::between(res,obsRed-obsETRed,obsRed+obsETRed)))
           diff2 <- diff1 * Penalty
           diff3 <-  replace(diff2, diff2 == 0,1)
           diff <- ((((obsRed - res)/obsRed)^2)*diff3)*coeff
           return(sum(sqrt((colSums(diff, na.rm=T))/(colSums(!is.na(diff)))),na.rm=T))
         },
         "RTEST" = {
           diff1 <- abs(1-(data.table::between(res,obsRed-(2*obsETRed),obsRed+(2*obsETRed))))
           diff2 <- diff1 * Penalty
           diff3 <- abs(1-(data.table::between(res,obsRed-obsETRed,obsRed+obsETRed)))
           diff4 <-  replace(diff3, diff3 == 0,0.5)
           diff5 <- pmax(diff2,diff4)
           diff <- ((((obsRed - res)/obsRed)^2)*diff5)*coeff
           return(sum(sqrt((colSums(diff, na.rm=T))/(colSums(!is.na(diff)))),na.rm=T))
         }
  )
}
simulatedAnnealing <- function(func, start_par, lower, upper, itermax = 1000, step = 0.1, printlvl = 100) {
  best_par <- current_par <- neighbor_par <- start_par
  best_value <- current_value <- neighbor_value <- func(start_par)

  message("It\tBest\tCurrent")
  message(sprintf("%i\t%.4f\t%.4f", 0L, best_value, current_value))

  for (k in 1:itermax) {
    temp <- (1 - step)^k
    neighbor_par <- rtruncnorm(length(start_par),lower, upper, current_par, 1)
    neighbor_value <- func(neighbor_par)

    if (neighbor_value < current_value || runif(1,0,1) < exp(-(neighbor_value - current_value) / temp)) {
      current_par <- neighbor_par
      current_value <- neighbor_value
    }

    if (neighbor_value < best_value) {
      best_par <- neighbor_par
      best_value <- neighbor_value
    }

    if(k%%printlvl == 0) {
      message(sprintf("%i\t%.4f\t%.4f", k, best_value, current_value))
    }
  }
  return(list(iterations = itermax, value = best_value, par = best_par))
}
optimisation <- function(Optimizer, MaxIter, SolTol, NbParam, Bounds, NbEnv, SDate, EDate) {
  coeff <<- coefCompute()
  switch(Optimizer,
         "D" = {
           if(ACluster && detectCores() >= 4) {
             resOptim <- DEoptim(optimEcomeristem, lower=Bounds[,1], upper=Bounds[,2], DEoptim.control(VTR=SolTol,itermax=MaxIter, strategy=2, cluster=cl, packages=c("recomeristem"), parVar=c("meteo","vObs","obsRed", "paramInitTrans", "ParamOfInterest", "ParamList", "RmseM", "obsETRed","obsCoef","Penalty","coeff")))

           } else {
             resOptim <- DEoptim(optimEcomeristem, lower=Bounds[,1], upper=Bounds[,2], DEoptim.control(VTR=SolTol,itermax=MaxIter, strategy=2))
           }
           res$optimizer <- "Diffential Evolution Optimization"
           res$par <- resOptim$optim$bestmem
           res$value <- resOptim$optim$bestval
           res$iter <- resOptim$optim$iter
         },
         "G" = {
           if(ACluster && detectCores() >= 4) {
             resOptim <- genoud(optimEcomeristem, NbParam, max=FALSE, pop.size=100, max.generations=MaxIter, wait.generations=max(10,MaxIter/10), hard.generation.limit=TRUE, MemoryMatrix=TRUE, starting.values=NULL, Domains=Bounds, print.level=1, boundary.enforcement=2, gradient.check=FALSE, solution.tolerance=SolTol,cluster=cl)
           } else {
             resOptim <- genoud(optimEcomeristem, NbParam, max=FALSE, pop.size=100, max.generations=MaxIter, wait.generations=max(10,MaxIter/10), hard.generation.limit=TRUE, MemoryMatrix=TRUE, starting.values=NULL, Domains=Bounds, print.level=1, boundary.enforcement=2, gradient.check=FALSE, solution.tolerance=SolTol)
           }
           res$optimizer <- "Genetic Optimization Using Derivatives"
           TPath <- gsub("\\","/",tempdir(), fixed=TRUE)
           pop <- read.table(paste(TPath,paste("/","genoud.pro", sep=""),sep=""), comment.char='G')
           best <- pop[which.min(pop$V2),, drop = FALSE]
           res$par <- as.matrix(best[nrow(best), 3:ncol(best)])
           res$value <- as.character(best[2])
           res$iter <- nrow(best)
         },
         "GL" = {
           if(ACluster && detectCores() >= 4) {
             resOptim <- genoud(optimEcomeristem, NbParam, max=FALSE, pop.size=1000, max.generations=MaxIter, wait.generations=max(10,MaxIter/10), hard.generation.limit=TRUE, MemoryMatrix=TRUE, starting.values=NULL, Domains=Bounds, print.level=1, boundary.enforcement=2, gradient.check=FALSE, solution.tolerance=SolTol, lexical=TRUE, cluster=cl)
           } else {
             resOptim <- genoud(optimEcomeristem, NbParam, max=FALSE, pop.size=1000, max.generations=MaxIter, wait.generations=max(10,MaxIter/10), hard.generation.limit=TRUE, MemoryMatrix=TRUE, starting.values=NULL, Domains=Bounds, print.level=1, boundary.enforcement=2, gradient.check=FALSE, solution.tolerance=SolTol, lexical=TRUE)
           }
           res$optimizer <- "Genetic Optimization Using Derivatives"
           TPath <- gsub("\\","/",tempdir(), fixed=TRUE)
           pop <- read.table(paste(TPath,paste("/","genoud.pro", sep=""),sep=""), comment.char='G')
           best <- pop[which.min(pop$V2),, drop = FALSE]
           res$par <- as.matrix(best[nrow(best), 9:ncol(best)])
           res$value <- as.character(best[2:8])
           res$iter <- nrow(best)
         },
         "A" = {
           resOptim <- simulatedAnnealing(optimEcomeristem, start_par=Bounds[,1], lower=Bounds[,1], upper=Bounds[,2], itermax=MaxIter, step=SolTol, printlvl=max(1,MaxIter/20))
           res$optimizer <- "Simulated Annealing"
           res$par <- resOptim$par
           res$value <- resOptim$value
           res$iter <- MaxIter
         }
  )
  print("End of estimation, type fnList() to see function possibilities")
  return(list(res,resOptim))
}
coefCompute <- function() {
  tmp <- c()
  for(i in 0:(obsLength-1)) {
    tmp <- c(tmp, (obsCoef+(i*(obsCoef/(100/coefIncrease)))))
  }
  return(data.frame(matrix(unlist(tmp),nrow=obsLength,byrow=T)))
}
dePlot <- function() {
  if(Optimizer == "D") {
    par(mfrow=c(1,1))
    plot(resOptim, type="b")
    par(mfrow=c(1,1))
    plot(resOptim, plot.type="bestvalit",type="l")
  }
}
allPlot <- function() {
  bestp <- as.vector(res$par)
  paramInitTrans[ParamOfInterest] <- bestp
  parameters <- data.frame(Name=ParamList, Values=unlist(paramInitTrans[1,]), row.names=NULL)
  Res_ecomeristem <- recomeristem::rcpp_run_from_dataframe(parameters,meteo)
  obsRed$day <- vObs$day
  plotF <- function(x) {
    plot(Res_ecomeristem[[x]], type="l", xlab="DAS", ylab=x,ylim=c(1.5*min(Res_ecomeristem[[x]], na.rm=T),1.5*max(Res_ecomeristem[[x]], na.rm=T)))
    return("Done")
  }
  sapply(names(Res_ecomeristem), plotF)
}
saveERes <- function() {
  bestp <- as.vector(res$par)
  paramInitTrans[ParamOfInterest] <- bestp
  parameters <- data.frame(Name=ParamList, Values=unlist(paramInitTrans[1,]), row.names=NULL)
  Res_ecomeristem <- recomeristem::rcpp_run_from_dataframe(parameters,meteo)
  write.csv(Res_ecomeristem, file="resultESimu.csv")
}
saveRes <- function() {
  meteo <- recomeristem::getMeteo_from_files(MPath)
  param <- recomeristem::getParameters_from_files(PPath)
  Res_ecomeristem <- recomeristem::rcpp_run_from_dataframe(paramInit, meteo)
  write.csv(Res_ecomeristem, file="resultSimu.csv")
}
recoverR <- function() {
  if(Optimizer == "G") {
    res$optimizer <<- "Genetic Optimization Using Derivatives"
    TPath <- gsub("\\","/",tempdir(), fixed=TRUE)
    pop <- read.table(paste(TPath,paste("/","genoud.pro", sep=""),sep=""), comment.char='G')
    best <- pop[which.min(pop$V2),, drop = FALSE]
    res$par <<- as.matrix(best[nrow(best), 3:ncol(best)])
    res$value <<- as.character(best[2])
  } else {
    res$optimizer <- "Genetic Optimization Using Derivatives"
    TPath <- gsub("\\","/",tempdir(), fixed=TRUE)
    pop <- read.table(paste(TPath,paste("/","genoud.pro", sep=""),sep=""), comment.char='G')
    best <- pop[which.min(pop$V2),, drop = FALSE]
    res$par <- as.matrix(best[nrow(best), 9:ncol(best)])
    res$value <- as.character(best[2:8])
    res$iter <- nrow(best)
  }
}
runModel <- function() {
  meteo <- recomeristem::getMeteo_from_files(MPath)
  param <- recomeristem::getParameters_from_files(PPath)
  assign("resultModel", recomeristem::rcpp_run_from_dataframe(paramInit, meteo), envir=.GlobalEnv)
  print("Done, check resultModel for simulation results")
}
runEModel <- function(par = res$par) {
  bestp <- as.vector(par)
  paramInitTrans[ParamOfInterest] <- bestp
  parameters <- data.frame(Name=ParamList, Values=unlist(paramInitTrans[1,]), row.names=NULL)
  assign("resultEModel", recomeristem::rcpp_run_from_dataframe(parameters,meteo), envir=.GlobalEnv)
  print("Done, check resultEModel for simulation results")
}
saveEModel <- function(par = res$par) {
  runEModel(par)
  write.table(resultEModel, "resultEModel.csv", sep=",", dec=".", quote=F, row.names=F, col.names=T)
}
resEPlot <- function(par = res$par) {
  bestp <- as.vector(par)
  paramInitTrans[ParamOfInterest] <- bestp
  parameters <- data.frame(Name=ParamList, Values=unlist(paramInitTrans[1,]), row.names=NULL)
  Res_ecomeristem <- recomeristem::rcpp_run_from_dataframe(parameters,meteo)
  resRed <- recomeristem::rcpp_reduceResults(Res_ecomeristem,vObs)
  obsRed$day <- vObs$day
  plotF <- function(x) {
    if(RmseM != "RS" && RmseM != "RX") {
      plot(Res_ecomeristem[[x]], type="l", xlab="DAS", ylab=x,ylim=c(min(min(Res_ecomeristem[[x]], na.rm=T),min(obsRed[[x]]-obsETRed[[x]], na.rm=T)),max(max(Res_ecomeristem[[x]], na.rm=T),max(obsRed[[x]]+obsETRed[[x]], na.rm=T))))
      points(obsRed$day, obsRed[[x]], type="p", col="red")
      if(!is.null(obsETRed[[x]])) {
        arrows(obsRed$day,obsRed[[x]]-obsETRed[[x]],obsRed$day,obsRed[[x]]+obsETRed[[x]], code=3, length=0.02, angle = 90)
      }
      if(!is.null(obsRed[[x]])) {
        diff1 = abs(1-between(resRed[[x]],obsRed[[x]]-obsETRed[[x]],obsRed[[x]]+obsETRed[[x]]))
        diff2 = ((abs(obsRed[[x]] + obsETRed[[x]] - resRed[[x]]))/obsRed[[x]])^2
        diff3 = ((abs(obsRed[[x]] - obsETRed[[x]] - resRed[[x]]))/obsRed[[x]])^2
        return(sum(pmin(diff1*999,diff2,diff3),na.rm=T))
      } else {
        return(NULL)
      }
    } else {
      plot(Res_ecomeristem[[x]], type="l", xlab="DAS", ylab=x,ylim=c(min(Res_ecomeristem[[x]], obsRed[[x]], na.rm=T),max(Res_ecomeristem[[x]], obsRed[[x]], na.rm=T)))
      points(obsRed$day, obsRed[[x]], type="p", col="red")
      return(paste((sum((abs(obsRed[[x]]-resRed[[x]]))/(obsRed[[x]]),na.rm=T))*100,"%", sep=""))
    }

  }
  sapply(VarList, plotF)
}
savePar <- function(name = Sys.time()) {
  resPar <- matrix(as.vector(res$par), ncol=11)
  write.table(resPar, file="par.csv", sep=",", append=T, dec=".",col.names = F,row.names = name)
}
savePlots <- function(name = Sys.Date()) {
  bestp <- as.vector(res$par)
  paramInitTrans[ParamOfInterest] <- bestp
  parameters <- data.frame(Name=ParamList, Values=unlist(paramInitTrans[1,]), row.names=NULL)
  Res_ecomeristem <- recomeristem::rcpp_run_from_dataframe(parameters,meteo)
  resRed <- recomeristem::rcpp_reduceResults(Res_ecomeristem,vObs)
  obsRed$day <- vObs$day
  pdf(paste(name,".pdf",sep=""))
  plotF <- function(x) {
    if(RmseM != "RS" && RmseM != "RX") {
      plot(Res_ecomeristem[[x]], type="l", xlab="DAS", ylab=x,ylim=c(min(min(Res_ecomeristem[[x]], na.rm=T),min(obsRed[[x]]-obsETRed[[x]], na.rm=T)),max(max(Res_ecomeristem[[x]], na.rm=T),max(obsRed[[x]]+obsETRed[[x]], na.rm=T))))
      points(obsRed$day, obsRed[[x]], type="p", col="red")
      if(!is.null(obsETRed[[x]])) {
        arrows(obsRed$day,obsRed[[x]]-obsETRed[[x]],obsRed$day,obsRed[[x]]+obsETRed[[x]], code=3, length=0.02, angle = 90)
      }
      if(!is.null(obsRed[[x]])) {
        return(paste((abs(sum((abs((obsRed[[x]] - resRed[[x]]))-abs(obsETRed[[x]]))/(obsRed[[x]]), na.rm=T))*100),"%", sep=""))
      } else {
        return(NULL)
      }
    } else {
      plot(Res_ecomeristem[[x]], type="l", xlab="DAS", ylab=x,ylim=c(min(Res_ecomeristem[[x]], obsRed[[x]], na.rm=T),max(Res_ecomeristem[[x]], obsRed[[x]], na.rm=T)))
      points(obsRed$day, obsRed[[x]], type="p", col="red")
      return(paste((sum((abs(obsRed[[x]]-resRed[[x]]))/(obsRed[[x]]),na.rm=T))*100,"%", sep=""))
    }
  }
  sapply(VarList, plotF)
  dev.off()
}
saveAPlots <- function(name= Sys.Date()) {
  bestp <- as.vector(res$par)
  paramInitTrans[ParamOfInterest] <- bestp
  parameters <- data.frame(Name=ParamList, Values=unlist(paramInitTrans[1,]), row.names=NULL)
  Res_ecomeristem <- recomeristem::rcpp_run_from_dataframe(parameters,meteo)
  obsRed$day <- vObs$day
  pdf(paste("ALL",paste(name,".pdf",sep=""),sep=""))
  plotF <- function(x) {
    plot(Res_ecomeristem[[x]], type="l", xlab="DAS", ylab=x,ylim=c(1.5*min(Res_ecomeristem[[x]], na.rm=T),1.5*max(Res_ecomeristem[[x]], na.rm=T)))
    return("Done")
  }
  sapply(names(Res_ecomeristem), plotF)
  dev.off()
}
saveParF <- function() {
  bestp <- as.vector(res$par)
  paramInitTrans[ParamOfInterest] <- bestp
  parameters <- data.frame(Name=ParamList, Values=unlist(paramInitTrans[1,]), row.names=NULL)
  write.table(parameters, "ECOMERISTEM_parameters.txt", sep="=", dec=".", quote=F, row.names=F, col.names=F)
}

saveErrParF <- function(parameters) {
  write.table(parameters, "ECOMERISTEM_parameters_err.txt", sep="=", dec=".", quote=F, row.names=F, col.names=F)
}
resPlot <- function() {
  bestp <- as.vector(res$par)
  paramInitTrans[ParamOfInterest] <- bestp
  parameters <- data.frame(Name=ParamList, Values=unlist(paramInitTrans[1,]), row.names=NULL)
  Res_ecomeristem <- recomeristem::rcpp_run_from_dataframe(parameters,meteo)
  resRed <- recomeristem::rcpp_reduceResults(Res_ecomeristem,vObs)
  obsRed$day <- vObs$day
  plotF <- function(x) {
    if(RmseM != "RS" && RmseM != "RX") {
      plot(Res_ecomeristem[[x]], type="l", xlab="DAS", ylab=x,ylim=c(min(min(Res_ecomeristem[[x]], na.rm=T),min(obsRed[[x]]-obsETRed[[x]], na.rm=T)),max(max(Res_ecomeristem[[x]], na.rm=T),max(obsRed[[x]]+obsETRed[[x]], na.rm=T))))
      points(obsRed$day, obsRed[[x]], type="p", col="red")
      if(!is.null(obsETRed[[x]])) {
        arrows(obsRed$day,obsRed[[x]]-obsETRed[[x]],obsRed$day,obsRed[[x]]+obsETRed[[x]], code=3, length=0.02, angle = 90)
      }
    } else {
      plot(Res_ecomeristem[[x]], type="l", xlab="DAS", ylab=x,ylim=c(min(Res_ecomeristem[[x]], obsRed[[x]], na.rm=T),max(Res_ecomeristem[[x]], obsRed[[x]], na.rm=T)))
      points(obsRed$day, obsRed[[x]], type="p", col="red")
      diff <- ((obsRed[[x]] - resRed[[x]])/obsRed[[x]])^2
    }
    diff <- ((obsRed[[x]] - resRed[[x]])/obsRed[[x]])^2
    return(sqrt((sum(diff, na.rm=T))/(sum(!is.na(diff)))))
  }
  sapply(VarList, plotF)
}
fnList <- function() {
  flist <- c("resPlot() : Plot estimated simulation restults and observations",
             "resEPlot(''vector of parameters'') : Plot simulation results and observations with redefined estimated parameter values",
             "saveRes() : Save estimated simulation results in csv file",
             "recoverR() : Recover best estimation after intentional RGENOUD stop",
             "dePlot() : Plot DE convergence and estimation plots",
             "savePar() : Save estimated parameter values in csv file",
             "saveParF() : Save all parameter (fixed + estimated) in ECOMERISTEM_parameters.txt",
             "savePlots() : Save resPlot() in pdf file",
             "runModel() : Run model with parameter and weither files",
             "runEModel() : Run model with estimated parameters (also possible to give parameter values as argument)",
             "allPlot() : Plot all output variables of model (without observations)")
  print(sort(flist))
}

#Optimisation run
if(ACluster && detectCores() >= 4) {
  nbCores <- detectCores() - 2
  cl <- makeCluster(nbCores, outfile="clusterlog.txt")
  clusterEvalQ(cl, library(recomeristem,rgenoud))
  clusterExport(cl, varlist = c("meteo","vObs","obsRed", "paramInitTrans", "ParamOfInterest", "ParamList", "RmseM", "obsETRed", "obsCoef","Penalty","coeff"))
} #parallel
time <- system.time(resOptim <- optimisation(Optimizer, MaxIter, SolTol, NbParam, Bounds, NbEnv, SDate, EDate))
res <- resOptim[[1]]
res$time <- time
resOptim <- resOptim[[2]]
stopCluster(cl)
print(paste("Elapsed time :", paste(time[[3]],"s")))
