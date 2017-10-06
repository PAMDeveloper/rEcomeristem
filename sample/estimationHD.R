#Script estimating parameters for EcoMeristem model
#Authors : Florian Larue, Gregory Beurier, Lauriane Rouan, Delphine Luquet
#-- (PAM, AGAP, BIOS, CIRAD)

###############################################################################
Path <- "D:/Workspace/R/estim/estimlisa"
VName <- "vobs.txt"
VETName <- "vobsET.txt"
ParamOfInterest <- c("Epsib", "Ict", "MGR_init", "plasto_init", "SLAp", "leaf_length_to_IN_length", "coef_MGR_PI", "slope_length_IN", "slope_LL_BL_at_PI", "density_IN1", "density_IN2")
MinValue <- c(1, 1, 1, 20, 10, 0.01, -0.5, 0.0, 0.0, 0.01, 0.1)
MaxValue <- c(10, 10, 20, 60, 120, 0.5, 0.5, 2, 0.5, 0.1, 0.5)
obsCoef <- c(1,1,1,1,1,1,1,1)
Optimizer <- "D" #(D = DE, G = RGenoud, A = Simulated Annealing)
RmseM <- "RET" #(RS = RSME-sum, RET = RMSE-ET, RC = RMSE-coef, RETC = RMSE-ET-coef)
MaxIter <- 1
SolTol <- 0.05
ACluster <- TRUE  #Active la parallelisation pour les machines a minimum 4 coeurs
genPar <- c("density")
###############################################################################
#Install and load packages
list.of.packages <- c("Rcpp","recomeristem","parallel","pso", "rgenoud", "ABCoptim", "DEoptim", "truncnorm","data.table")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, library, character.only=TRUE)
lct <- Sys.getlocale("LC_TIME"); Sys.setlocale("LC_TIME", "C")
###############################################################################
#functions
Dat2Jul <- function( yr=-4712, mo=1, dy=1, hr=12 ) { #
  if ( mo <= 2 )  { yr <- yr-1; mo <-  mo+12 } # Jan, Feb
  H <-  hr/24
  if ( (yr*100 + mo )*100 + dy >= 15821015 ) {  # i.e. 1582-10-15
    # Gregorian calendar:
    A <- yr %/% 100  # centuries
    B <- 2 - A + A %/% 4 # correction for yr=k*400
  } else if (  (yr*100 + mo )*100 + dy <= 15821004 )  {
    # Julian calendar:
    #    if (yr < 0) yr <- yr + 1  # no yr == 0
    B <- 0
  } else  {  # "no legal date, hiatus of calendars
    JD <- NA
    return(JD)
  }
  JD <- floor(365.25*(yr+4716))  # days of julian centuries since 4715 BC
  JD <- JD + floor(30.6001*(mo+1))  # days since March 1, except for Jan, Feb
  JD <- JD + dy + H + B
  JD <- JD - 1524.5
  return( JD )
}
Optim_Ecomeristem_funct <- function(p){
  paramInitTrans[ParamOfInterest] <- p
  parameters <- data.frame(Name=ParamList, Values=unlist(paramInitTrans[1,]))
  Res_ecomeristem <- recomeristem::rcpp_run_from_dataframe(parameters,meteo)
  res <- recomeristem::rcpp_reduceResults(Res_ecomeristem, vObs)
  switch(RmseM,
         "RS" = {
           diff <- ((obsRed - res)/obsRed)^2
           return(sum(sqrt((colSums(diff, na.rm=T))/(colSums(!is.na(diff)))),na.rm=T))
         },
         "RET" = {
           diff1 <- abs(1-(data.table::between(res,obsRed-obsETRed,obsRed+obsETRed)))
           diff2 <- (abs(obsRed+obsETRed-res)/obsRed)^2
           diff3 <- (abs(obsRed-obsETRed-res)/obsRed)^2
           diff <- pmin(as.matrix(diff1*999),as.matrix(diff2),as.matrix(diff3))
           return(sum(sqrt((colSums(diff, na.rm=T))/(colSums(!is.na(diff)))),na.rm=T))
         },
         "RC" = {
           diff <- ((obsRed - res)/obsRed)^2
           rmse = sqrt((colSums(diff, na.rm=T))/(colSums(!is.na(diff))))
           return(sum(rmse*obsCoef, na.rm=T))
         },
         "RETC" = {
           diff1 <- abs(1-(data.table::between(res,obsRed-obsETRed,obsRed+obsETRed)))
           diff2 <- (abs(obsRed+obsETRed-res)/obsRed)^2
           diff3 <- (abs(obsRed-obsETRed-res)/obsRed)^2
           diff <- pmin(as.matrix(diff1*999),as.matrix(diff2),as.matrix(diff3))
           rmse = sqrt((colSums(diff, na.rm=T))/(colSums(!is.na(diff))))
           return(sum(rmse*obsCoef, na.rm=T))
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
optimisation <- function(Optimizer, MaxIter, SolTol, NbParam, Bounds, vobs) {
  writeLines(paste("\nEstimation de",genList[[i]],sep=" "))
  paramInit[paramInit$Name=="BeginDate","Values"] <<- getBdate(i)
  paramInit[paramInit$Name=="EndDate","Values"] <<- getEdate(i)
  lapply(genPar, function(x) paramUpdate(x,i))
  print(paramInit[paramInit$Name == "density", "Values"])
  assign("vObs", vobs, envir = .GlobalEnv)
  assign("obsET", vobsETList[[i]], envir = .GlobalEnv)
  assign("resulttmp", recomeristem::rcpp_run_from_dataframe(paramInit,meteo), envir = .GlobalEnv)
  assign("obsRed", recomeristem::rcpp_reduceVobs(vObs,resulttmp), envir = .GlobalEnv)
  assign("obsETRed", recomeristem::rcpp_reduceVobs(obsET,resulttmp), envir = .GlobalEnv)
  switch(Optimizer,
         "D" = {
           if(ACluster && detectCores() >= 4) {
             resOptim[[genList[[i]]]] <<- DEoptim(Optim_Ecomeristem_funct, lower=Bounds[,1], upper=Bounds[,2], DEoptim.control(VTR=SolTol,itermax=MaxIter, strategy=2, cluster=cl, packages=c("recomeristem"), parVar=c("meteo","vObs","obsRed", "paramInitTrans", "ParamOfInterest", "ParamList", "RmseM", "obsETRed","obsCoef")))

           } else {
             resOptim[[genList[[i]]]] <<- DEoptim(Optim_Ecomeristem_funct, lower=Bounds[,1], upper=Bounds[,2], DEoptim.control(VTR=SolTol,itermax=MaxIter, strategy=2))
           }
           resT$optimizer <- "Diffential Evolution Optimization"
           resT$par <- resOptim[[genList[[i]]]]$optim$bestmem
           resT$value <- resOptim[[genList[[i]]]]$optim$bestval
           resT$iter <- resOptim[[genList[[i]]]]$optim$iter
         },
         "G" = {
           if(ACluster && detectCores() >= 4) {
             resOptim[[genList[[i]]]] <<- genoud(Optim_Ecomeristem_funct, NbParam, max=FALSE, pop.size=100, max.generations=MaxIter, wait.generations=max(10,MaxIter/10), hard.generation.limit=TRUE, MemoryMatrix=TRUE, starting.values=NULL, Domains=Bounds, print.level=1, boundary.enforcement=2, gradient.check=FALSE, solution.tolerance=SolTol,cluster=cl)
           } else {
             resOptim[[genList[[i]]]] <<- genoud(Optim_Ecomeristem_funct, NbParam, max=FALSE, pop.size=100, max.generations=MaxIter, wait.generations=max(10,MaxIter/10), hard.generation.limit=TRUE, MemoryMatrix=TRUE, starting.values=NULL, Domains=Bounds, print.level=1, boundary.enforcement=2, gradient.check=FALSE, solution.tolerance=SolTol)
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
             resOptim[[genList[[i]]]] <<- genoud(Optim_Ecomeristem_funct, NbParam, max=FALSE, pop.size=1000, max.generations=MaxIter, wait.generations=max(10,MaxIter/10), hard.generation.limit=TRUE, MemoryMatrix=TRUE, starting.values=NULL, Domains=Bounds, print.level=1, boundary.enforcement=2, gradient.check=FALSE, solution.tolerance=SolTol, lexical=TRUE, cluster=cl)
           } else {
             resOptim[[genList[[i]]]] <<- genoud(Optim_Ecomeristem_funct, NbParam, max=FALSE, pop.size=1000, max.generations=MaxIter, wait.generations=max(10,MaxIter/10), hard.generation.limit=TRUE, MemoryMatrix=TRUE, starting.values=NULL, Domains=Bounds, print.level=1, boundary.enforcement=2, gradient.check=FALSE, solution.tolerance=SolTol, lexical=TRUE)
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
           resOptim[[genList[[i]]]] <<- simulatedAnnealing(Optim_Ecomeristem_funct, start_par=Bounds[,1], lower=Bounds[,1], upper=Bounds[,2], itermax=MaxIter, step=SolTol, printlvl=max(1,MaxIter/20))
           res$optimizer <- "Simulated Annealing"
           res$par <- resOptim[[genList[[i]]]]$par
           res$value <- resOptim[[x]]$value
           res$iter <- MaxIter
         }
  )
  writeLines(paste("\nFin estimation de",genList[[i]],sep=" "))
  assign("i", i + 1, envir = .GlobalEnv)
  return(resT)
}
paramUpdate <- function(x,i) {
  paramInit[paramInit$Name==x,"Values"] <<- getPValue(x,i)
}
getBdate <- function(i) {
  tmp <- read.csv("vobs.txt", header=T, sep="")
  return(Dat2Jul(as.numeric(as.list(strsplit(as.character(tmp[tmp$ligne==i,]$StartingDate[[1]]), split="/")[[1]])[[3]]),as.numeric(as.list(strsplit(as.character(tmp[tmp$ligne==i,]$StartingDate[[1]]), split="/")[[1]])[[2]]),as.numeric(as.list(strsplit(as.character(tmp[tmp$ligne==i,]$StartingDate[[1]]), split="/")[[1]])[[1]])))
}
getEdate <- function(i) {
  tmp <- read.csv("vobs.txt", header=T, sep="")
  return(Dat2Jul(as.numeric(as.list(strsplit(as.character(tmp[tmp$ligne==i,]$EndingDate[[1]]), split="/")[[1]])[[3]]),as.numeric(as.list(strsplit(as.character(tmp[tmp$ligne==i,]$EndingDate[[1]]), split="/")[[1]])[[2]]),as.numeric(as.list(strsplit(as.character(tmp[tmp$ligne==i,]$EndingDate[[1]]), split="/")[[1]])[[1]])))
}
getPValue <- function(x,i) {
  return(vobsGlobal[vobsGlobal$ligne==i,x][1])
}
#result functions
resE <- function(i) {return(res[[i]])}
dePlot <- function(i) {
  if(Optimizer == "D") {
    par(mfrow=c(1,1))
    plot(resOptim[[i]], type="b")
    par(mfrow=c(1,1))
    plot(resOptim[[i]], plot.type="bestvalit",type="l")
  }
}
resPlot <- function(i) {
  bestp <- as.vector(res[[i]]$par)
  paramInitTrans[ParamOfInterest] <- bestp
  parameters <- data.frame(Name=ParamList, Values=unlist(paramInitTrans[1,]), row.names=NULL)
  vobsI <- as.data.frame(vobsList[[i]])
  Res_ecomeristem <- recomeristem::rcpp_run_from_dataframe(parameters,meteo)
  obsRed <- recomeristem::rcpp_reduceVobs(vobsI,Res_ecomeristem)
  resRed <- recomeristem::rcpp_reduceResults(Res_ecomeristem,vobsI)
  obsETRed <- recomeristem::rcpp_reduceVobs(vobsETList[[i]], Res_ecomeristem)
  VarList <- names(obsRed)
  obsRed$day <- vobsI$day
  plotF <- function(x) {
    if(RmseM != "RS" && RmseM != "RC") {
      plot(Res_ecomeristem[[x]], type="l", xlab="DAS", ylab=x,ylim=c(min(min(Res_ecomeristem[[x]], na.rm=T),min(obsRed[[x]]-obsETRed[[x]], na.rm=T)),max(max(Res_ecomeristem[[x]], na.rm=T),max(obsRed[[x]]+obsETRed[[x]], na.rm=T))), main=paste("geno",i,sep="_"))
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
}
savePlot <- function(i) {
  pdf(paste(paste("geno",i,sep="_"),".pdf",sep=""))
  resPlot(i)
  dev.off()
}
savePar <- function(i) {
  resPar <- matrix(as.vector(res[[i]]$par), ncol=11)
  write.table(resPar, file="par.csv", sep=",", append=T, dec=".", col.names = F, row.names = paste("geno",i,sep="_"))
}
simPlot <- function(i) {
  bestp <- as.vector(res[[i]]$par)
  paramInitTrans[ParamOfInterest] <- bestp
  parameters <- data.frame(Name=ParamList, Values=unlist(paramInitTrans[1,]), row.names=NULL)
  Res_ecomeristem <- recomeristem::rcpp_run_from_dataframe(parameters,meteo)
  plotF <- function(x) {
    plot(Res_ecomeristem[[x]], type="l", xlab="DAS", ylab=x,ylim=c(1.5*min(Res_ecomeristem[[x]], na.rm=T),1.5*max(Res_ecomeristem[[x]], na.rm=T)), main=paste("geno",i,sep="_"))
    return("Done")
  }
  sapply(names(Res_ecomeristem), plotF)
}
saveRes <- function(i) {
  bestp <- as.vector(res[[i]]$par)
  paramInitTrans[ParamOfInterest] <- bestp
  parameters <- data.frame(Name=ParamList, Values=unlist(paramInitTrans[1,]), row.names=NULL)
  Res_ecomeristem <- recomeristem::rcpp_run_from_dataframe(parameters,meteo)
  write.csv(Res_ecomeristem, file=paste(paste("geno",i,sep="_"),"_res.csv",sep=""))
}
resAPlot <- function() {
  lapply(1:length(lineList), function(x) resPlot(x))
}
deAPlot <- function() {
  lapply(1:length(lineList), function(x) dePlot(x))
}
saveAPlot <- function() {
  lapply(1:length(lineList), function(x) savePlot(x))
}
saveAPar <- function() {
  lapply(1:length(lineList), function(x) savePar(x))
}
simAPlot <- function() {
  lapply(1:length(lineList), function(x) simPlot(x))
}
saveARes <- function() {
  lapply(1:length(lineList), function(x) saveRes(x))
}

###############################################################################

###############################################################################
setwd(Path)
NbParam <- length(ParamOfInterest)
Bounds <- matrix(c(MinValue,MaxValue),ncol=2)
meteo <- recomeristem::getMeteo_from_files(Path)
paramInit <- recomeristem::getParameters_from_files(Path)
ParamList <- as.character(paramInit$Name)
paramInitTrans <- as.data.frame(matrix(ncol=length(ParamList)))
colnames(paramInitTrans) <- ParamList
paramInitTrans[1,] <- paramInit$Values
vobsGlobal <- recomeristem::rcpp_get_vObs(paste(Path,paste("/",VName, sep=""), sep=""))
lineList <- as.list(unique(vobsGlobal$ligne))
genList <- as.list(paste("geno",lineList,sep="_"))
vobsList <- lapply(lineList, function(x) vobsGlobal[vobsGlobal$ligne == x,])
vobsETGlobal <- recomeristem::rcpp_get_vObs(paste(Path,paste("/",VETName, sep=""), sep=""))
vobsETList <- lapply(lineList, function(x) vobsETGlobal[vobsETGlobal$ligne == x,])
i <- 1
resOptim <- list()
resT <- list()
###############################################################################
#First model run
resulttmp <- recomeristem::rcpp_run_from_dataframe(paramInit, meteo)
obsRed <- recomeristem::rcpp_reduceVobs(vobsList[[1]], resulttmp)
obsETRed <- recomeristem::rcpp_reduceVobs(vobsETList[[1]], resulttmp)
resRed <- recomeristem::rcpp_reduceResults(resulttmp, vobsList[[1]])
VarList <- names(obsRed)
#Paraller
if(ACluster && detectCores() >= 4) {
  nbCores <- detectCores() - 2
  cl <- makeCluster(nbCores, outfile="clusterlog.txt")
  clusterEvalQ(cl, library(recomeristem,rgenoud))
  clusterExport(cl, varlist = c("meteo", "paramInitTrans", "ParamOfInterest", "ParamList", "RmseM","obsETRed", "obsCoef"))
}
###############################################################################
#Estimation run
res <- lapply(vobsList, function(x) optimisation(Optimizer, MaxIter, SolTol, NbParam, Bounds, x))
stopCluster(cl)
writeLines("\nFin d'estimation des differents geno, tapez resE(i) pour les resultats d'estimation du geno i (remplacer i par le numero du geno)")
###############################################################################
