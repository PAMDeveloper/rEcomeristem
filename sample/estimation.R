#Estimation EcoMeristem

# Gestion des dates
lct <- Sys.getlocale("LC_TIME"); Sys.setlocale("LC_TIME", "C")

library(rgenoud)
library(recomeristem)
library(parallel)
library(pso)
library(ABCoptim)
library(DEoptim)

meteo = recomeristem::getMeteo_from_files("D:/PAMStudio_dev/data/ecomeristem/estimlisa/g1lisa")
paramInit = recomeristem::getParameters_from_files("D:/PAMStudio_dev/data/ecomeristem/estimlisa/g1lisa")
vObs = recomeristem::rcpp_get_vObs("D:/PAMStudio_dev/data/ecomeristem/estimlisa/g1lisa/VobsG1.txt")
ParamList <- as.character(paramInit$Name)
paramInitTrans <- as.data.frame(matrix(ncol=length(ParamList)))
colnames(paramInitTrans) <- ParamList

paramInit[which(paramInit=="nbleaf_pi"),2] = 18
paramInit[which(paramInit=="nb_leaf_stem_elong"),2] = 9
paramInit[which(paramInit=="phenostage_to_end_filling"),2] = 30
paramInit[which(paramInit=="phenostage_to_maturity"),2] = 31

param <- paramInit

paramInitTrans[1,] <- paramInit$Values

resulttmp = recomeristem::rcpp_run_from_dataframe(paramInit, meteo)
obsRed = recomeristem::rcpp_reduceVobs(vObs, resulttmp)
resRed = recomeristem::rcpp_reduceResults(resulttmp, vObs)

ParamOfInterest <- c("Epsib", "Ict", "MGR_init", "plasto_init", "SLAp", "leaf_length_to_IN_length", "coef_MGR_PI", "slope_length_IN", "slope_LL_BL_at_PI", "density_IN1", "density_IN2")

NbParam <- length(ParamOfInterest)
MinValue <- c(6, 1, 5, 20, 10, 0.01, -0.5, 0.0, 0.0, 0.01, 0.1)
MaxValue <- c(9, 10, 14, 60, 120, 0.5, 0.5, 1, 0.5, 0.1, 0.5)
Bounds <- matrix(c(MinValue,MaxValue),ncol=2)
Obs <- length(obsRed)

Optim_Ecomeristem_funct <- function(p){
    paramInitTrans[ParamOfInterest] <- p
    parameters <- data.frame(Name=ParamList, Values=unlist(paramInitTrans[1,]))
    Res_ecomeristem <- recomeristem::rcpp_run_from_dataframe(parameters,meteo)
    res = recomeristem::rcpp_reduceResults(Res_ecomeristem, vObs)
    diff = ((obsRed - res)/obsRed)^2
    #return(sqrt((sum(diff^2, na.rm=T))/sum(colSums(!is.na(diff))))) #RMSE
    return(sqrt(sum(colSums(diff, na.rm=T)/colSums(!is.na(diff)),na.rm=T)))
}

##Parallel
no_cores <- max(1, detectCores() - 2)
cl <- makeCluster(no_cores, outfile="clusterlog.txt")
clusterEvalQ(cl, library(recomeristem,rgenoud))
clusterExport(cl, varlist = c("param","meteo","vObs","obsRed", "paramInitTrans", "ParamOfInterest", "ParamList"))

##DE
resOptimDE3 <- DEoptim(Optim_Ecomeristem_funct, lower=Bounds[,1], upper=Bounds[,2], DEoptim.control(VTR=0,itermax=2000, strategy=2, cluster=cl, packages=c("recomeristem"), parVar=c("param","meteo","vObs","obsRed", "paramInitTrans", "ParamOfInterest", "ParamList") ))

##GENOUD
#resOptimGAFULL <- genoud(Optim_Ecomeristem_funct, NbParam, max=FALSE, pop.size=1000, max.generations=1000, wait.generations=100,
#                    hard.generation.limit=TRUE, MemoryMatrix=TRUE, starting.values=NULL, Domains=Bounds, print.level = 2, boundary.enforcement = 0, gradient.check=FALSE, cluster=cl)

stopCluster(cl)

##PSO
#resOptimPSO <- psoptim(rep(NA,NbParam), Optim_Ecomeristem_funct, lower=Bounds[,1], upper=Bounds[,2], control=list(maxit=10000,abstol=0))

##ABC
#resOptimABC <- abc_optim(Bounds[,1], Optim_Ecomeristem_funct, lb=Bounds[,1], ub=Bounds[,2], maxCycle = 1000, criter=10000, limit = 500)

#resOptimDE2 <- DEoptim(Optim_Ecomeristem_funct, lower=Bounds[,1], upper=Bounds[,2], DEoptim.control(VTR=0,itermax=10000, strategy=2))


##plots result

bestp <- as.vector(resOptimDE3$optim$bestmem)
paramInitTrans[ParamOfInterest] <- bestp
parameters <- data.frame(Name=ParamList, Values=unlist(paramInitTrans[1,]), row.names=NULL)
Res_ecomeristem <- recomeristem::rcpp_run_from_dataframe(parameters,meteo)
res = recomeristem::rcpp_reduceResults(Res_ecomeristem, vObs)
obsRed$day <- vObs$day
par(mfrow=c(2, 2))
plot(Res_ecomeristem$pht, type="l")
points(obsRed$day, obsRed$pht, type="p", col="red")
plot(Res_ecomeristem$lig, type="l")
points(obsRed$day, obsRed$lig, type="p", col="red")
plot(Res_ecomeristem$arealfel, type="l")
points(obsRed$day, obsRed$arealfel, type="p", col="red")
plot(Res_ecomeristem$tillernb_1, type="l")
points(obsRed$day, obsRed$tillernb_1, type="p", col="red")
plot(Res_ecomeristem$biomaero2, type="l")
points(obsRed$day, obsRed$biomaero2, type="p", col="red")
plot(Res_ecomeristem$biomleafmainstem, type="l")
points(obsRed$day, obsRed$biomleafmainstem, type="p", col="red")
plot(Res_ecomeristem$biominmainstem, type="l")
points(obsRed$day, obsRed$biominmainstem, type="p", col="red")
plot(Res_ecomeristem$nbleaf, type="l")
points(obsRed$day, obsRed$nbleaf, type="p", col="red")
plot(Res_ecomeristem$ic, type="l")
plot(Res_ecomeristem$mainstem_stock_in, type="l")
points(obsRed$day, obsRed$mainstem_stock_in, col="red")

###PSOplots
##plot(resOptim4$stats$it, resOptim3$stats$error, log="y", xlab="It", ylab="Erreur", type="l")
#
#DEoptimplots
#plot(resOptim5, type="b")
#plot(resOptim5, plot.type="bestvalit",type="l")

#ABCplots
#par(mfrow=c(2, 1))
#for (i in 1:ncol(resOptim1$hist)) {
#  plot(resOptim1$hist[,i], ylab=i, xlab="It", type="b")
#}


### RGENOUD POP ###
# MODIFIER LE CHEMIN SELON LE RETOUR DE tempdir() DANS READ.TABLE, ne pas oublier genoud.pro à la fin
#tempdir()
#pop <- read.table('C:/Users/U108-J~1/AppData/Local/Temp/RtmpOsHrKz/genoud.pro', comment.char = 'G')
#best <- pop[pop$V1 == 1,, drop = FALSE]
#very.best <- as.matrix(best[nrow(best), 3:ncol(best)])
