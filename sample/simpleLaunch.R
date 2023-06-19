detach("package:recomeristem", unload=TRUE)
library(ecomeristem)
require(Rcpp)
require(stats)
library(compare)

### TEST perf and launch
path <- "D:/Samples/_Estimation/G1"
VName <- "vobs_moy.txt"
meteo <- recomeristem::getMeteo_from_files(path)
param <- recomeristem::getParameters_from_files(path)
obs <- recomeristem::get_clean_obs(paste(path,paste("/",VName, sep=""), sep=""))
recomeristem::init_simu(param, meteo, obs, "test")
launch <- function(){
  for (i in seq(1, 1000)) {
    result <- recomeristem::launch_simu("test")
  }
}
system.time(launch())


### TEST repro
resRef <- recomeristem::launch_simu("test")
meteo2 <- recomeristem::getMeteo_from_files(path)
param2 <- recomeristem::getParameters_from_files(path)
recomeristem::init_simu(param2, meteo2, obs, "test2")
for (i in seq(1, 200)) {
  res2 <- recomeristem::launch_simu("test2")
  if(!isTRUE(compare(resRef, res2 ,TRUE))) {
    print(paste("ERREUR at step:",i))
    break
  } else {
    print (paste("OK at step:",i))
  }
}

### DEOptim single proc
require(DEoptim)
names <- c("Epsib", "Ict", "MGR_init", "plasto_init", "phyllo_init",
                     "ligulo_init", "coef_MGR_PI", "slope_length_IN", "density_IN2", "coef_plasto_PI",
                     "coef_phyllo_PI", "coef_ligulo_PI", "slope_LL_BL_at_PI",
                     "thresAssim","thresINER","thresLER","thresLEN","stressBP","stressBP2")
MinValue <- c(3, 0.5, 6, 20, 20, 20, -0.5, 0.5, 0.08, 1, 1, 1, 0.0,1,1,1,1,0,1)
MaxValue <- c(8, 2.5, 14, 45, 45, 45, 0.5, 1, 0.3, 3.0, 3.0, 3.0, 0.4,20,20,20,20,10,10)

Bounds <- matrix(c(MinValue,MaxValue),ncol=2)
res <- recomeristem::launch_simu("test")
optim <- function(p){
  res <- recomeristem::launch_simu("test", names, p)
  diff <- ((obs - res)/obs)^2
  return(sum(sqrt((colSums(diff, na.rm=T))/(colSums(!is.na(diff)))),na.rm=T))
}

system.time(resOptim <- DEoptim(optim, lower=Bounds[,1], upper=Bounds[,2], DEoptim.control(VTR=0.1,itermax=1000, strategy=2)))
