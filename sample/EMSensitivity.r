#Script sensitivity analysis for EcoMeristem
#Authors : Florian Larue, Gregory Beurier, Jean-Christophe Soulie, Lauriane Rouan
#(PAM, AGAP, BIOS, CIRAD)

##############################Input/Output management##########################
PPath <- "D:/Workspace/estimworkspace/AS"  #path to parameter and valparam files
MPath <- "D:/Workspace/estimworkspace/AS"  #path to weither files
OVariable <- "stockin"              #name of observed variable
Method <- "F"                         #SA method (M = morris, C = complete, F = fast)
ValParamName <- "ValParam.txt"        #name of valparam file
###############################################################################


#Install and load packages
list.of.packages <- c("recomeristem","sensitivity","rgl")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only=TRUE)

#Files
parameterFilename <- paste(PPath, "ECOMERISTEM_parameters.txt", sep = "/")
ValParamFilename <- paste(PPath, ValParamName, sep = "/")
meteo <- recomeristem::getMeteo_from_files(MPath)
parameterFileContents <- recomeristem::getParameters_from_files(PPath)
ValParam <- read.table(ValParamFilename, blank.lines.skip = TRUE, header = TRUE)
nbFact <- nrow(ValParam)
ParamAFaireVarier <- as.vector(ValParam[,1])

#obs for init simu
obs <- recomeristem::get_clean_obs(paste(PPath,"/vobs_moy.txt",sep=""))
VarList <- names(obs)
if(!OVariable %in% VarList) {
  VarList <- c(VarList,OVariable)
}
VarList
obsT <- data.frame(day = seq(1,nrow(meteo)))
for(v in VarList) {
  if(v != "day") {
    obsT[,v] <- 1
  }
}
recomeristem::init_simu(parameterFileContents, meteo, obsT, "envAS")

#Bounds
BornesParamAEstimer <- as.matrix(ValParam[,2:4])
B1 <- BornesParamAEstimer[,1]
B2 <- BornesParamAEstimer[,2]
Step <- BornesParamAEstimer[,3]
NbParam <- length(as.vector(B1))
bounds <- apply(cbind(B1, B2), 1, function(x){list(min = x[1], max = x[2])})

EcoMerist <- function(x) {
  res <- recomeristem::launch_simu("envAS", ParamAFaireVarier, x)
  ValSimTot <- res[[OVariable]]
  ValSim <- as.numeric(ValSimTot[length(ValSimTot)-1])
  return(ValSim)
}

EcoMerist.fun <- function(m) {toto=apply(m, 1, EcoMerist); return(toto)}

##########################Sensitivity Analysis#################################
#Morris
SensitivityEcoMerist <- function() {
  sa<<-morris(model=NULL,
             factors = ParamAFaireVarier,
             r = 100,
             binf=B1,
             bsup=B2,
             design = list(type = "oat", levels = 10, grid.jump = 5))
  y <- EcoMerist.fun(sa$X)
  tell(sa,y)
  print(sa)
  plot(sa)
  plot3d.morris(sa)
  plot3d(sa)
  return(sa)
}

#Complete
PlanComplet <- function() {
  ListeFacteurs <- vector("list", NbParam)
  for (i in 1:NbParam){
    ListeFacteurs[[i]] <- seq(B1[i], B2[i], Step[i])
  }
  PlanComplet <- expand.grid(ListeFacteurs)
  y <- EcoMerist.fun(PlanComplet)
  sortie <- cbind(PlanComplet, y)
  write.table(sortie,
              paste(PPath, "sortie_plan_complet.txt", sep = "/"),
              sep="\t",
              row.names = FALSE,
              col.names = FALSE,
              dec = ".")
  return(sortie)
}

#Fast
IndicesSensitivityEcoMerist <- function() {
  sa <- fast99(model = NULL,
             factors = nbFact,
             n = 50000,
             q = rep("qunif",nbFact),
             q.arg=bounds)
  y <- EcoMerist.fun(sa$X)

  tell(sa, y)

  sortie <- cbind(sa$X, sa$y)
  write.table(sortie,
              paste(PPath, "sortie_plan_fast.txt", sep = "/"),
              sep="\t",
              row.names = FALSE,
              col.names = FALSE,
              dec = ".")
  print(sa)
  plot(sa)
  return(sa)
}

################################Start Analyse##################################
switch(Method,
       "M" = {
         result <<- SensitivityEcoMerist()
       },
       "F" = {
         result <<- IndicesSensitivityEcoMerist()
       },
       "C" = {
         result <<- PlanComplet()
       }
)
