#Script sensitivity analysis for EcoMeristem
#Authors : Florian Larue, Gregory Beurier, Jean-Christophe Soulie, Lauriane Rouan
#(PAM, AGAP, BIOS, CIRAD)

##############################Input/Output management##########################
PPath <- "D:/Workspace/PAM/estim/estimlisa"  #path to parameter and valparam files
MPath <- "D:/Workspace/PAM/estim/estimlisa"  #path to weither files
OVariable <- "biomaero2"              #name of observed variable
Method <- "F"                         #SA method (M = morris, C = complete, F = fast)
ValParamName <- "ValParam.txt"        #name of valparam file
###############################################################################


#Install and load packages
list.of.packages <- c("recomeristem","sensitivity","rgl")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only=TRUE)

#Files
parameterFilename <- paste(PPath, "ECOMERISTEM_parameters.txt", sep = "")
ValParamFilename <- paste(PPath, ValParamName, sep = "")
meteo <- recomeristem::getMeteo_from_files(MPath)
parameterFileContents <- recomeristem::getParameters_from_files(PPath)
ValParam <- read.table(ValParamFilename, blank.lines.skip = TRUE, col.names = c("Name", "BInf", "BSup", "Step"))
nbFact <- nrow(ValParam)
ParamAFaireVarier <- as.vector(ValParam[,1])

#Bounds
BornesParamAEstimer <- as.matrix(ValParam[,2:4])
B1 <- BornesParamAEstimer[,1]
B2 <- BornesParamAEstimer[,2]
Step <- BornesParamAEstimer[,3]
NbParam <- length(as.vector(B1))

bounds <- apply(cbind(B1, B2), 1, function(x){list(min = x[1], max = x[2])})

EcoMerist <- function(x) {
  ind <- match(ParamAFaireVarier,parameterFileContents$Name)
  parameterFileContents$Values[ind] <- x
  write.table(parameterFileContents,
            file = parameterFilename,
            quote = FALSE,
            sep = "=",
            row.names = FALSE,
            col.names = FALSE)
  EcoSim <- recomeristem::rcpp_run_from_dataframe(parameterFileContents, meteo)
  ValSimTot <- EcoSim[[OVariable]]
  ValSim <- as.numeric(ValSimTot[length(ValSimTot)])
  return(ValSim)
}

EcoMerist.fun <- function(m) {toto=apply(m, 1, EcoMerist); return(toto)}

##########################Sensitivity Analysis#################################
#Morris
SensitivityEcoMerist <- function() {
  sa<-morris(model=NULL,
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
              paste(PPath, "sortie_plan_complet.txt", sep = ""),
              sep="\t",
              row.names = FALSE,
              col.names = FALSE,
              dec = ".")
}

#Fast
IndicesSensitivityEcoMerist <- function() {
  sa <- fast99(model = NULL,
             factors = nbFact,
             n = 450,
             q = rep("qunif", nbFact),
             q.arg=bounds)
  y <- EcoMerist.fun(sa$X)

  tell(sa, y)

  sortie <- cbind(sa$X, sa$y)
  write.table(sortie,
              paste(PPath, "sortie_plan_fast.txt", sep = ""),
              sep="\t",
              row.names = FALSE,
              col.names = FALSE,
              dec = ".")
  print(sa)
  plot(sa)
}

################################Start Analyse##################################
switch(Method,
       "M" = {
         SensitivityEcoMerist()
       },
       "F" = {
         IndicesSensitivityEcoMerist()
       },
       "C" = {
         PlanComplet()
       }
)
