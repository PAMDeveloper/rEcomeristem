#Test EcoMeristem R

library(recomeristem)

# Gestion des dates
lct <- Sys.getlocale("LC_TIME"); Sys.setlocale("LC_TIME", "C")

library(rgenoud)
library(stringr)

setwd("D:/Modeles/Test modele R/ecomeristem")

meteo <<- recomeristem::getMeteo_from_files("sample")
#meteo <<- read.table("meteo.txt",sep='\t', header=TRUE)
paramInit <- recomeristem::getParameters_from_files("sample")
#paramInit <- read.table("ecoparams.txt",sep='\t', header=TRUE)

# ParamList <- as.character(paramInit$Name)
# paramInitTrans <- as.data.frame(matrix(ncol=length(ParamList)))
# colnames(paramInitTrans) <- ParamList
# paramInitTrans[1,] <- paramInit$Values


recomeristem::rcpp_run_from_dataframe(paramInit,meteo)

# #Observations
#
# DataObsInit <- read.table("DataObs.txt",sep='\t', header=TRUE)
#
# DataObsVarietyITK <- subset(DataObsInit, idvariete=="06SB15_12fev13"& iditinerairetechnique=="Sotuba-2011-D1")
#
# ListObs <- list()
#
# VarOfInterest <- c("grainyieldpop")
#
# for (v in VarOfInterest){
#   IdData <- !is.na(DataObsVarietyITK[,v])
#   Obs <- DataObsVarietyITK[IdData,v]
#   DateObs <- as.Date(as.character(DataObsVarietyITK[IdData,"jour"]),format="%d/%m/%Y")
#   ListObs[[v]] <- list(Obs,DateObs)
# }
#
#
#
#
# # Param?tres ? calibrer	Observations	V1	V2	V3
# # LeafLengthMax (mm)	Longueur de la plus grande feuille	1000	1000	900
# # CoefLeafWLRatio	Largeur/Longueur de la +G feuille	0.08	0.08	0.12
# # CoefLeafDeath	Vari?t? STAY GREEN	-	-	0.002
#
# ParamOfInterest <<- c("leaflengthmax","coeffleafwlratio","coeffleafdeath")
# NbParam <- length(ParamOfInterest)
# MinValue <- c(500,0.05,0.001)
# MaxValue <- c(1000,0.2,0.003)
# Bounds <- matrix(c(MinValue,MaxValue),ncol=2)
# Obs <- 5
#
# param <- paramInit
#
# Optim_Samara_funct <- function(p){
#   error_v <- c()
#   paramInitTrans[ParamOfInterest] <- p
#   param <- data.frame(Name=ParamList,Values=unlist(paramInitTrans[1,]))
#   Res_Samara <- rsamara::rcpp_run_from_dataframe(param,meteo)[[1]]
#   colnames(Res_Samara) <- str_to_lower(colnames(Res_Samara))
#   Res_Samara$time <- as.Date(as.character(Res_Samara$time),format="%Y-%b-%d  %H:%M:%S")
#   for (v in VarOfInterest){
#     IdDate <- pmatch(as.character(ListObs[[v]][[2]]),as.character(Res_Samara$time))
#     Simul_V_Date <- Res_Samara[IdDate,v]
#     error_v <- c(error_v,ListObs[[v]][[1]]-Simul_V_Date)
#   }
#   return(sum(error_v^2))
# }
#
# # test <- c(900,0.12,0.002)
# #
# #
# # Optim_Samara_funct(test)
#
#
#
#
# resOptim <- genoud(Optim_Samara_funct, NbParam, max=FALSE, pop.size=1000, max.generations=100, wait.generations=10,
#                    hard.generation.limit=TRUE, starting.values=NULL, MemoryMatrix=TRUE,
#                    Domains=Bounds)
#
#
