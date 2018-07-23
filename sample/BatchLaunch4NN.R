library(recomeristem)

meteoDataPath <- "D:/Workspace_NN/meteo.txt"
meteoDataDF <- read.csv(file=meteoDataPath, header=FALSE, sep="\t")
param <- recomeristem::getParameters_from_files("D:/Workspace_NN/")
param$Values[param$Name=="BeginDate"] <- 1
param$Values[param$Name=="EndDate"] <- 150
fakeobs <- as.data.frame(matrix(0:149, ncol=2, nrow=150))
names(fakeobs) <- c("day", "nbleaf")
meteo <- data.frame(matrix(vector(), 150, 5,
                           dimnames=list(c(), c("Temperature", "Par", "ETP", "Irrigation", "P"))),
                           stringsAsFactors=F)
meteo[,3:5] <- 0
recomeristem::init_simu(param, meteo, fakeobs, "env1")

res <- as.data.frame(matrix(vector(),0,6))
for(i in seq(from=1, to=nrow(meteoDataDF), by=150)){
  meteo[1:150,1:2] <- meteoDataDF[i:(i+149),c(3,2)]
  r <- recomeristem::launch_simu_meteo("env1", meteo)
  res <- rbind(res,cbind(meteo,r[,2]))
  if (i%%10000 == 1) {
    print(i)
  }
}

write.table(res[c(1,2,6)], file = "D:/Workspace_NN/batch.txt",col.names = F, row.names = F, sep = ";")
