detach("package:recomeristem", unload=TRUE)
library(recomeristem)
library(compare)
path <- "D:/Samples/ecomError"
meteo <- recomeristem::getMeteo_from_files(path)
param <- recomeristem::getParameters_from_files(path)
paramBis <- param
meteoBis <- meteo
result <- recomeristem::rcpp_run_from_dataframe(param, meteo)
res1 <- result$tillernb_1

for (i in seq(1, 20)) {
  meteo2 <- recomeristem::getMeteo_from_files(path)
  param2 <- recomeristem::getParameters_from_files(path)
  result2 <- recomeristem::rcpp_run_from_dataframe(param2, meteo2)
  res2 <- result2$tillernb_1
  #print(compare(res1,res2,TRUE))
  if(!isTRUE(compare(result$tillernb_1,result2$tillernb_1,TRUE))) {
    print(paste("ERREUR at step:",i))
    print(res1)
    print(res2)
    print(compare(meteo,meteo2))
    print(compare(param,param2))
    print(compare(result,result2))
    break
  }
}
