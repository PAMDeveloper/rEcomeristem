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
result <- recomeristem::rcpp_run_from_dataframe(paramBis, meteoBis)
res2 <- result$tillernb_1
print(compare(res1,res2))
