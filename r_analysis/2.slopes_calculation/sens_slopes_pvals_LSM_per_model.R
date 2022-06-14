library(raster)
library(tictoc)
library(trend)
library(zoo)
library(zyp) # fast Sens slope calculation
library(snow)


fun_sens_pval = function(y) {
  if (sum(is.na(y)) >= 10) NA  else 
    tryCatch(zyp.trend.vector(y)[6], 
             error=function(err) zyp.trend.vector(na.aggregate(y))[6])
} 

fun_sens = function(y, x = c(1:15)) {
  if (sum(is.na(y)) >= 10) NA  else zyp.sen(y~x)$coefficients[2]
} 

modnames = list.dirs("/home/cluster/eplekh/scratch/climate_data/15Jul/0.1/tiffs/",
                     full.names = F)[-1]


firstmod = raster("/home/cluster/eplekh/scratch/climate_data/15Jul/0.1/result_tiffs/
  mean_climate_models.tif")


for (modname in modnames) {
  print(modname)
  ims = stack(list.files(
    paste0("/home/cluster/eplekh/scratch/climate_data/15Jul/0.1/albedo/",
           modname), full.names = T))

  tic()
  ims = resample(ims, firstmod, method = "ngb")
  toc()
  
  tic("sens on 16 nodes")
  beginCluster(n = 16)
  slps = clusterR(ims, calc, args=list(fun=fun_sens))
  endCluster()
  toc()
  
  writeRaster(slps,paste0("./scratch/climate_data/0.1/result_tiffs/slopes_per_model/", 
                          "slope_", modname, ".tif"), overwrite = T)
  
  print("Calculated pvals, starting pvalue")
  tic("sens pvalue on 16 nodes")
  beginCluster(n = 16)
  pvals = clusterR(ims, calc, args=list(fun=fun_sens_pval))
  endCluster()
  toc()
  
  writeRaster(pvals,paste0("./scratch/climate_data/0.1/result_tiffs/slopes_per_model/", 
                           "pval_", modname, ".tif"), overwrite = T)
}


