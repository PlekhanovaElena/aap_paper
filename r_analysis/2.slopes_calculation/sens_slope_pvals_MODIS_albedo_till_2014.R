library(raster)
library(tictoc)
library(trend)
library(zoo)
library(zyp) # fast Sens slope calculation


fun_sens_pval = function(y) {
  if (sum(is.na(y)) >= 10) NA  else 
    tryCatch(zyp.trend.vector(y)[6], 
             error=function(err) zyp.trend.vector(na.aggregate(y))[6])
} 

fun_sens = function(y, x = c(1:15)) {
  if (sum(is.na(y)) >= 10) NA  else zyp.sen(y~x)$coefficients[2]
} 

path_ims = list.files("~/data/modis/snowfree_qa_filtered/sw/", full.names = T)
ims = stack(path_ims[1:15])
print("Start calculations")
tic("sens on 16 nodes")
beginCluster(n = 16)
slps = clusterR(ims, calc, args=list(fun=fun_sens))
endCluster()
toc()

writeRaster(slps,
            "./data/regression_slopes_061/slope_sw2014.tif", 
            overwrite = T)

tic("sens pvalue on 16 nodes")
beginCluster(n = 16)
pvals = clusterR(ims, calc, args=list(fun=fun_sens_pval))
endCluster()
toc()

writeRaster(pvals,
            "./data/regression_slopes_061/pval_sw2014.tif",
            overwrite = T)
