library(raster)
library(tictoc)
library(trend)
library(zoo)
library(zyp) # fast Sens slope calculation


fun_sens_pval = function(y) {
  if (sum(is.na(y)) >= 18) NA  else 
    tryCatch(zyp.trend.vector(y)[6], 
             error=function(err) zyp.trend.vector(na.aggregate(y))[6])
} 

fun_sens = function(y, x = c(1:22)) {
  if (sum(is.na(y)) >= 18) NA  else zyp.sen(y~x)$coefficients[2]
} 

print("Visible")

ims = stack(list.files(
  "/home/cluster/eplekh/data/modis/snowfree_qa_filtered/vis/", 
  full.names = T))

print("Start calculations")
tic("sens on 16 nodes")
beginCluster(n = 16)
slps = clusterR(ims, calc, args=list(fun=fun_sens))
endCluster()
toc()

writeRaster(slps,
            "~/data/regression_slopes_061/sens_slope_vis.tif", 
            overwrite = T)

tic("sens pvalue on 16 nodes")

beginCluster(n = 16)
pvals = clusterR(ims, calc, args=list(fun=fun_sens_pval))
endCluster()
toc()

writeRaster(pvals,
            "~/data/regression_slopes_061/sens_pvals_vis.tif",
            overwrite = T)

print("Visible")

ims = stack(list.files(
  "/home/cluster/eplekh/data/modis/snowfree_qa_filtered/nir/", 
  full.names = T))

print("Start calculations")
tic("sens on 16 nodes")
beginCluster(n = 16)
slps = clusterR(ims, calc, args=list(fun=fun_sens))
endCluster()
toc()

writeRaster(slps,
            "~/data/regression_slopes_061/sens_slope_nir.tif", 
            overwrite = T)

tic("sens pvalue on 16 nodes")

beginCluster(n = 16)
pvals = clusterR(ims, calc, args=list(fun=fun_sens_pval))
endCluster()
toc()

writeRaster(pvals,
            "~/data/regression_slopes_061/sens_pvals_nir.tif",
            overwrite = T)

print("Visible")

ims = stack(list.files(
  "/home/cluster/eplekh/data/modis/snowfree_qa_filtered/sw/", 
  full.names = T))
print(ims)

print("Start calculations")
tic("sens on 16 nodes")
beginCluster(n = 16)
slps = clusterR(ims, calc, args=list(fun=fun_sens))
endCluster()
toc()

writeRaster(slps,
            "~/data/regression_slopes_061/sens_slope_sw.tif", 
            overwrite = T)

tic("sens pvalue on 16 nodes")

beginCluster(n = 16)
pvals = clusterR(ims, calc, args=list(fun=fun_sens_pval))
endCluster()
toc()

writeRaster(pvals,
            "~/data/regression_slopes_061/sens_pvals_sw.tif",
            overwrite = T)
