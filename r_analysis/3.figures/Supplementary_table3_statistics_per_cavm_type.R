library(raster)
library(tictoc)
library(dplyr)
library(rgdal)
library(maptools)

stats_dir = "./data/regression_slopes_061/" # input and output folder
cavm = raster("./data/cavm/cavm_raster_resampled3.tif") ### Raster CAVM
types = na.omit(unique(values(cavm)))
types = types[types < 90]
tr = read.csv("./data/cavm/tr_codes_to_veg_types.csv")


### Test

outp_stats = function(slps, pvals) {
  sslps_range = slps
  sslps_range[pvals > 0.05] = NA
  sign_prop = sum(!is.na(values(sslps_range)))/sum(!is.na(values(slps)))
  vals = na.omit(values(sslps_range))
  pasit = sum(vals > 0)/length(vals)
  cat("There are", round(sign_prop*100), "% of significant slopes", 
      "(5% occuring by chance).")
  cat(" From them,", round(pasit*100), "% are positive.")
}

slps = raster(paste0(stats_dir,"sens_slope_sw.tif"))
pvals = raster(paste0(stats_dir,"sens_pvals_sw.tif"))
outp_stats(slps, pvals)


##########################  Shortwave

slps = raster(paste0(stats_dir,"sens_slope_sw.tif"))
pvals = raster(paste0(stats_dir,"sens_pvals_sw.tif"))

mvals = mean(
  stack(list.files("./data/modis/snowfree_shortwave/", full.names = T)), 
  na.rm = T)
mvals = mvals/1000

sslps_range = slps
sslps_range[pvals > 0.05] = NA
list_vtype = lapply(types, function(type) {
  print(type)
  type_slps = slps
  sslps_range_type = sslps_range
  indx = cavm!=type
  mvals[indx] = NA
  type_slps[indx] = NA
  sslps_range_type[indx] = NA
  
  alb = na.omit(values(mvals))
  vals = na.omit(values(sslps_range_type))
  sign_prop = length(vals)/sum(!is.na(values(type_slps)))
  posit = sum(vals > 0)/length(vals)
  
  return(data.frame(code = type, mean_val = mean(alb), std_val = sd(alb),
                    mean_slp = mean(vals), std_slp = sd(vals), 
                    sign_prop = sign_prop, posit = posit))
})

stats_vtypes = do.call(rbind, list_vtype)
stats_vtypes$vtype = tr$vtype[match(stats_vtypes$code, tr$code)]
outp = stats_vtypes[order(stats_vtypes$vtype),]
outp$mean_slp = outp$mean_slp/1000
outp$std_slp = outp$std_slp/1000
outp$sign_prop = round(outp$sign_prop, 2)*100
outp$posit = round(outp$posit, 2)*100
outp$mean_change = round(outp$mean_slp*22, 4)
write.csv(outp, paste0(stats_dir,"stats_shortwave.csv"), 
          row.names = F)



##########################  VIS

slps = raster(paste0(stats_dir,"sens_slope_vis.tif"))
pvals = raster(paste0(stats_dir,"sens_pvals_vis.tif"))
mvals = mean(
  stack(list.files("./data/modis/snowfree_vis/", full.names = T)), 
  na.rm = T)
mvals = mvals/1000

sslps_range = slps
sslps_range[pvals > 0.05] = NA
list_vtype = lapply(types, function(type) {
  print(type)
  type_slps = slps
  sslps_range_type = sslps_range
  indx = cavm!=type
  mvals[indx] = NA
  type_slps[indx] = NA
  sslps_range_type[indx] = NA
  
  alb = na.omit(values(mvals))
  vals = na.omit(values(sslps_range_type))
  sign_prop = length(vals)/sum(!is.na(values(type_slps)))
  posit = sum(vals > 0)/length(vals)
  
  return(data.frame(code = type, mean_val = mean(alb), std_val = sd(alb),
                    mean_slp = mean(vals), std_slp = sd(vals), 
                    sign_prop = sign_prop, posit = posit))
})

stats_vtypes = do.call(rbind, list_vtype)
stats_vtypes$vtype = tr$vtype[match(stats_vtypes$code, tr$code)]
outp = stats_vtypes[order(stats_vtypes$vtype),]
outp$mean_slp = outp$mean_slp/1000
outp$std_slp = outp$std_slp/1000
outp$sign_prop = round(outp$sign_prop, 2)*100
outp$posit = round(outp$posit, 2)*100
outp$mean_change = round(outp$mean_slp*22, 4)
write.csv(outp,paste0(stats_dir,"stats_vis.csv"), 
          row.names = F)



##########################  NIR

slps = raster(paste0(stats_dir,"sens_slope_nir.tif"))
pvals = raster(paste0(stats_dir,"sens_pvals_nir.tif"))
mvals = mean(
  stack(list.files("./data/modis/snowfree_nir/", full.names = T)), 
  na.rm = T)
mvals = mvals/1000

sslps_range = slps
sslps_range[pvals > 0.05] = NA
list_vtype = lapply(types, function(type) {
  print(type)
  type_slps = slps
  sslps_range_type = sslps_range
  indx = cavm!=type
  mvals[indx] = NA
  type_slps[indx] = NA
  sslps_range_type[indx] = NA
  
  alb = na.omit(values(mvals))
  vals = na.omit(values(sslps_range_type))
  sign_prop = length(vals)/sum(!is.na(values(type_slps)))
  posit = sum(vals > 0)/length(vals)
  
  return(data.frame(code = type, mean_val = mean(alb), std_val = sd(alb),
                    mean_slp = mean(vals), std_slp = sd(vals), 
                    sign_prop = sign_prop, posit = posit))
})

stats_vtypes = do.call(rbind, list_vtype)
stats_vtypes$vtype = tr$vtype[match(stats_vtypes$code, tr$code)]
outp = stats_vtypes[order(stats_vtypes$vtype),]
outp$mean_slp = outp$mean_slp/1000
outp$std_slp = outp$std_slp/1000
outp$sign_prop = round(outp$sign_prop, 2)*100
outp$posit = round(outp$posit, 2)*100
outp$mean_change = round(outp$mean_slp*22, 4)
write.csv(outp,paste0(stats_dir,"stats_nir.csv"), 
          row.names = F)


############################################################



### Shape CAVM

shcavm = readOGR("./data/cavm/shape/cp_veg_la.shp") 

cs1 = subset(shcavm, PHYSIOG==1)
slp_s1 = extract(slps, cs1) #Note: it takes 1 min
slp_s1 = unlist(slp_s1)
s1 = extract(sslps_range, cs1) #Note: it takes 1 min
s1 = unlist(s1)
vals = na.omit(s1)
sign_prop = length(vals)/sum(!is.na(slp_s1))
pasit = sum(vals > 0)/length(vals)
cat("There are", round(sign_prop*100), "% of significant slopes", 
    "(5% occuring by chance).")
cat(" From them,", round(pasit*100), "% are positive.")

shcavm$PHYSIOG %>% unique()



