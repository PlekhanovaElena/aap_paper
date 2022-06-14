library(raster)
library(tictoc)
par(mar = c(0.5,1,0.5,0.5))



## Preparing Climate models

modnames = list.dirs("/home/cluster/eplekh/scratch/climate_data/15Jul/0.1/tiffs/",
                   full.names = F)[-1]



orig_path <- "/home/cluster/eplekh/scratch/climate_data/15Jul/0.1/tiffs/"
output_path = "/home/cluster/eplekh/scratch/climate_data/15Jul/0.1/albedo/"
dir.create(output_path)

modname = modnames[2]
yr = 2000
flname <- paste0(yr, "_", modname, "_") 
fname_d <- paste0(orig_path, modname, "/", flname, "snc_", ".tiff")
d <- raster(fname_d)
plot(d)
plot(d>=1)

for (modname in modnames) {
  output_dir = paste0(output_path, modname)
  dir.create(output_dir, showWarnings = FALSE)
  
  for (yr in c(2000:2014)) {
    flname <- paste0(yr, "_", modname, "_") 
    fname_d <- paste0(orig_path, modname, "/", flname, "rsds", ".tiff")
    d <- raster(fname_d)
    fname_u <- paste0(orig_path, modname, "/", flname, "rsus", ".tiff")
    u <- raster(fname_u)
    a = u/d
    fname_s <- paste0(orig_path, modname, "/", flname, "snc_", ".tiff")
    s <- raster(fname_s)
    a[s>=1] = NA
    #ar = resample(a, agg_alb, method = "ngb")
    writeRaster(a, paste0(output_dir, "/", yr, ".tif"), overwrite=TRUE)
  }
}


plot(a)

stm = stack(list.files(paste0(output_path,modnames[1]), full.names = T))
firstmod = mean(stm) # average across years albedo for this model
writeRaster(firstmod, paste0(
  "/home/cluster/eplekh/scratch/climate_data/15Jul/0.1/average_albedo/",
  modnames[1], ".tif"))

dir.create("/home/cluster/eplekh/scratch/climate_data/15Jul/0.1/average_albedo/",
           showWarnings = FALSE)

for (modname in modnames) {
  stm = stack(list.files(paste0(output_path,modname), full.names = T))
  out = mean(stm) # average across years albedo for this model
  tic()
  out = resample(out, firstmod, method = "ngb")
  toc()
  writeRaster(out, paste0(
    "/home/cluster/eplekh/scratch/climate_data/15Jul/0.1/average_albedo/",
    modname, ".tif"), overwrite = T)
}



st = stack(list.files(
  "/home/cluster/eplekh/scratch/climate_data/15Jul/0.1/average_albedo/",
  full.names = T))

mcm = mean(st, na.rm = T)
sdcm = calc(st, fun = sd, na.rm = T)
plot(mcm)
plot(sdcm)

par(mar = c(0.1,0.1,0.1,0.1))
plot(st)



## Preparing MODIS

stmod = stack(list.files(
  "/home/cluster/eplekh/data/modis/snowfree_qa_filtered/sw/",
                         full.names = T)[1:15])

#stmod = stack(list.files("/home/cluster/eplekh/scratch/modis_july/mean_year/",
#                      full.names = T))

modf = mean(stmod, na.rm = T)
modr = projectRaster(modf, mcm, res = res(mcm))
mod = resample(modr, mcm)
plot(mod)
mod = mod/1000
dif =  mcm - mod

cm_minus_mod = st - mod
stda = calc(cm_minus_mod, fun = sd, na.rm = T)
plot(stda)

plot(stda)

writeRaster(mcm, paste0(
  "/home/cluster/eplekh/scratch/climate_data/15Jul/0.1/result_tiffs/
  mean_climate_models.tif"))
writeRaster(mod, paste0(
  "/home/cluster/eplekh/scratch/climate_data/15Jul/0.1/result_tiffs/
  mean_modis.tif"))
writeRaster(stda, paste0(
  "/home/cluster/eplekh/scratch/climate_data/15Jul/0.1/result_tiffs/
  standard_deviation_climate_models.tif"))
writeRaster(dif, paste0(
  "/home/cluster/eplekh/scratch/climate_data/15Jul/0.1/result_tiffs/
  dif_betw_modis_and_climate_models.tif"), overwrite = T)

# Preparation of standard deviation data for slopes

stk = stack(list.files(
  "/home/cluster/eplekh/scratch/climate_data/15Jul/0.1/result_tiffs/slopes_per_model/",
  full.names = T)[7:12])
stk
plot(stk)

std_modslp = calc(stk, fun = sd, na.rm = T)

mcm = raster("/home/cluster/eplekh/scratch/climate_data/15Jul/0.1/result_tiffs/
  mean_climate_models.tif")

rel_std_modslp = std_modslp*100/mcm

plot(rel_std_modslp)
writeRaster( rel_std_modslp,
  "/home/cluster/eplekh/scratch/climate_data/15Jul/0.1/result_tiffs/rel_sd_cm_slopes.tif")

writeRaster( std_modslp,
  "/home/cluster/eplekh/scratch/climate_data/15Jul/0.1/result_tiffs/sd_cm_slopes.tif")












