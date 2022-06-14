library(raster)

broadband = "shortwave"

path_init = paste0("~/data/modis/cropped_mosaics_", broadband, "/")
path_output = paste0("~/data/modis/snowfree_", broadband, "/")
path_snow = "~/data/cropped_snow/"
path_snow_flag = "~/data/qa2/cropped_mosaics_snowflag/"


dir.create(path_output)

path_ims = list.files(path_init)
path_ims_snow = list.files(path_snow, full.names = T)
path_ims_flag = list.files(path_snow_flag, full.names = T)

for (i in 1:length(path_ims)) {
  print(i)
  s = raster(path_ims_snow[i])
  s[s>100] = NA
  s[s>1] = 100 # excluding pixels with snow contamination more than 1%
  
  r = raster(paste0(path_init, path_ims[i]))
  a = r
  a[a > 300] = 300
  s[is.na(s) & (a == 300)] = 100 # excluding values larger than 0.3
  
  f = raster(path_ims_flag[i]) # excluding snow flagged pixels
  s[f == 1] = 100
  
  r[s == 100] = NA
  
  writeRaster(r, paste0(path_output, path_ims[i]), overwrite = T)
}

broadband = "vis"

path_init = paste0("~/data/modis/cropped_mosaics_", broadband, "/")
path_output = paste0("~/data/modis/snowfree_", broadband, "/")
dir.create(path_output)

path_ims = list.files(path_init)
path_ims_snow = list.files(path_snow, full.names = T)

for (i in 1:length(path_ims)) {
  print(i)
  s = raster(path_ims_snow[i])
  s[s>100] = NA
  s[s>1] = 100
  
  r = raster(paste0(path_init, path_ims[i]))
  a = r
  a[a > 300] = 300
  s[is.na(s) & (a == 300)] = 100
  
  f = raster(path_ims_flag[i])
  s[f == 1] = 100
  
  r[s == 100] = NA
  
  writeRaster(r, paste0(path_output, path_ims[i]), overwrite = T)
}

broadband = "nir"

path_init = paste0("~/data/modis/cropped_mosaics_", broadband, "/")
path_output = paste0("~/data/modis/snowfree_", broadband, "/")
dir.create(path_output)

path_ims = list.files(path_init)
path_ims_snow = list.files(path_snow, full.names = T)

for (i in 1:length(path_ims)) {
  print(i)
  s = raster(path_ims_snow[i])
  s[s>100] = NA
  s[s>1] = 100
  
  r = raster(paste0(path_init, path_ims[i]))
  a = r
  a[a > 300] = 300
  s[is.na(s) & (a == 300)] = 100
  
  f = raster(path_ims_flag[i])
  s[f == 1] = 100
  
  r[s == 100] = NA
  
  writeRaster(r, paste0(path_output, path_ims[i]), overwrite = T)
}