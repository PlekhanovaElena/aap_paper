library(raster)

# this script goes after filter_quality_vis.R and filter_quality_nir.R
fnames = list.files("/home/cluster/eplekh/data/modis/snowfree_shortwave/",
                    full.names = T)
for (fname in fnames) {
  bn = basename(fname)
  print(bn)
  yr = as.numeric(substr(bn, 1,4))
  v = raster(paste0("/home/cluster/eplekh/data/modis/snowfree_vis/", 
                    yr, "_vis.tif"))
  n = raster(paste0("/home/cluster/eplekh/data/modis/snowfree_nir/", 
                    yr, "_nir.tif"))
  im = raster(fname)
  im[is.na(v)] = NA
  im[is.na(n)] = NA
  writeRaster(im, paste0(
    "/home/cluster/eplekh/data/modis/snowfree_qa_filtered/sw/", bn),
    overwrite = T)
  
}
