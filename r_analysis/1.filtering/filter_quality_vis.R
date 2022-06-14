library(raster)
library(tictoc)
par(mar = c(0.5,0.5,0.5,0.5))

dir.create("~/data/modis/snowfree_qa_filtered/")
dir.create("~/data/modis/snowfree_qa_filtered/vis/")
dir.create("~/data/modis/snowfree_qa_filtered/nir/")
dir.create("~/data/modis/snowfree_qa_filtered/sw/")

HOME_QA2 = "~/data/qa2/"
qa2to10 = function(year, band, homedir = HOME_QA2) {
  r = raster(paste0(homedir,"cropped_mosaics_",band,"/",year,"_",band,".tif"))
  r[r<3] = 0
  r[r >= 3] = 1
  return(r)
}



fnames = list.files("~/data/modis/snowfree_vis/",
                    full.names = T)
for (fname in fnames) {
  bn = basename(fname)
  print(bn)
  yr = as.numeric(substr(bn, 1,4))
  tic()
  r = qa2to10(yr, "band1")
  ra = qa2to10(yr, "band3")
  rb = qa2to10(yr, "band4")
  r[ra == 1] = 1
  r[rb == 1] = 1
  
  im = raster(fname)
  im[r==1] = NA
  writeRaster(im, paste0(
    "~/data/modis/snowfree_qa_filtered/vis/", bn),
    overwrite = T)
  toc()
  
}
