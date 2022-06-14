library(raster)
library(rasterVis)
library("RColorBrewer")
library("viridis")  
par(mar = rep(0.5, 4))

library(rgdal)
cavm_border = readOGR("./data/cavm/cp_coast_la_shp/cp_coast_la.shp")

restif = "/home/cluster/eplekh/scratch/climate_data/15Jul/0.1/result_tiffs/"

mdif = raster(paste0(
"/home/cluster/eplekh/scratch/climate_data/15Jul/0.1/result_tiffs/
  dif_betw_modis_and_climate_models.tif"))

msd = raster(
"/home/cluster/eplekh/scratch/climate_data/15Jul/0.1/result_tiffs/
  standard_deviation_climate_models.tif")


ex = raster(
  "/home/cluster/eplekh/data/regression_slopes_061/sens_pvals_sw.tif")
ex = resample(ex, mdif)


mdif[mdif > 0.2] = 0.2
mdif[mdif < -0.2] = -0.2

mdif[is.na(ex)] = NA
msd[is.na(ex)] = NA

cmm = raster(
  "/home/cluster/eplekh/scratch/climate_data/15Jul/0.1/result_tiffs/
  mean_climate_models.tif") 

mdm = mean(stack(list.files(
  "/home/cluster/eplekh/data/modis/snowfree_qa_filtered/sw/", full.names = T)[1:15]), 
  na.rm = T)
cmslp = raster(
  "/home/cluster/eplekh/scratch/climate_data/15Jul/0.1/result_tiffs/slope_ensemble.tif") 
cmpval = raster(
  "/home/cluster/eplekh/scratch/climate_data/15Jul/0.1/result_tiffs/pval_ensemble.tif") 

mdslp = raster(
  "~/data/regression_slopes_061/slope_sw2014.tif") 
mdpval = raster(
  "~/data/regression_slopes_061/pval_sw2014.tif") 
mdm = mdm/1000
mdslp = mdslp/1000

#mdm = aggregate(mdm, 8)
#mdm = resample(mdm, mdslp)
rel_mdslp = (mdslp/mdm)*100
rel_cmslp = (cmslp/cmm)*100

## align extents of the 2 rasters
plot(rel_mdslp)
plot(rel_cmslp)
rel_mdslp = projectRaster(rel_mdslp, ex, method = "ngb")
rel_cmslp = projectRaster(rel_cmslp, ex, method = "ngb")
mdpval = projectRaster(mdpval, ex, method = "ngb")
cmpval = projectRaster(cmpval, ex, method = "ngb")

rel_mdslp[is.na(ex)] = NA
rel_cmslp[is.na(ex)] = NA

pb = function(x) { #proportion of pixels beyond +-2 percent
  vls = na.omit(values(x))
  pr = sum(vls > 2)/length(vls)
  return(round(pr*100, 2))
}

pb(rel_cmslp)

rel_cmslp[rel_cmslp > 2] = 2
rel_cmslp[rel_cmslp < -2] = -2
rel_mdslp[rel_mdslp > 2] = 2
rel_mdslp[rel_mdslp < -2] = -2

plot(rel_mdslp)
plot(rel_cmslp)
srel_mdslp = rel_mdslp
srel_mdslp[mdpval > 0.05] = NA

srel_cmslp = rel_cmslp
srel_cmslp[cmpval > 0.05] = NA

#those are missing
rsdslp = raster( 
            "/home/cluster/eplekh/scratch/climate_data/15Jul/0.1/result_tiffs/rel_sd_cm_slopes.tif")

plot(rsdslp)

sdslp = raster("/home/cluster/eplekh/scratch/climate_data/15Jul/0.1/result_tiffs/sd_cm_slopes.tif")


output_folder = "/home/cluster/eplekh/data/data_for_figures/Fig1/"
writeRaster(mdif, paste0(output_folder, "mdif", ".tif"), overwrite = T)
writeRaster(msd, paste0(output_folder, "msd", ".tif"), overwrite = T)
writeRaster(rel_mdslp, paste0(output_folder, "rel_mdslp", ".tif"), overwrite = T)
writeRaster(rel_cmslp, paste0(output_folder, "rel_cmslp", ".tif"), overwrite = T)  
writeRaster(rsdslp, paste0(output_folder, "rsdslp", ".tif"), overwrite = T)
writeRaster(srel_mdslp, paste0(output_folder, "srel_mdslp", ".tif"), overwrite = T)
writeRaster(srel_cmslp, paste0(output_folder, "srel_cmslp", ".tif"), overwrite = T)



pal=colorRampPalette(c("#984EA3",  "#009FFF", "white", "#FF7F00", "#A65628"))
pal=colorRampPalette(c("#8018b8",  "#4d9fe7", "#eaeaea", "#f5812c", "#b82418"))
#299bfa 49acff

# Difference and std absolute values
pal=colorRampPalette(c("#8018b8",  "#4d9fe7", "#eaeaea", "#f5812c", "#b82418"))
two_plots = stack(mdif, msd)
names(two_plots) = c("dif", "std")

levelplot(two_plots, 
          pretty = T,
          colorkey=list(space='bottom', # plot legend at bottom
                        labels=list(at=(-2:2)/10, font=3)),  # legend ticks 
          scales=list(draw=FALSE),            # suppress axis labels
          col.regions=pal,   
          at=c(-0.2,-0.15,-0.1,-0.05,-0.01,0.01,0.05,0.1,0.15,0.2)) +
  latticeExtra::layer(sp.polygons(cavm_border["COAST_"], lwd = 0.4))

# Slopes for MODIS and models
pal=colorRampPalette(c("#8018b8",  "#4d9fe7", "#eaeaea", "#f5812c", "#b82418"))
two_plots = stack(rel_mdslp, rel_cmslp)
names(two_plots) = c("smodis", "scm")

levelplot(two_plots, 
          pretty = T,
          colorkey=list(space='bottom', # plot legend at bottom
                        labels=list(at=(-2:2),
                                    labels =paste0(c(-2:2), "%"), font=3)),  # legend ticks 
          scales=list(draw=FALSE),            # suppress axis labels
          col.regions=pal,   
          at=c(-2,-1,-0.1,0.1,1,2)) +
  latticeExtra::layer(sp.polygons(cavm_border["COAST_"], lwd = 0.4))

# SD of the climate models slopes
rsdslp = projectRaster(rsdslp, ex, method = "ngb")
rsdslp[rsdslp>2] = 2
rsdslp[is.na(ex)] = NA
pal=colorRampPalette(c("#eaeaea", "#f5812c", "#b82418"))
two_plots = stack(rsdslp, rsdslp)
names(two_plots) = c("relsdscm", "relsdscm2")

levelplot(two_plots, 
          pretty = T,
          colorkey=list(space='bottom', # plot legend at bottom
                        labels=list(at=(0:4)/2, font=3)),  # legend ticks 
          scales=list(draw=FALSE),            # suppress axis labels
          col.regions=pal,   
          at=(0:8)/4) +
  latticeExtra::layer(sp.polygons(cavm_border["COAST_"], lwd = 0.4))


# Significant slopes for MODIS and models
pal=colorRampPalette(c("#8018b8",  "#4d9fe7", "#eaeaea", "#f5812c", "#b82418"))
two_plots = stack(srel_mdslp, srel_cmslp)
names(two_plots) = c("smodis", "scm")

levelplot(two_plots, 
          pretty = T,
          colorkey=list(space='bottom', # plot legend at bottom
                        labels=list(at=(-2:2),
                                    labels =paste0(c(-2:2), "%"), font=3)),  # legend ticks 
          scales=list(draw=FALSE),            # suppress axis labels
          col.regions=pal,   
          at=c(-2,-1,-0.1,0.1,1,2)) +
  latticeExtra::layer(sp.polygons(cavm_border["COAST_"], lwd = 0.4))



pal=colorRampPalette(c("#984EA3",  "#2832C2", "#009FFF", "white",
                       "#FF7F00", "#E41A1C", "#A65628"))
pal=colorRampPalette(c("#8018b8",  "#4d9fe7","#eaeaea","#f5812c", "#b82418"))





two_plots = stack(mdif, msd)
names(two_plots) = c("dif", "std")

levelplot(two_plots, 
          pretty = T,
          colorkey=list(space='bottom', # plot legend at bottom
                        labels=list(at=(-2:2)/10, font=3)),  # legend ticks 
          scales=list(draw=FALSE),            # suppress axis labels
          col.regions=pal,   
          at=c(-0.2,-0.1,-0.015,0.015,0.1,0.2)) #+
#latticeExtra::layer(sp.polygons(cavm_border["COAST_"], lwd = 0.5))
four_plots_p = stack(ssa, ss1, ss2, ss3)
names(four_plots_p) = c("MODIS", "CNRM-CM6", "IPSL-CM6", "MIROC6_0")



d1 = m1 - ma
d2 = m2 - ma
d3 = m3 - ma
pal=colorRampPalette(c("#984EA3", "#377EB8", "#A6CEE3", 
                       "lightgray",
                       "#FF7F00", "#E41A1C", "#A65628"))

three_plots = stack(d1, d2, d3)
three_plots[three_plots > 0.2] = 0.2
three_plots[three_plots < -0.2] = -0.2

names(three_plots) = c("CNRM-CM6", "IPSL-CM6", "MIROC6_0")

pal=colorRampPalette(c("#984EA3",  "#2832C2", "#009FFF", "lightgray",
                       "#FF7F00", "#E41A1C", "#A65628"))

levelplot(three_plots, 
          pretty = T,
          colorkey=list(space='bottom', # plot legend at bottom
                        labels=list(at=(-2:2)/10, font=3)),  # legend ticks 
          scales=list(draw=FALSE),            # suppress axis labels
          col.regions=pal,   
          at=seq(-0.2, 0.2, len=8)) +
  latticeExtra::layer(sp.polygons(cavm_border["COAST_"], lwd = 0.5))


pal=colorRampPalette(c(viridis_pal()(10), "#FEC44F", "#FE9929", 
                       "#EC7014", "#CC4C02"))

levelplot(stack(ma, ma), 
          pretty = T,
          colorkey=list(space='bottom', # plot legend at bottom
                        labels=list(at=(0:5)/10, font=3)),  # legend ticks 
          scales=list(draw=FALSE),            # suppress axis labels
          col.regions=pal,   
          at=seq(0, 0.5, len=11)) +
  latticeExtra::layer(sp.polygons(cavm_border["COAST_"], lwd = 0.5))




two_plots = stack(vslp, sign_slp)
names(two_plots) = c("all_slopes", "significant_slopes")

cls = brewer.pal(n = 9, name = "Set1")
pal=colorRampPalette(c("#984EA3", "#377EB8", "#A6CEE3", 
                       "#FF7F00", "#E41A1C", "#A65628"))


levelplot(two_plots, 
          pretty = T,
          colorkey=list(space='bottom', # plot legend at bottom
                        labels=list(at=(-3:3)/1000, font=3)),  # legend ticks 
          scales=list(draw=FALSE),            # suppress axis labels
          col.regions=pal,   
          at=seq(-0.003, 0.003, len=7))


st = stack(list.files(
  "/home/cluster/eplekh/scratch/climate_data/15Jul/0.1/result_tiffs/slopes_per_model/slope/",
  full.names = T)) 
stm = stack(list.files(
  "/home/cluster/eplekh/scratch/climate_data/15Jul/0.1/average_albedo/",
  full.names = T)) 
relst = (st/stm)*100
rsdslp = calc(relst, fun = sd, na.rm = T)
writeRaster(rsdslp, 
"/home/cluster/eplekh/scratch/climate_data/15Jul/0.1/result_tiffs/rel_sd_cm_slopes.tif")


#st = stack(list.files(
#  "/home/cluster/eplekh/scratch/climate_data/15Jul/0.1/result_tiffs/slopes_per_model/slope/",
#  full.names = T)) 
#sdslp = calc(st, fun = sd, na.rm = T)
#writeRaster(sdslp, 
#"/home/cluster/eplekh/scratch/climate_data/15Jul/0.1/result_tiffs/sd_cm_slopes.tif")
#display.brewer.pal(n = 9, name = "Set1")
