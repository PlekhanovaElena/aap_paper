library(raster)
library(rasterVis)
library("viridis")  
library(rgdal)
cavm_border = readOGR("./data/cavm/cp_coast_la_shp/cp_coast_la.shp")

input_folder = "/home/cluster/eplekh/data/data_for_figures/Fig1/"
mdif = raster(paste0(input_folder, "mdif", ".tif"))
msd = raster(paste0(input_folder, "msd", ".tif"))
rel_mdslp = raster(paste0(input_folder, "rel_mdslp", ".tif"))
rel_cmslp = raster(paste0(input_folder, "rel_cmslp", ".tif"))  
rsdslp = raster(paste0(input_folder, "rsdslp", ".tif"))
srel_mdslp = raster(paste0(input_folder, "srel_mdslp", ".tif"))
srel_cmslp = raster(paste0(input_folder, "srel_cmslp", ".tif"))

stk = stack(list.files(
  "/home/cluster/eplekh/scratch/climate_data/15Jul/0.1/average_albedo/", 
  full.names = T))

cmm = calc(stk, mean, na.rm = T)
plot(cmm)
writeRaster(cmm, paste0(
  "/home/cluster/eplekh/scratch/climate_data/15Jul/0.1/result_tiffs/", 
  "mean_climate_models", ".tif"), overwrite = T)

cmm = raster(
  paste0(
    "/home/cluster/eplekh/scratch/climate_data/15Jul/0.1/result_tiffs/", 
    "mean_climate_models", ".tif")) 

stk = stack(list.files(
  "/home/cluster/eplekh/data/modis/snowfree_qa_filtered/sw/", 
  full.names = T)[1:15])

modm = calc(stk, mean, na.rm = T)
writeRaster(modm, paste0(input_folder, "mean_modis_15y", ".tif"))
modm = raster(paste0(input_folder, "mean_modis_15y", ".tif"))
modm = modm/1000

# Absolute albedo values
pal=colorRampPalette(c(viridis_pal()(10), "#FEC44F", "#FE9929", 
                       "#EC7014", "#CC4C02"))
#plot(modm)
modm = resample(modm, cmm)
cmm[cmm > 0.3] = 0.3
#cmm[is.na(modm)] = NA
#modm[is.na(cmm)] = NA
two_plots = stack(modm, cmm)
names(two_plots) = c("m1", "m2")
pal=colorRampPalette(c(viridis_pal()(7), "#FEC44F", "#FE9929", 
                       "#EC7014", "#CC4C02"))

levelplot(two_plots, 
          pretty = T,
          colorkey=list(space='bottom', # plot legend at bottom
                        labels=list(at=(0:3)/10, font=6, cex = 0.7)),  # legend ticks 
          scales=list(draw=FALSE),            # suppress axis labels
          col.regions=pal,   
          at=seq(0, 0.3, len=10)) +
  latticeExtra::layer(sp.polygons(cavm_border["COAST_"], lwd = 0.4))

# Absolute albedo values all models
pal=colorRampPalette(c(viridis_pal()(10), "#FEC44F", "#FE9929", 
                       "#EC7014", "#CC4C02"))

mods = stack(list.files("~/scratch/climate_data/15Jul/0.1/average_albedo/", 
                        full.names = T))


plot(mods)

CNRM_CM6 = mean(stack(list.files(
  "~/scratch/climate_data/15Jul/0.1/albedo/CNRM-CM6/", 
  full.names = T)), na.rm = T)
CNRM_ESM = mean(stack(list.files(
  "~/scratch/climate_data/15Jul/0.1/albedo/CNRM-ESM/", 
  full.names = T)), na.rm = T)
EC_Earth = mean(stack(list.files(
  "~/scratch/climate_data/15Jul/0.1/albedo/EC-Earth/", 
  full.names = T)), na.rm = T)
IPSL_CM6 = mean(stack(list.files(
  "~/scratch/climate_data/15Jul/0.1/albedo/IPSL-CM6/", 
  full.names = T)), na.rm = T)
MIROC6_2 = mean(stack(list.files(
  "~/scratch/climate_data/15Jul/0.1/albedo/MIROC6_2/", 
  full.names = T)), na.rm = T)
UKESM1_0 = mean(stack(list.files(
  "~/scratch/climate_data/15Jul/0.1/albedo/UKESM1-0/", 
  full.names = T)), na.rm = T)

av_folder = "~/scratch/climate_data/15Jul/0.1/average_albedo/"
writeRaster(CNRM_CM6, paste0(av_folder, "CNRM-CM6", ".tif"), overwrite = T)
writeRaster(CNRM_ESM, paste0(av_folder, "CNRM-ESM", ".tif"), overwrite = T)
writeRaster(EC_Earth, paste0(av_folder, "EC-Earth", ".tif"), overwrite = T)
writeRaster(IPSL_CM6, paste0(av_folder, "IPSL-CM6", ".tif"), overwrite = T)
writeRaster(MIROC6_2, paste0(av_folder, "MIROC6_2", ".tif"), overwrite = T)
writeRaster(UKESM1_0, paste0(av_folder, "UKESM1-0", ".tif"), overwrite = T)

av_mod = stack(list.files(av_folder, full.names = T))
plot(av_mod)
#modm = resample(modm, cmm)
av_mod[av_mod > 0.3] = 0.3

pal=colorRampPalette(c(viridis_pal()(7), "#FEC44F", "#FE9929", 
                       "#EC7014", "#CC4C02"))
names(av_mod) = c("CNRM-CM6", "CNRM-ESM", "EC-Earth", "IPSL-CM6", "MIROC6_2", "UKESM1-0")
levelplot(av_mod, 
          pretty = T,
          colorkey=list(space='bottom', # plot legend at bottom
                        labels=list(at=(0:3)/10, font=6, cex = 0.7)),  # legend ticks 
          scales=list(draw=FALSE),            # suppress axis labels
          col.regions=pal,   
          at=seq(0, 0.3, len=10)) +
  latticeExtra::layer(sp.polygons(cavm_border["COAST_"], lwd = 0.4))

# Difference and std absolute values
av_folder = "~/scratch/climate_data/15Jul/0.1/average_albedo/"
av_mod = stack(list.files(av_folder, full.names = T))
msd = calc(av_mod, sd, na.rm = T)
pal=colorRampPalette(c("#8018b8",  "#4d9fe7", "#eaeaea", "#f5812c", "#b82418"))
cmm = raster(
  paste0(
    "/home/cluster/eplekh/scratch/climate_data/15Jul/0.1/result_tiffs/", 
    "mean_climate_models", ".tif")) 
mdif = cmm - modm
mdif[mdif >= 0.2] = 0.2
mdif[mdif <= -0.2] = -0.2
plot(mdif)
two_plots = stack(mdif, msd)
names(two_plots) = c("dif", "std")

levelplot(two_plots, 
          pretty = T,
          colorkey=list(space='bottom', # plot legend at bottom
                        labels=list(at=(-2:2)/10, font=3)),  # legend ticks 
          scales=list(draw=FALSE),            # suppress axis labels
          col.regions=pal,   
          at=c(-0.2,-0.1,-0.01,0.01,0.1,0.2)) +
  latticeExtra::layer(sp.polygons(cavm_border["COAST_"], lwd = 0.3))

levelplot(two_plots, 
          pretty = T,
          colorkey=list(space='bottom', # plot legend at bottom
                        labels=list(at=(-1:1)/20, font=3)),  # legend ticks 
          scales=list(draw=FALSE),            # suppress axis labels
          col.regions=pal,   
          at=c(-0.1,-0.05,-0.01,0.01,0.05,0.1)) +
  latticeExtra::layer(sp.polygons(cavm_border["COAST_"], lwd = 0.4))

levelplot(two_plots, 
          pretty = T,
          colorkey=list(space='bottom', # plot legend at bottom
                        labels=list(at=(-1:1)/10, font=3)),  # legend ticks 
          scales=list(draw=FALSE),            # suppress axis labels
          col.regions=pal,   
          at=c(-0.1,-0.05,-0.01,0.01,0.05,0.1)) +
  latticeExtra::layer(sp.polygons(cavm_border["COAST_"], lwd = 0.4))

levelplot(two_plots, 
          pretty = T,
          colorkey=list(space='bottom', # plot legend at bottom
                        labels=list(at=(-2:2)/10, font=3)),  # legend ticks 
          scales=list(draw=FALSE),            # suppress axis labels
          col.regions=pal,   
          at=c(-0.2, -0.15, -0.1,-0.05, -0.01,0.01, 0.05,0.1,0.15,0.2)) +
  latticeExtra::layer(sp.polygons(cavm_border["COAST_"], lwd = 0.4))

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


# Difference and std relative values
pal=colorRampPalette(c("#8018b8",  "#4d9fe7", "#eaeaea", "#f5812c", "#b82418"))
rel_diff = rel_cmslp-rel_mdslp
plot(rel_diff)
rel_diff[rel_diff > 2] = 2
rel_diff[rel_diff < -2] = -2
two_plots = stack(rel_diff, rel_diff)
names(two_plots) = c("dif1", "dif2")

levelplot(two_plots, 
          pretty = T,
          colorkey=list(space='bottom', # plot legend at bottom
                        labels=list(at=(-2:2),
                                    labels =paste0(c(-2:2), "%"), font=3)),  # legend ticks 
          scales=list(draw=FALSE),            # suppress axis labels
          col.regions=pal,   
          at=c(-2,-1,-0.1,0.1,1,2)) +
  latticeExtra::layer(sp.polygons(cavm_border["COAST_"], lwd = 0.4))


rel_diff = rel_mdslp-rel_cmslp
two_plots = stack(rel_diff, rel_diff)
names(two_plots) = c("dif1", "dif2")

levelplot(two_plots, 
          pretty = T,
          colorkey=list(space='bottom', # plot legend at bottom
                        labels=list(at=(-4:4),
                                    labels =paste0(c(-4:4), "%"), font=3)),  # legend ticks 
          scales=list(draw=FALSE),            # suppress axis labels
          col.regions=pal,   
          at=c(-4,-3,-2,-1,-0.1,0.1,1,2,3,4))

# Bias LSM-MODIS abs
pal=colorRampPalette(c("#8018b8",  "#4d9fe7", "#eaeaea", "#f5812c", "#b82418"))
two_plots = stack(mdif, msd)
names(two_plots) = c("dif", "std")

levelplot(two_plots, 
          pretty = F,
          colorkey=list(space='bottom', # plot legend at bottom
                        labels=list(at=(-2:2)/10, font=1, cex=0.6)),  # legend ticks 
          scales=list(draw=FALSE),            # suppress axis labels
          col.regions=pal,   
          at=c(-0.2,-0.15,-0.1,-0.05,-0.01,0.01,0.05,0.1,0.15,0.2)) +
  latticeExtra::layer(sp.polygons(cavm_border["COAST_"], lwd = 0.4))


# Bias LSM-MODIS slopes
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
                        labels=list(at=(0:4)/2,
                                    labels =paste0((0:4)/2, "%"), font=6)),  # legend ticks 
          scales=list(draw=FALSE),            # suppress axis labels
          col.regions=pal,   
          at=(0:8)/4) +
  latticeExtra::layer(sp.polygons(cavm_border["COAST_"], lwd = 0.4))

# Difference and std absolute values
pal=colorRampPalette(c("#eaeaea", "#f5812c", "#b82418"))
two_plots = stack(msd, msd)
names(two_plots) = c("std1", "std")

levelplot(two_plots, 
          pretty = T,
          colorkey=list(space='bottom', # plot legend at bottom
                        labels=list(at=(0:1)/10, font=3)),  # legend ticks 
          scales=list(draw=FALSE),            # suppress axis labels
          col.regions=pal,   
          at=c(0,0.01,0.05,0.1)) +
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

# Four plots with slopes and significant
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


mdm = mean(stack(list.files(
  "/home/cluster/eplekh/data/modis_shortwave_snowfree/", full.names = T)[1:15]), 
  na.rm = T)

mdm = aggregate(mdm, fact = 10)
ma = mean(mdm, na.rm = T)
ma = ma/1000
plot(ma)
#ex = raster("/home/cluster/eplekh/scratch/regression_slopes_snowfree/sens_pvals_qa_nir.tif")
#ex = resample(ex, ma)
#ma[is.na(ex)] = NA
ma[ma > 0.3] = 0.3
plot(ma)
writeRaster(ma, paste0(input_folder, "mean_modis_15y", ".tif"), overwrite=TRUE)

# Absolute albedo values for all models
pal=colorRampPalette(c(viridis_pal()(10), "#FEC44F", "#FE9929", 
                       "#EC7014", "#CC4C02"))

mds = stack("~/scratch/climate_data/15Jul/0.1/average_albedo/")
modm = resample(modm, cmm)
cmm[cmm > 0.3] = 0.3
cmm[is.na(modm)] = NA
two_plots = stack(modm, cmm)
names(two_plots) = c("m1", "m2")

#levelplot(two_plots, pretty = T)

pal=colorRampPalette(c(viridis_pal()(7), "#FEC44F", "#FE9929", 
                       "#EC7014", "#CC4C02"))

levelplot(two_plots, 
          pretty = T,
          colorkey=list(space='bottom', # plot legend at bottom
                        labels=list(at=(0:3)/10, font=6, cex = 0.7)),  # legend ticks 
          scales=list(draw=FALSE),            # suppress axis labels
          col.regions=pal,   
          at=seq(0, 0.3, len=10)) +
  latticeExtra::layer(sp.polygons(cavm_border["COAST_"], lwd = 0.3))
