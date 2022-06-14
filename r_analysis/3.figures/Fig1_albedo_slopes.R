# load libraries
library(tictoc)
library(ggplot2)

library(raster)
library(rgdal)
cavm_border = readOGR("./data/cavm/cp_coast_la_shp/cp_coast_la.shp")
cavmr = raster( "./data/cavm/cavm_raster_resampled3.tif")


## read
stats_dir = "./data/regression_slopes_061/" # input folder

slps = raster(paste0(stats_dir, "sens_slope_vis.tif"))
pv1 = raster(paste0(stats_dir, "sens_pvals_vis.tif"))
slps[cavmr >= 90] = NA
ss = slps
ss[pv1 > 0.05] = NA
slpv = slps
ssv = ss

median(values(slps), na.rm = T)*22/1000
median(values(ss), na.rm = T)*22/1000

slps = raster(paste0(stats_dir, "sens_slope_nir.tif"))
pv1 = raster(paste0(stats_dir, "sens_pvals_nir.tif"))
slps[cavmr >= 90] = NA
ss = slps
ss[pv1 > 0.05] = NA
slpn = slps
ssn = ss

median(values(ss), na.rm = T)*22/1000

slps = raster(paste0(stats_dir, "sens_slope_sw.tif"))
pv1 = raster(paste0(stats_dir, "sens_pvals_sw.tif"))
slps[cavmr >= 90] = NA
ss = slps
ss[pv1 > 0.05] = NA



ms = mean(stack(list.files(
  "/home/cluster/eplekh/data/modis/snowfree_shortwave/", full.names = T)), 
  na.rm = T)
mv = mean(stack(list.files(
  "/home/cluster/eplekh/data/modis/snowfree_vis/", full.names = T)), 
  na.rm = T)
mn = mean(stack(list.files(
  "/home/cluster/eplekh/data/modis/snowfree_nir/", full.names = T)), 
  na.rm = T)

rm(cavmr)

rslpn = (slpn/mn)*100
rslpv = (slpv/mv)*100
rslps = (slps/ms)*100

rssn = (ssn/mn)*100
rssv = (ssv/mv)*100
rsss = (ss/ms)*100
output_folder = "/home/cluster/eplekh/data/data_for_figures/Fig1/"
writeRaster(rslpn, paste0(output_folder, "rslpn", ".tif"), overwrite=T)
writeRaster(rslpv, paste0(output_folder, "rslpv", ".tif"), overwrite=T)
writeRaster(rslps, paste0(output_folder, "rslps", ".tif"), overwrite=T)  
writeRaster(rssn, paste0(output_folder, "rssn", ".tif"), overwrite=T)
writeRaster(rssv, paste0(output_folder, "rssv", ".tif"), overwrite=T)
writeRaster(rsss, paste0(output_folder, "rsss", ".tif"), overwrite=T)

#####################################################################

input_folder = "/home/cluster/eplekh/data/data_for_figures/Fig1/"

rslpn = raster(paste0(input_folder, "rslpn", ".tif"))
rslpv = raster(paste0(input_folder, "rslpv", ".tif"))
rslps = raster(paste0(input_folder, "rslps", ".tif"))
rssn = raster(paste0(input_folder, "rssn", ".tif"))  
rssv = raster(paste0(input_folder, "rssv", ".tif"))
rsss = raster(paste0(input_folder, "rsss", ".tif"))
## visualize
library(rasterVis)






x = 2
rslpn[rslpn > x] = x
rslpv[rslpv > x] = x
rslps[rslps > x] = x
rssn[rssn > x] = x
rssv[rssv > x] = x
rsss[rsss > x] = x

rslpn[rslpn < -x] = -x
rslpv[rslpv < -x] = -x
rslps[rslps < -x] = -x
rssn[rssn < -x] = -x
rssv[rssv < -x] = -x
rsss[rsss < -x] = -x

six_plots = stack(rslps, rslpv, rslpn, rsss, rssv, rssn)

names(six_plots) = c("Shortwave", "VIS", "NIR+SWIR",
                     "si1", "si2","si3")


#pal=colorRampPalette(c("#984EA3", "#377EB8", "#A6CEE3", "white",
#                       "#FF7F00", "#E41A1C", "#A65628"))
pal=colorRampPalette(c("#8018b8",  "#4d9fe7", "#eaeaea", "#f5812c", "#b82418"))
levelplot(six_plots, 
          #margin=FALSE,                       # suppress marginal graphics
          pretty = T,
          colorkey=list(
            space='bottom',                   # plot legend at bottom
            labels=list(at=(-x:x),
                        labels =paste0(c(-x:x), "%"), font=5)),  # legend ticks and labels 
    
          scales=list(draw=FALSE),            # suppress axis labels
          col.regions=pal,   
          at=c(-2,-1,-0.1,0.1,1,2)) +
  latticeExtra::layer(sp.polygons(cavm_border["COAST_"], lwd = 0.4))



six_plots = stack(slps, slpv, slpn, ss, ssv, ssn)
names(six_plots) = c("Shortwave", "Visible", "Near-infrared",
                     "si1", "si2","si3")


pal=colorRampPalette(c("#984EA3", "#377EB8", "#A6CEE3", 
                       "#FF7F00", "#E41A1C", "#A65628"))
levelplot(six_plots, 
          #margin=FALSE,                       # suppress marginal graphics
          pretty = T,
          colorkey=list(
            space='bottom',                   # plot legend at bottom
            labels=list(at=(-3:3)/1000, font=3)  # legend ticks and labels 
          ),    

          scales=list(draw=FALSE),            # suppress axis labels
          col.regions=pal,   
          at=seq(-0.003, 0.003, len=7)) +
  latticeExtra::layer(sp.polygons(cavm_border["COAST_"], lwd = 0.4))




levelplot(six_plots, 
          #margin=FALSE,                       # suppress marginal graphics
          pretty = T,
          colorkey=list(
            space='bottom',                   # plot legend at bottom
            labels=list(at=(-3:3)/1000, font=3)  # legend ticks and labels 
          ),    
          
          scales=list(draw=FALSE),            # suppress axis labels
          col.regions=rainbow(16, rev = T),   
          at=seq(-0.003, 0.003, len=16)) +
  latticeExtra::layer(sp.polygons(cavm_border["COAST_"], lwd = 0.5))




ggplot() + layer(sp.polygons(cavm_border["COAST_"]), stat="identity", 
               position = "identity", geom="line")


cst = cavm_border["COAST_"]
ggplot() + layer(cavm_border["COAST_"]@polygons, stat="identity", 
                 position = "identity", geom="line")

ggplot() + geom_polygon(cst, aes(x = ))
ggplot(sp.polygons(cst)) + layer(aes(x = lat, y = lon), stat="identity", 
                    position = "identity", geom="line")



levelplot(slps, 
          #margin=FALSE,                       # suppress marginal graphics
          pretty = T,
          colorkey=list(
            space='bottom',                   # plot legend at bottom
            labels=list(at=(-3:3)/1000, font=3)  # legend ticks and labels 
          ),    
          
          scales=list(draw=FALSE),            # suppress axis labels
          col.regions=pal,   
          at=seq(-0.003, 0.003, len=7)) + 
  layer(cavm_border["COAST_"])


levelplot(s, contour=TRUE, par.settings = list(axis.line = list(col = "transparent"), 
                                               strip.background = list(col = 'transparent'), 
                                               strip.border = list(col = 'transparent')), 
          scales = list(col = "black"))`

four_plots = stack(slpv, ssv, slpn, ssn)
names(four_plots) = c("VIS_slopes", "VIS_significant_slopes",
                     "NIR_slopes", "NIR_significant_slopes")

#plot(cavm, lwd = 0.01, col = "grey90")
levelplot(four_plots, 
          #margin=FALSE,                       # suppress marginal graphics
          pretty = T,
          colorkey=list(
            space='bottom',                   # plot legend at bottom
            labels=list(at=(-3:3)/1000, font=3)  # legend ticks and labels 
          ),    

          scales=list(draw=FALSE),            # suppress axis labels
          col.regions=rainbow(16, rev = T),   
          at=seq(-0.003, 0.003, len=16))



two_plots = stack(slps, ss)
names(two_plots) = c("all_slopes", "significant_slopes")

#plot(cavm, lwd = 0.01, col = "grey90")
levelplot(two_plots, 
          #margin=FALSE,                       # suppress marginal graphics
          pretty = T,
          colorkey=list(
            space='bottom',                   # plot legend at bottom
            labels=list(at=(-3:3)/1000, font=3)  # legend ticks and labels 
          ),    

          scales=list(draw=FALSE),            # suppress axis labels
          col.regions=rainbow(16, rev = T),   
          at=seq(-0.003, 0.003, len=16))


###############################################################################
####################  Calculating statistics of overlap  #####################

pb = function(x) { #proportion of pixels beyond +-2 percent
  vls = na.omit(values(x))
  pr = sum(vls < -2)/length(vls)
  return(round(pr*100, 2))
}

pb(rslpn)
pb(rslpv)
pb(rslps)
pb(rssn)
pb(rssv)
pb(rsss)

av = cavmr
av[av>=90] = NA
av = values(av)
av = av[!is.na(av)]
a = length(av)
percstats = function(all_slopes, sign_slopes, all_notna = a) {
  a = all_notna
  bnv = values(all_slopes)
  cnv = values(sign_slopes)
  
  bnv = bnv[!is.na(bnv)]
  cnv = cnv[!is.na(cnv)]
  cn = length(cnv)
  bn = length(bnv)
  cat("percent significant slopes\n")
  print(round(cn/a, 2)*100)
  cat("percent positive significant slopes\n")
  print(round(sum(cnv > 0)/a, 2)*100)
  cat("median slope\n")
  print(round(median(bnv), 2))
  cat("median significant slope\n")
  print(round(median(cnv), 2))
}

### NIR
percstats(rslpn, rssn)

### VIS
percstats(rslpv, rssv)

### SW
percstats(rslps, rsss)


### Overlap
ovrlap = rssn + rssv*0
ov = values(ovrlap)
ov = ov[!is.na(ov)]
cat("percent overlap significant slopes NIR and VIS\n")
print(round(length(ov)/a, 3)*100)

ovrlap = rssn*0 + rsss
ov = values(ovrlap)
ov = ov[!is.na(ov)]
cat("percent overlap significant slopes SW and NIR\n")
print(round(length(ov)/a, 3)*100)

ovrlap = rssv*0 + rsss
ov = values(ovrlap)
ov = ov[!is.na(ov)]
cat("percent overlap significant slopes SW and VIS\n")
print(round(length(ov)/a, 3)*100)

ovrlap = rssv*0 + rsss + rssn
ov = values(ovrlap)
ov = ov[!is.na(ov)]
cat("percent overlap significant slopes NIR and VIS and SW\n")
print(round(length(ov)/a, 3)*100)



### NIR
bnv = values(rslpn)
cnv = values(rssn)

bnv = bnv[!is.na(bnv)]
cnv = cnv[!is.na(cnv)]
cn = length(cnv)
bn = length(bnv)
cn/bn
cn/a
sum(cnv > 0)/a
median(cnv)/1000

### VIS

slpvalsv = values(rslpv)
slpvalsv = slpvalsv[!is.na(slpvalsv)]
valsv = values(rssv)
valsv = valsv[!is.na(valsv)]
length(valsv)
print("VIS")
length(valsv)/length(slpvalsv)
sum(valsv > 0)/length(slpvalsv)
median(valsv)


### Overlap
ovrlap = rssn + rssv*0
plot(ovrlap, col = "black")
valso = values(ovrlap)
valso = valso[!is.na(valso)]
length(valso)/length(slpvalsn)
sum(valso > 0)/length(slpvalsn)

ovrlap = rssv + rssn*0
plot(ovrlap, col = "black")
valso2 = values(ovrlap)
valso2 = valso2[!is.na(valso2)]
length(valso2)/length(slpvalsn)
sum(valso2 > 0)/length(slpvalsn)

v = rssv
v[is.na(ovrlap)] = NA

plot(v, col = "black")
v[(v > 0) & (rssn > 0)] = NA 
plot(v, col = "black")


### Shortwave
slpvals = values(rslps)
slpvals = slpvals[!is.na(slpvals)]
vals = values(rsss)
vals = vals[!is.na(vals)]
length(vals)
length(vals)/length(slpvals)
sum(vals > 0)/length(slpvals)

median(vals)

### Overlap sw
ovrlap = rssn*0 + rslps
plot(ovrlap)
valso = values(ovrlap)
valso = valso[!is.na(valso)]
length(valso)/length(slpvalsn)

ovrlap = rssv*0 + rslps
plot(ovrlap)
valso = values(ovrlap)
valso = valso[!is.na(valso)]
length(valso)/length(slpvalsn)

ovrlap = rssv*0 + rslps + rssn
plot(ovrlap)
valso = values(ovrlap)
valso = valso[!is.na(valso)]
length(valso)/length(slpvalsn)

ovrlap = rslps
plot(ovrlap)
valso = values(ovrlap)
valso = valso[!is.na(valso)]
length(valso)/length(slpvalsn)








