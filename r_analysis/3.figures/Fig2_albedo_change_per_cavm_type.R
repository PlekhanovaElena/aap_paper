library(ggplot2)
library(reshape)
library(RColorBrewer)
theme_set(theme_bw())
sw = read.csv("~/data/regression_slopes_061/stats_shortwave.csv")
dv = read.csv("~/data/regression_slopes_061/stats_vis.csv")
dn = read.csv("~/data/regression_slopes_061/stats_nir.csv")

dsave = sw
dsave$vis_mean_slope = dv$mean_slp
dsave$vis_std_slope = dv$std_slp
dsave$nir_mean_slope = dn$mean_slp
dsave$nir_std_slope = dn$std_slp
dsave$vis_mean_val = dv$mean_val
dsave$vis_std_val = dv$std_val
dsave$nir_mean_val = dn$mean_val
dsave$nir_std_val = dn$std_val

dat = as.data.frame(rbind(sw, dv, dn))
dat$band = c(rep("Shortwave", nrow(sw)), 
             rep("VIS", nrow(dv)), 
             rep("NIR+SWIR", nrow(dn)))
#dat$band = as.factor(dat$band)
dat$band = factor(dat$band, levels = c("Shortwave", "VIS", "NIR+SWIR"))

dat$veg_group = substr(dat$vtype, 1, 1)

ggplot(dat, aes(x=vtype, y=mean_slp, color=band)) + 
  geom_hline(yintercept = 0, col = "grey") +
  scale_color_brewer(palette = "Dark2") +
  geom_pointrange(aes(ymin=mean_slp-std_slp, ymax=mean_slp+std_slp),
                  position = position_dodge(0.4), shape = 18) +
  xlab("CAVM vegetation type") + ylab("mean albedo slope")




ggplot(dat, aes(x=vtype, y=mean_val, color=band)) + 
  geom_hline(yintercept = 0, col = "grey") +
  scale_color_brewer(palette = "Dark2") +
  geom_pointrange(aes(ymin=mean_val-std_val, ymax=mean_val+std_val),
                  position = position_dodge(0.4), shape = 18) +
  xlab("CAVM vegetation type") + ylab("mean albedo")


dat$percent_change_per_year = dat$mean_slp*100/dat$mean_val
dat$shape = 24
dat$shape[dat$mean_slp < 0] = 25
#dat$shape[abs(dat$percent_change_per_year) < 0.1] = NA
ggplot(dat, aes(x=vtype, color=band, fill = band)) + 
  scale_color_manual(values = brewer.pal(3, "Dark2")[c(1,3,2)]) +
  scale_fill_manual(values = brewer.pal(3, "Dark2")[c(1,3,2)]) +
  
  geom_linerange(aes(y=mean_val, ymin=mean_val-std_val, ymax=mean_val+std_val),
                 position = position_dodge(0.4), 
                 alpha = 0.4, size = 1.5) +
  
  geom_pointrange(aes(y=mean_val-mean_slp*10, ymin=mean_val-mean_slp*10-std_slp, 
                      ymax=mean_val-mean_slp*10+std_slp),
                  position = position_dodge(0.4), shape = 20, size = 0.4) +
  

  
  xlab("CAVM vegetation type") + ylab("albedo") +
  geom_pointrange(aes(y=mean_val+mean_slp*10, ymin=mean_val+mean_slp*10-std_slp, 
                      ymax=mean_val+mean_slp*10+std_slp),
                  position = position_dodge(0.4), shape = dat$shape, size = 0.25) +
  theme(panel.grid.major = element_line(colour = "lightgrey"),
        panel.grid.minor = element_blank())

gg = ggplot(dat, aes(x=vtype, color=band, fill = band)) + 
  scale_color_manual(values = brewer.pal(3, "Dark2")[c(1,3,2)]) +
  scale_fill_manual(values = brewer.pal(3, "Dark2")[c(1,3,2)]) +
  
  geom_linerange(aes(y=mean_val, ymin=mean_val-std_val, ymax=mean_val+std_val),
                 position = position_dodge(0.4), 
                 alpha = 0.4, size = 1.5) +
  
  geom_pointrange(aes(y=mean_val-mean_slp*10, ymin=mean_val-mean_slp*10-std_slp, 
                      ymax=mean_val-mean_slp*10+std_slp),
                  position = position_dodge(0.4), shape = 20, size = 0.4) +
  
  
  
  xlab("CAVM vegetation type") + ylab("albedo") +
  geom_pointrange(aes(y=mean_val+mean_slp*10, ymin=mean_val+mean_slp*10-std_slp, 
                      ymax=mean_val+mean_slp*10+std_slp),
                  position = position_dodge(0.4), shape = dat$shape, size = 0.25)


gg

gg + theme(panel.grid.major = element_line(colour = "lightgrey", size = 0.2),
        panel.grid.minor = element_blank())
 

ggplot(dat, aes(x=vtype, y=mean_slp*100/mean_val, color=band)) + 
  geom_hline(yintercept = 0, col = "grey") +
  scale_color_brewer(palette = "Dark2") +
  geom_pointrange(aes(ymin=mean_slp-std_slp, ymax=mean_slp+std_slp),
                  position = position_dodge(0.4), shape = 18) +
  xlab("CAVM vegetation type") + ylab("mean albedo slope")


ggplot(dat, aes(x=veg_group, y=mean_slp, color=band)) + 
  geom_hline(yintercept = 0, col = "grey") +
  geom_pointrange(aes(ymin=mean_slp-std_slp, ymax=mean_slp+std_slp))

mlt = melt(sw, id.vars = c())

ggplot(sw, aes(x=vtype, y=mean_slp, color=veg_group)) + 
  geom_hline(yintercept = 0, col = "grey") +
  geom_pointrange(aes(ymin=mean_slp-std_slp, ymax=mean_slp+std_slp))


