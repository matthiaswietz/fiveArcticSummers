
###################################################
   ### COMMUNITY RDA -- Fig. 3
###################################################

amp_ordinate(ampvis, 
  type = "RDA",  
  transform = "hellinger", 
  distmeasure = "bray",
  filter_species = 0.01,
  constrain = c("depth","lon"),
  envfit_numeric = c("depth","lon","temp","iceConc"),
  envfit_numeric_arrows_scale = 0.6,
  sample_colorframe = F,
  sample_point_size = 0.2) +
  geom_point(aes(#size=lon, 
    shape=region, fill=depth, color=depth), size=2.2) +
  scale_shape_manual(values=c(17,19)) +
  scale_color_gradient2(
    low = "mediumturquoise",mid="orchid4",
    high="gray22", midpoint=50) +
  scale_fill_gradient2(
    low = "mediumturquoise",mid="orchid4",
    high="gray22", midpoint=50) +
  theme_bw() + 
  theme(axis.text = element_blank(),
        #axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right") 

## Stats
adonis2(
  dist.asv ~ lon + depth + jday, 
  data = ENV, 
  sqrt.dist = F)

adonis2(
  dist.bac ~ lon + depth + jday, 
  data = ENV, 
  sqrt.dist = F)

adonis2(
  dist.arch ~ lon + depth + jday, 
  data = ENV, 
  sqrt.dist = F)

adonis2(
  dist.env ~ lon + depth + jday, 
  data = subset(ENV, 
  sample_title %in% rownames(dist.env)), 
  sqrt.dist = F)


