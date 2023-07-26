######################################################
  ### INTERANNUAL PATTERNS -- Fig. 6
######################################################

## MICROBES -- RDA by year; surface only 

amp_ordinate(amp_subset_samples(
  ampvis, layer %in% c("surface","chl-max")),
  type = "RDA",  
  transform = "hellinger", 
  distmeasure = "bray",
  filter_species = 0.01,
  constrain = c("year"),
  sample_shape_by = "layer",
  #sample_label_size = "lon",
  sample_color_by = "year", 
  sample_point_size = 2) +
  geom_point(aes(fill=year)) +
  scale_shape_manual(values=c(21,21)) +
  scale_fill_manual(values=c(
    "2015"="azure4",
    "2016"="cyan3",
    "2017"="midnightblue",
    "2018"="magenta3",
    "2019"="orange")) +
  scale_color_manual(values=c(
    "2015"="azure4",
    "2016"="cyan3",
    "2017"="midnightblue",
    "2018"="magenta3",
    "2019"="orange")) +
  theme_bw() + 
  theme(axis.text = element_blank(),
        #axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right") 

adonis2(
  dist~year, 
  data=subset(ENV, layer %in% c(
    "surface","chl-max")), 
  sqrt.dist=F)

###############################################

## SUBSTRATES -- PCA by year

# Select only variables with continuous data
# Omit TDN: missing 2015
PCA1 <- ENV %>%
  filter(sample_title!="2018-HG4-35" & depth < 50) %>%
  select_if(is.numeric) %>%
  dplyr::select(-c(
    "lat","lon","O2_sat","TDN","O2_conc","temp",
    "depth","iceDist","chlPast","chlConc","chl",
    "chl_0.4-3um","chl_larger3um","iceConc",
    "IceConcPast","IceDistPast",
    #,"DHAA-C","DHAA-N","year","jday",
    "sal", "Bacterial cells")) %>% 
  drop_na() 

PCA2 <- ENV %>% 
  mutate(year = str_extract(sample_title, "\\d{4}")) %>%
  mutate(region=case_when(
    lon<0~"EGC",TRUE~"WSC")) %>%
  mutate(layer=case_when(
    depth <50 ~"surface/chl-max",
    TRUE~"bel.chl-max/deep")) %>%
 semi_join(PCA1)

#color by year
autoplot(
  prcomp(PCA1, center=T, scale=T),
  data = PCA2, label=F, 
  #frame=T, frame.colour = 'season', 
  loadings=T, loadings.colour="lemonchiffon4",
  loadings.label=T, loadings.label.size=2,
  loadings.label.vjust=1.5, residuals=T,
  loadings.label.colour="black") +
  geom_point(aes(
    color=year, fill=year), size=2) +
  scale_shape_manual(values=c(25,21)) +
  scale_size(range = c(0.5,5)) + 
  scale_fill_manual(values=c(
    "2015"="azure4",
    "2016"="cyan3",
    "2017"="midnightblue",
    "2018"="magenta",
    "2019"="orange")) +
  scale_color_manual(values=c(
    "2015"="azure4",
    "2016"="cyan3",
    "2017"="midnightblue",
    "2018"="magenta",
    "2019"="orange")) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        #panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        #legend.key.size = unit(0.05, "cm"),
        legend.position = "right")

###############################################

## Corresponding bacterial differences

plot.year1 <- amp_heatmap(amp_subset_samples(
  ampvis, layer %in% c("surface","chl-max")),
  tax_aggregate = "Order",
  tax_show = c(
    "Opitutales","SAR86","SAR11",
    "Oceanospirillales","Flavobacteriales"),
  plot_values = F, 
  plot_colorscale = "sqrt",
  plot_legendbreaks = c(5, 10, 20, 40),
  min_abundance = 1, 
  max_abundance = 55,
  round = 0,
  normalise = T,
  group_by = "year",
  facet_by = "region",
  color_vector = c( 
    "white","lightsteelblue3",
    "skyblue4","midnightblue",
    "darkorange3")) +
  geom_text(aes(
    label = round(Abundance),
    color = ifelse(
      Abundance < 8,
      "black","white"))) +
  scale_color_identity()+
  theme(
    axis.ticks = element_blank(),
    axis.text.x = element_text(
      angle=45,hjust=1,vjust=1),
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.position = "right") 

###############################################

## Corresponding CHO differences

plot.year2 <- ENV %>%
  reshape2::melt(id.vars=c(
    "region","layer","year"), 
    measure.vars=c("Carbohydrates")) %>% 
  filter(!is.na(value)) %>%
  filter(layer %in% c("surface","chl-max")) %>%
  # filter(grepl("Carbohydrates", variable)) %>%
  group_by(year, layer, region) %>%
  summarize(value=mean(value)) %>%
  ungroup %>%
  ggplot() +
  geom_bar(aes(x=year, y=value), 
    position="stack", stat="identity", fill="gray22") +
  scale_y_continuous(n.breaks=4) +  
  facet_wrap(region~.,ncol=2) +
  ylab("nM CHO [0-20m]") +
  theme_classic2() + 
  theme(
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x= element_blank(),
    legend.position = "right")

# export size 4x5
cowplot::plot_grid(
  plot.year2,
  plot.year1,
  ncol=1,
  rel_heights = c(1,1.1),
  align="v",
  axis="tblr")

# for plot: calculate average julian day of annual sampling
avgJulian <- ENV %>%
  group_by(year) %>%
  summarize(jmean=mean(jday))

# table avg values; specify numeric columns as measure variables
measure_vars <- colnames(ENV)[sapply(ENV, is.numeric)]

ENV %>%
  reshape2::melt(
    id.vars=c(
    "region","year"),
    measure.vars=c(
      "temp","Amines","Carbohydrates",
      "Amino acids","DOC","Sugar acids",
      "iceConc","iceDist","Bacterial cells")) %>% 
filter(!is.na(value)) %>%
  #filter(layer %in% c("surface","chl-max")) %>%
  # filter(grepl("Carbohydrates", variable)) %>%
  group_by(year, region, variable) %>%
  summarize(value=mean(value)) %>%
  ungroup %>%
  ggplot() +
  geom_bar(aes(x=year, y=value), 
           position="stack", stat="identity", fill="gray22") +
  scale_y_continuous(n.breaks=4) +  
  facet_grid(variable~region, scales = "free") +
  ylab("nM CHO [0-20m]") +
  theme_classic2() + 
  theme(
    axis.text.x = element_text(size=10),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "right")

