
###########################################################
  ### LATITUDINAL PATTERNS -- ENV.PARAMETERS & SUBSTRATES
###########################################################

## Environmental parameters -- Fig. 2a
# Round cell numbers (10e6)
# Export size 4.5 x 4.5
ENV %>% 
 mutate(Cells=`Bacterial cells`/1000000)  %>% 
  reshape2::melt(id.vars=c(
    "lon_short","layer"), 
    measure.vars=c(
      "temp","sal","chl",
      "IceConcPast","Cells")) %>% 
  filter(!is.na(value) & layer!="bel.chl-max") %>%
  filter(!(layer=="lower-photic" & variable=="IceConcPast")) %>%
  mutate(variable=factor(
    variable, levels = c(
      "IceConcPast","temp","sal","chl","Cells"))) %>%
  ggplot(aes(
    x=lon_short, y=value, 
    group=layer, color=layer)) +
  stat_summary(
    fun = mean,
    fun.min = function(x) mean(x) - sd(x), 
    fun.max = function(x) mean(x) + sd(x), 
    geom = "pointrange") +
  stat_summary(
    fun = mean,
    geom = "line", 
    linewidth=0.8) +
  scale_color_manual(values=c(
    "surface"="mediumturquoise",
    "chl-max"="mediumorchid4",
    "lower-photic"="gray22")) +
  facet_grid(
    variable~., scale="free") +
  theme_classic2() +
  xlab("Longitude") +
  theme(
    axis.text.x = element_text(angle = 0),
    axis.title.y = element_blank())

#################################################

## Full lineplot -- Fig S2
# Export size 6.5 x 4
ENV %>% reshape2::melt(id.vars=c(
  "sample_title","layer","lon_short"), 
  measure.vars=c(
    "Carbohydrates","Amino acids","Sugar acids",
    "Amines","TDN","DOC")) %>%
  filter(!is.na(value) & layer!="bel.chl-max") %>%
  drop_na(value) %>%
  group_by(sample_title, lon_short, layer, variable) %>%
  summarize(value=sum(value)) %>%
  ungroup %>%
  mutate(layer=factor(layer, level=c(
    "surface","chl-max","lower-photic")))  %>%
  #group_by(lon_short, layer, variable) %>%
  # summarize(value=mean(value)) %>%
  ungroup() %>%
  drop_na %>%
  ggplot(aes(
    x=lon_short, y=value, 
    group=layer, color=layer)) +
  stat_summary(
    fun = mean,
    fun.min = function(x) mean(x) - sd(x), 
    fun.max = function(x) mean(x) + sd(x), 
    geom = "pointrange") +
  stat_summary(
    fun = mean,
    geom = "line", 
    linewidth=0.8) +
  xlab("Longitude") +
  scale_color_manual(values=c(
    "surface"="mediumturquoise",
    "chl-max"="mediumorchid4",
    "lower-photic"="gray22")) +
  facet_grid(
    variable~., scale="free") +
  theme_classic2() + 
  theme(
    axis.text.x = element_text(angle = 0),
    axis.title.y = element_blank())

#################################################

## Substrates -- Fig. 2b
# export size 5x5
ENV %>%
  reshape2::melt(id.vars=c(
    "sample_title","layer","lon_short"), measure.vars=c(
    "Carbohydrates","Amino acids","Sugar acids",
    "Amines","TDN","DOC")) %>%
  filter(!is.na(value) & layer!="bel.chl-max") %>%
  drop_na(value) %>%
  group_by(sample_title, lon_short, layer, variable) %>%
  summarize(value=sum(value)) %>%
  ungroup %>%
  mutate(layer=factor(layer, level=c(
    "surface","chl-max","lower-photic")))  %>%
  mutate(variable=factor(variable, level=c(
    "Carbohydrates","Amino acids","Sugar acids",
    "Amines","DOC","TDN")))  %>%
  group_by(lon_short, layer, variable) %>%
  summarize(value=mean(value)) %>%
  ungroup() %>%
  drop_na %>%
  ggplot() +
  geom_bar(aes(
    x=lon_short, y=value, fill=variable), 
    position="stack", stat="identity") +
  scale_fill_fish_d(option="Lepomis_megalotis") +
  scale_y_continuous(n.breaks = 4) +  
  facet_grid(variable~layer, scale="free") +
  xlab("Longitude") +
  theme_classic() + 
  theme(
    axis.text.x = element_text(
      angle = 90, vjust=0.5),
    axis.ticks = element_blank(),
    axis.title.y = element_blank())

#################################################

## Fractionated chlorophyll
# export size 3x5
ENV %>%
  filter(depth<50) %>%
  reshape2::melt(
    id.vars=c("region","year"),
    measure.vars=c("chl_larger3um","chl_0.4-3um","chl")) %>% 
  filter(!is.na(value)) %>%
  group_by(year, region, variable) %>%
  mutate(variable=factor(variable, levels=c(
    "chl","chl_larger3um","chl_0.4-3um"))) %>%
  summarize(value=mean(value)) %>%
  ungroup %>%
  ggplot() +
  geom_bar(aes(
    x=variable, y=value, fill=variable), 
    position="stack", stat="identity") +
  scale_fill_manual(values=c(
    "chl_larger3um"="aquamarine3",
    "chl_0.4-3um"="lightcyan3",
    "chl"="seagreen4")) +
  scale_y_continuous(n.breaks = 4) +  
  facet_wrap(region~.,ncol=2) +
  ylab("chl-a [µM]") +
  theme_classic2() + 
  theme(
    axis.text.x = element_text(angle=45,vjust=1,hjust=1),
    axis.ticks = element_blank(),
    axis.title.x= element_blank(),
    legend.position = "none")

#################################################

## Stats

vars <- select_if(ENV, is.numeric)
id <- setdiff(names(ENV), names(vars))

stats <- ENV %>%
  reshape2::melt(
    value.vars=vars,
    id.vars=id) %>% 
  filter(layer!="bel.chl-max" & !variable %in% c(
    "simpson","richness")) %>%
  mutate_at(c("variable"), as.character) %>%
  drop_na(value) %>%
  group_by(sample_title, variable, region, layer) %>%
  summarize(value=mean(value)) %>%
  ungroup 

## REGION
# significant: CHO-WSC-chlmax, AA-WSC-chlmax, 
# Amines_WSC-chlmax, SUgarAcids-WSC-chlmax 
# DOC EGC-surface
with(subset(
  stats, layer==c("chl-max") & variable==c("Carbohydrates")),
  kruskal.test(value~region))

## LAYER
# significant: CHO,AA,Amines,SugarAcids,DOC,TDN 
with(subset(
  stats, variable==c("DOC")),
  kruskal.test(value~layer))

