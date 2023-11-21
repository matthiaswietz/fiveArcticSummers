
###################################################
  ## BROAD COMMUNITY SHIFTS -- Fig. 4a
###################################################

# Reformat taxdata
tax <- TAX %>%
  rownames_to_column("asv")

# Plot classes
ASV.rel %>% 
  rownames_to_column("asv") %>%
  reshape2::melt(
    variable.name="sample_title",
    value.name="Abundance") %>%
  left_join(ENV) %>%
  left_join(tax) %>%
  group_by(sample_title, Class, lon_short, layer) %>%
  summarise(Abundance=sum(Abundance)) %>%
  ungroup %>%
  group_by(Class, lon_short, layer) %>%
  summarize(Abundance=mean(Abundance)) %>% 
  ungroup %>%
  mutate(layer=factor(
    layer, levels = c(
      "surface","chl-max",
      "bel.chl-max","lower-photic"))) %>%
  filter(layer!="bel.chl-max") %>%
  mutate(Class=case_when(!Class %in% c(
    "Gammaproteobacteria","Bacteroidia",
    "Alphaproteobacteria","Verrucomicrobiae",
    "Planctomycetes")~"other",
    TRUE~Class)) %>%
  mutate(Class=factor(Class, levels=c(
    "Alphaproteobacteria",
    "Gammaproteobacteria","Bacteroidia",
    "Verrucomicrobiae","other"))) %>%
  ggplot() +
  geom_bar(aes(
    x=lon_short, y=Abundance, fill=Class), 
    stat="identity",position="stack") +
  facet_grid(
    layer~., scales="free") +
  scale_fill_manual(values=c(
    "Alphaproteobacteria"="lightseagreen",
    "Gammaproteobacteria"="cadetblue2",
    "Bacteroidia"="seashell2",
    "Verrucomicrobiae"="firebrick3",
    #"Planctomycetes"="mediumpurple3",
    "other"="lightgoldenrod2"))+
  theme_classic2() +
  theme(axis.ticks = element_blank(),
        axis.title = element_blank(),
        strip.text.y = element_text(
          angle = 0),
        strip.background = element_blank(),
        strip.text = element_text(hjust=0))


###################################################
  ## ORDERS BY LAYER + LAT -- Fig. 4b
###################################################

abundances <- ASV.hel %>%
  rownames_to_column("asv") %>%
  reshape2::melt(value.name="Abundance") %>%
  dplyr::rename(sample_title=variable) %>%
  left_join(ENV) 

abundances <- TAX %>%
  rownames_to_column("asv") %>%
  right_join(abundances)

# Find top100 lat/depth correlated (PLS)
max.lon <- corr %>% data.frame() %>%
  rownames_to_column("asv") %>%
  dplyr::select(c("asv","lon")) %>%
  slice_max(
    order_by = abs(lon), 
    n = 100) %>%
  dplyr::rename(cor=lon) %>%
  left_join(abundances)

max.depth <- corr %>% data.frame() %>%
  rownames_to_column("asv") %>%
  dplyr::select(c("asv","depth")) %>%
  slice_max(
    order_by = abs(depth), 
    n = 100) %>%
  dplyr::rename(cor=depth) %>%
  left_join(abundances)

# Combine 
max.deplon <- bind_rows(
  max.depth, max.lon) %>%
  mutate(layer=factor(layer, levels = c(
    "lower-photic","bel.chl-max","chl-max","surface"))) %>%
  group_by(lon_short, layer, Class, Order) %>%
  summarize(Abundance=mean(Abundance))

# Plot; export size 8x4
max.deplon %>% filter(layer!="bel.chl-max" & Order %in% c(
  "Cytophagales","Pirellulales",
  "Dadabacteriales","SAR324 uc",
  "Nitrosopumilales","SAR11",
  "Nitrospinales","SAR86","OM182",
  "Thiotrichales","Puniceispirillales",
  #"Marinimicrobia_SAR406 uc","Thiomicrospirales",
  "Flavobacteriales")) %>%
  mutate(Order=factor(Order, levels = c(
    "Cytophagales","Pirellulales",
    "Dadabacteriales","SAR324 uc",
    #"Marinimicrobia_SAR406 uc",
    "Nitrospinales","Nitrosopumilales",
    "SAR11",#"Rhodobacterales", 
    "Flavobacteriales","SAR86",
    "Thiomicrospirales", "Thiotrichales","OM182",
    "Puniceispirillales"))) %>%
  ggplot() +
  geom_point(aes(
    x=lon_short, y=layer, group=Order,
    size=Abundance, color=Abundance)) +
  scale_size(
    trans="log",
    breaks = c(0.0001,0.001,0.005,0.06),
    range=c(0.2,4)) +
  facet_grid(
    Order~., scales="free") +
  scale_color_gradientn(colors = c(
    "aliceblue","lightyellow1","lightyellow3",
    "darkseagreen3","purple4","gray8"), 
    trans="log",
    breaks = c(0.0001,0.001,0.005,0.06)) +
  xlab("Longitude") +
  theme_classic2() +
  theme(axis.ticks = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(
          angle = 45, vjust=0.5),
        strip.text.y = element_text(
          angle = 0),
        strip.background = element_blank(),
        strip.text = element_text(hjust=0))
