#############################################

## SUPPLEMENTARY FIGURES
## Supplementary Figures for DOI


#############################################
  ## SUBSTRATES --- Fig S1
#############################################

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


#############################################
  ## DIVERSITY --- Fig S2
#############################################

# export size 6x7
AlphaDiv %>%
  drop_na(region) %>%
  reshape2::melt(
    id.vars=c("region","layer"),
    measure.vars=c("simpson","richness")) %>%
  ggplot() + 
  geom_boxplot(aes(
    x=layer, y=value, fill=layer)) +
  scale_fill_manual(values=c(
    "surface"="paleturquoise",
    "chl-max"="darkcyan",
    "bel.chl-max"="orchid4",
    "lower-photic"="gray22")) +
  facet_grid(
    variable~region, scales="free") +
  theme_classic() + 
  theme(
    axis.text.x = element_text(
      angle = 45, hjust = 1),
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.ticks = element_blank())


#############################################
  ## CORE-ASVs -- Fig S3
#############################################

# filter by location & depth
# only ASVs >0.1% in 90% of samples 
core.srf.WSC <- ASV.hel %>%
  rownames_to_column("asv") %>%
  reshape2::melt(value.name="Abundance") %>%
  dplyr::rename(sample_title=variable) %>%
  left_join(ENV) %>%
  filter(layer=="surface" & lon > 0.5) %>%
  group_by(asv) %>% 
  filter(Abundance > 0.1 & n() >= 0.9) %>%
  add_column(type="eastern Fram Strait")

core.cmax.WSC <- ASV.hel %>%
  rownames_to_column("asv") %>%
  reshape2::melt(value.name="Abundance") %>%
  dplyr::rename(sample_title=variable) %>%
  left_join(ENV) %>%
  filter(layer=="chl-max" & lon > 0.5) %>%
  group_by(asv) %>% 
  filter(Abundance > 0.1 & n() >= 0.9) %>%
  add_column(type="eastern Fram Strait")

core.deep.WSC <- ASV.hel %>%
  rownames_to_column("asv") %>%
  reshape2::melt(value.name="Abundance") %>%
  dplyr::rename(sample_title=variable) %>%
  left_join(ENV) %>%
  filter(layer=="lower-photic" & lon > 0.5) %>%
  group_by(asv) %>% 
  filter(Abundance > 0.1 & n() >= 0.9) %>%
  add_column(type="eastern Fram Strait")

############

core.srf.EGC <- ASV.hel %>%
  rownames_to_column("asv") %>%
  reshape2::melt(value.name="Abundance") %>%
  dplyr::rename(sample_title=variable) %>%
  left_join(ENV) %>%
  filter(layer=="surface" & lon < 0.5) %>%
  group_by(asv) %>% 
  filter(Abundance > 0.1 & n() >= 0.9) %>%
  add_column(type="western Fram Strait")

core.cmax.EGC <- ASV.hel %>%
  rownames_to_column("asv") %>%
  reshape2::melt(value.name="Abundance") %>%
  dplyr::rename(sample_title=variable) %>%
  left_join(ENV) %>%
  filter(layer=="chl-max" & lon < 0.5) %>%
  group_by(asv) %>% 
  filter(Abundance > 0.1 & n() >= 0.9) %>%
  add_column(type="western Fram Strait")

core.deep.EGC <- ASV.hel %>%
  rownames_to_column("asv") %>%
  reshape2::melt(value.name="Abundance") %>%
  dplyr::rename(sample_title=variable) %>%
  left_join(ENV) %>%
  filter(layer=="lower-photic" & lon < 0.5) %>%
  group_by(asv) %>% 
  filter(Abundance > 0.1 & n() >= 0.9) %>%
  add_column(type="western Fram Strait")

##########################

core.all <- rbind(
  core.srf.WSC, core.srf.EGC,
  core.cmax.EGC, core.cmax.WSC,
  core.deep.WSC, core.deep.EGC)

# select most relevant ASVs
core.all <- TAX %>%
  rownames_to_column("asv") %>%
  right_join(core.all) %>%
  filter(asv %in% c(
    "asv175","asv131","asv16","asv95","asv1",
    "asv183","asv11","asv30","asv19","asv91",
    "asv5","asv77","asv127","asv4","asv72",
    "asv75","asv70","asv12","asv126","asv23",
    "asv3","asv223","asv122","asv223",
    "asv26","asv43","asv37","asv66")) %>%
  mutate(tax=case_when(Class %in% c(
    "Dadabacteriia","SAR324 uc",
    "Thermoplasmata","Acidimicrobiia",
    "Marinimicrobia_SAR406 uc",
    "Verrucomicrobiae")~"other",
    T~Class)) %>%
  mutate(asv=paste(Genus, asv, sep="_")) %>%
  group_by(lon_short, type, asv, tax, layer) %>%
  summarize(Abundance=mean(Abundance)) %>%
  ungroup %>%
  mutate(type=factor(type, levels = c(
    "western Fram Strait","eastern Fram Strait"))) 

# Plot -- export size 8x9
ggplot(core.major) +
  geom_point(aes(
    x=type, y=asv, group=asv,
    size=Abundance, color=tax)) +
  facet_grid(
    tax~layer, scales="free", space="free", 
    labeller = label_wrap_gen(width=1)) +
  scale_color_fish_d(
    option="Acanthurus_olivaceus",
    direction = -1) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(
          angle=45, hjust=1,vjust=1),
        strip.text.y = element_text(
          angle=-90, face="bold"),
        strip.background = element_rect(
          fill="white"),
        strip.text = element_text(
          face="bold"))


#############################################
  ## ASV pres-abs -- Fig. S4
#############################################

core1 <- ASV.rel %>%
  rownames_to_column("asv") %>%
  reshape2::melt(value.name="Abundance") %>%
  dplyr::rename(sample_title=variable) %>%
  left_join(ENV) %>%
  filter(layer %in% c("surface","chl-max") & lon>0.5) %>%
  group_by(asv) %>% 
  summarize(Abundance=mean(Abundance)) %>% 
  filter(Abundance >= 0.005) %>%
  distinct(asv, .keep_all = T)  %>%
  column_to_rownames("asv")

core2 <- ASV.rel %>%
  rownames_to_column("asv") %>%
  reshape2::melt(value.name="Abundance") %>%
  dplyr::rename(sample_title=variable) %>%
  left_join(ENV) %>%
  filter(layer=="lower-photic" & lon>0.5) %>%
  group_by(asv) %>% 
  summarize(Abundance=mean(Abundance)) %>% 
  filter(Abundance >= 0.005) %>%
  distinct(asv, .keep_all = T)  %>%
  column_to_rownames("asv")

core3 <- ASV.rel %>%
  rownames_to_column("asv") %>%
  reshape2::melt(value.name="Abundance") %>%
  dplyr::rename(sample_title=variable) %>%
  left_join(ENV) %>%
  filter(layer %in% c("surface","chl-max") & lon<0.5) %>%
  group_by(asv) %>% 
  group_by(asv) %>% 
  summarize(Abundance=mean(Abundance)) %>% 
  filter(Abundance >= 0.005) %>%
  distinct(asv, .keep_all = T)  %>%
  column_to_rownames("asv")

core4 <- ASV.rel %>%
  rownames_to_column("asv") %>%
  reshape2::melt(value.name="Abundance") %>%
  dplyr::rename(sample_title=variable) %>%
  left_join(ENV) %>%
  filter(layer=="lower-photic" & lon<0.5) %>%
  group_by(asv) %>% 
  summarize(Abundance=mean(Abundance)) %>%
  filter(Abundance >= 0.005) %>%
  distinct(asv, .keep_all = T)  %>%
  column_to_rownames("asv")

# combine in list
core <- list()
core[["surface/chl-max (east)"]] <- as.character(row.names(core1))
core[["subsurface (east)"]] <- as.character(row.names(core2))
core[["surface/chl-max (west)"]] <- as.character(row.names(core3))
core[["subsurface (west)"]] <- as.character(row.names(core4))

# Function for presence-absence matrix
presAbs <- function (input) {
  # Same as original fromList()...
  elements <- unique(unlist(input))
  data <- unlist(lapply(input, function(x) {
    x <- as.vector(match(elements, x))
  }))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(input), byrow = F))
  data <- data[which(rowSums(data) != 0), ]
  names(data) <- names(input)
  # ... Except now it conserves your original value names!
  row.names(data) <- elements
  return(data)
}

overlaps <- presAbs(core) %>%
  rownames_to_column("asv")

upset(
  overlaps,
  number.angles = 0, 
  nsets=4,
  main.bar.color = "gray22",
  sets.bar.color = "gray22",
  matrix.color = "gray22",
  point.size = 2.44, line.size = 0.8, text.scale = 1.2,
  mainbar.y.label = "Shared ASVs",
  order.by = "freq")


#################################################
## Fractionated chlorophyll  -- Fig. S5
#################################################

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
