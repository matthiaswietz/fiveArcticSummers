###############################################
  ### CORRELATIONS -- Fig. 5
###############################################

# Subset by abundance 
filter <- apply(ASV.rel, 1, function(row) {
  sum(row >= 0.5) >= 3})

# Filter accordingly
subset <- ASV.rel[filter, ] 

# Set filtered taxinfo
tax <- TAX[filter, ] %>%
  rownames_to_column("asv")

# Hellinger transform
subset = as.data.frame(
  apply(subset, 2, function(x) sqrt(x / sum(x)))) %>%
  rownames_to_column("asv")

# Join taxinfo + ENV
subset <- TAX %>%
  rownames_to_column("asv") %>%
  right_join(subset) 

# Subset: surface WSC
subset1 <- subset %>%
  reshape2::melt(value.name="Abundance") %>%
  dplyr::rename(sample_title=variable) %>%
  left_join(ENV) %>%
  filter(lon > 0.5 & layer %in% c("surface","chl-max")) %>% #
  pivot_wider(
    id_cols = c("sample_title"), 
    names_from = asv, 
    values_from = c("Abundance")) %>%
  column_to_rownames("sample_title")

# Subset: surface EGC
subset2 <- subset %>%
  reshape2::melt(value.name="Abundance") %>%
  dplyr::rename(sample_title=variable) %>%
  left_join(ENV) %>%
  filter(lon < 0.5 & layer %in% c("surface","chl-max")) %>% #
  pivot_wider(
    id_cols = c("sample_title"), 
    names_from = asv, 
    values_from = c("Abundance")) %>%
  column_to_rownames("sample_title")

# Select variables
meta <- ENV %>%
  remove_rownames() %>%
  column_to_rownames("sample_title") %>%
  select_if(is.numeric) %>%
  dplyr::select(-c(
    "sal","IceConcPast","IceDistPast",
    "iceConc","chl","chl_0.4-3um",
    "chl_larger3um","chlPast","chlConc",
    "Bacterial cells","DHAA-C","DHAA-N",
    "O2_sat","iceDist",#"temp","depth",
    "O2_conc","TDN","lat","lon"))

#########################################

## Correlations: WSC

# Match rows
meta1 <- meta[row.names(subset1),]

# Compute
cor <- corr.test(
  subset1, meta1, 
  use = "pairwise",
  method="spearman",
  adjust = "BH",
  alpha=.05, ci=T,
  minlength=5, normal=T)

# Extract, subset most significant; 
r <- cor$r 
r[abs(r) < -0.4 | abs(r) < 0.4 ] = NA 

# Remove taxa w/ depth+jday correlations
# Reformat; remove all-NA cols+rows
r <- as.data.frame(r) %>%
  filter(if_all(c(depth, jday), ~is.na(.x))) %>%
  dplyr::select(-c("depth","jday")) %>%
  filter(if_any(everything(), ~!is.na(.))) %>%
  rownames_to_column("asv") %>%
  distinct(asv, .keep_all=T) %>%
  reshape2::melt() %>%
  mutate_if(is.numeric, round, 2) 

# Reformat p-values
p <- as.data.frame(cor$p) %>%
  mutate(across(
    everything(), ~replace(.,.>0.05, NA))) %>%
  rownames_to_column("asv") %>%
  distinct(asv, .keep_all=T) %>%
  reshape2::melt() %>%
  mutate_if(is.numeric, ~case_when(
    . < 0.05 & . > 0.01 ~ "*",
    . < 0.01 & . > 0.001 ~ "**",
    . < 0.001 ~ "***")) %>%
  dplyr::rename(pvalue=value)

# Merge, reformat
cor1 <- left_join(r, p) %>%
  unite(sign, value, pvalue, sep="", remove=F) %>%
  mutate_at(c("sign"), ~gsub("NA", NA, .)) %>%
 # mutate_at(c("Class"), ~gsub(
  #  "proteobacteria|microbiae|eroidia", "", .)) %>%
  mutate(across(value, as.numeric)) %>%
  #drop_na(OTU) %>%
  add_column(region="eastern Fram Strait")

#########################################

# Correlations: EGC 

# Match rows
meta2 <- meta[row.names(subset2),]

# Compute
cor <- corr.test(
  subset2, meta2, 
  use = "pairwise",
  method="spearman",
  adjust = "BH",
  alpha=.05, ci=T,
  minlength=5, normal=T)

# Extract, subset most significant; 
r <- cor$r 
r[abs(r) < -0.4 | abs(r) < 0.4 ] = NA 

# Remove taxa w/ depth+jday correlations
# Reformat; remove all-NA cols+rows
r <- as.data.frame(r) %>%
  filter(if_all(c(depth, jday), ~is.na(.x))) %>%
  dplyr::select(-c("depth","jday")) %>%
  filter(if_any(everything(), ~!is.na(.))) %>%
  rownames_to_column("asv") %>%
  distinct(asv, .keep_all=T) %>%
  reshape2::melt() %>%
  mutate_if(is.numeric, round, 2) 

# Reformat p-values
p <- as.data.frame(cor$p) %>%
  mutate(across(
    everything(), ~replace(.,.>0.05, NA))) %>%
  rownames_to_column("asv") %>%
  distinct(asv, .keep_all=T) %>%
  reshape2::melt() %>%
  mutate_if(is.numeric, ~case_when(
    . < 0.05 & . > 0.01 ~ "*",
    . < 0.01 & . > 0.001 ~ "**",
    . < 0.001 ~ "***")) %>%
  dplyr::rename(pvalue=value)

# Merge, reformat
cor2 <- left_join(r, p) %>%
  unite(sign, value, pvalue, sep="", remove=F) %>%
  mutate_at(c("sign"), ~gsub("NA", NA, .)) %>%
  mutate(across(value, as.numeric)) %>%
  #drop_na(OTU) %>%
  add_column(region="western Fram Strait")

#########################################

## MERGE EVERYTHING

cor.micro <- rbind(cor1, cor2) %>%
  complete(asv, region, variable) %>%
  left_join(tax) %>%
 # mutate_at(c("Class"), ~gsub(
  #  "proteobacteria|microbiae|eroidia", "", .)) %>%
  mutate(across(value, as.numeric)) %>%
  mutate(pvalue = ifelse(is.na(value), NA, pvalue)) %>%
  #drop_na() %>%
  mutate(id = paste(asv,Genus,sep="-")) %>%
  mutate(region=factor(region, levels=c(
    "western Fram Strait","eastern Fram Strait")))

# Export as SI table
cor.micro %>% drop_na %>%
  group_by(id, variable, value) %>%
  summarise(n = dplyr::n(), .groups="drop") %>%
  filter(n > 1L) %>%
  pivot_wider(
    names_from = variable,
    id_cols = id,
    values_from = value) %>%
  write.table(file="../ASV_SpearmanCorr.txt")

#########################################

## Total correlations -- Fig. 5a

# Export size 3x4
cor.micro %>% 
  drop_na() %>%
    group_by(region, variable) %>%
  summarize(count=n()) %>%
  ungroup() %>%
  mutate(variable=factor(variable, levels=c(
    "DOC","Amines","Sugar acids","Amino acids",
    "Carbohydrates","temp"))) %>%
  ggplot() +
  geom_bar(aes(
    x=variable, y=count, fill=region), stat="identity",
    position=position_dodge(preserve="single")) +
  scale_fill_manual(values=c(
    "western Fram Strait"="deepskyblue2",
    "eastern Fram Strait"="khaki2")) +
  ylab("Number of correlations") +
  scale_y_continuous(
    expand=c(0.01,0.01)) +
  theme_classic() + 
  theme(
    axis.title.x = element_blank(),
    axis.text.x=element_text(
      angle=45, vjust=1,hjust=1),
    axis.ticks = element_blank())

# Correlation by taxon -- Fig. 5b
cor.micro %>% 
  drop_na() %>%
  mutate(Order=case_when(!Order %in% c(
    "Thiotrichales","Rhodobacterales","SAR86",
    "SAR11","Oceanospirillales","Cellvibrionales",
    "PeM15","Thiomicrospirales","Alteromonadales",
    "Verrucomicrobiales","Cytophagales",
    "Flavobacteriales")~"Other",
    TRUE~Order)) %>%
  group_by(region, variable, Order) %>%
  summarize(count=n()) %>%
  filter(variable %in% c("Carbohydrates","temp","DOC") & 
      count>1) %>%
  mutate(Order=factor(Order, levels=c(
    "Other","Cellvibrionales","Thiotrichales",
    "Rhodobacterales","SAR86","SAR11","Oceanospirillales",
    "Verrucomicrobiales","Thiomicrospirales",
    "Alteromonadales","PeM15","Cytophagales",
    "Flavobacteriales"))) %>%
  mutate(variable=factor(variable, levels=c(
    "DOC","Carbohydrates","temp"))) %>%
  ggplot(aes(
    x=Order, y=variable, 
    fill=count, label=count)) + 
  geom_tile() +
  geom_text(color="white", size=4) +
  scale_fill_gradientn(
    colors=colorRampPalette(c(
      "gray84","gray55","gray28",
      "gray16","gray2"))(10),
     na.value = "white",
    limits = c(1, 10),
    breaks = c(1, 5,10)) +
  coord_flip() +
  facet_grid(
    ~region,  
    #scales="free", 
    space="free") +
  coord_flip() +
  theme_classic() +
  theme(axis.ticks = element_blank(),
        axis.text.x = element_text(
          angle=45, hjust=1, vjust=1),
        # panel.grid = element_blank(),
        #panel.border = element_blank(),
        # strip.background = element_blank(),
        #legend.position = "none",
        # plot.background = theme_rect(colour = "blue"),
        axis.title = element_blank())

