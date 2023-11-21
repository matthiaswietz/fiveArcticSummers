
###############################################
   ###  MAG-ASV comparisons -- Fig. 8
################################################

setwd("./MAGs")

# Load MAG annotations 
files <- list.files(
  pattern="_gene", recursive=F)

# Reformat
anno <- lapply(
  files, fread, fill=T, quote=F, data.table=F)
names(anno) <- gsub(
  "_gene_table_complete.txt$", "", files)

# Add MAG name as new column
anno <- Map(cbind, anno, MAG=names(anno))

# Delete unnecessary cols 
anno <- lapply(anno, function(x) x[!(names(x) %in% c(
  "Contig","Start","Stop","Length"))])

# Load KEGG identifiers; remove leading space
KEGG.id <- read.delim(
  "KEGG_id.txt", sep="\t", header=T) %>% 
  group_by(Knumber) %>% 
  distinct(Knumber, .keep_all=T) %>%
  ungroup %>%
  mutate_all(str_trim)

# To dataframe; remove all-NA columns
anno <- anno %>%
  bind_rows() %>%
  #na_if("") %>% 
  select_if(~!all(is.na(.))) %>%
  #mutate(seq = gsub("FRAM18_", "", seq)) %>%
  left_join(KEGG.id, by=c(
    "KEGG_identifier"="Knumber")) %>%
  group_by(Pathway, dbCAN) %>%
  filter(n() >= 3) %>%
  mutate(KO = str_extract(Pathway, "ko[0-9]*")) %>%
  mutate(Pathway = gsub(" \\[.*","", Pathway))

################################

# Count genes per MAG
MAG.genes <- anno %>%
  group_by(MAG) %>%
  tally()

# Load ASV-MAG pairs
# select sPLS populations
MAG.hits <- read.table(
  "MAG_Hits.txt", header=T, 
  check.names=F, sep="\t") %>%
  left_join(clustCol) %>%
  left_join(MAG.genes) %>%
  distinct(MAG, .keep_all=T) %>%
  drop_na(cluster) %>%
  mutate(type=case_when(
    cluster %in% c("C1","C2","C6","C8","C9","C10")~"surface",
    cluster %in% c("C3","C4","C5","C7")~"lower-photic",
    TRUE ~ as.character(NA))) 

#######################################

# Combine CAZy subfamilies (e.g. PL6_1 / PL6_3 to PL6)
# calculate fraction by surface/deep
cazy <- anno %>%
  right_join(MAG.hits, by=c("MAG")) %>%
  drop_na(dbCAN) %>%
  group_by(MAG, dbCAN, n, cluster, Order) %>%
  summarise(count=n()) %>% 
  #drop_na() %>%
  mutate(frac=count/n*100) %>%
  ungroup() %>% 
  #replace(is.na(.), "unassigned") %>%
  ungroup() %>% as_tibble() %>%
  #distinct(MAG, dbCAN, .keep_all=T) %>%
  mutate(CAZY = case_when(
    grepl("GH", dbCAN)~"GH", 
    grepl("CE", dbCAN)~"CE",
    grepl("AA", dbCAN)~"AA",
    grepl("GT", dbCAN)~"GT",
    grepl("PL", dbCAN)~"PL",
    grepl("CB", dbCAN)~"CBM", 
    TRUE~"unassigned")) %>%
  filter(dbCAN!="unassigned") %>%
  mutate(dbCAN = gsub("_.*","", dbCAN)) %>%
  mutate(type=case_when(
    cluster %in% c("C1","C2","C6","C8","C9","C10")~"surface",
    cluster %in% c("C3","C4","C5","C7")~"lower-photic",
    TRUE ~ as.character(NA))) 

#######################################

## Pres/abs - Fig. 8a

plot.mag1 <- cazy %>%
#  mutate(dbCAN = gsub("[A-z]","", dbCAN)) %>%
  filter(CAZY %in% c("GH","PL") & type!="unassigned") %>%
  mutate(
    type=factor(type, levels=c("surface","lower-photic")),
    CAZY=factor(CAZY, levels=c("GH","PL"))) %>%
  distinct(type, dbCAN, .keep_all=T) %>%
  mutate(cazyOrdered = factor(
    dbCAN, levels=rev(mixedsort(unique(dbCAN))))) %>%
  arrange(cazyOrdered) %>%
  mutate(frac=case_when(
    frac > 0~"present", TRUE~"absent")) %>%
  ggplot() +
  geom_tile(aes(
    x=cazyOrdered, y=frac, fill=type)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0.02, 0.5)) + 
  facet_grid(
    CAZY~type, scales="free", space="free") +
  xlab("CAZY family") +
  coord_flip() +
  scale_fill_manual(values=c(
    "surface"="gray64",
    "lower-photic"="gray28")) +
  theme_classic2() +
  theme(
    axis.text.x = element_blank(),
    strip.background.y = element_blank(),
     strip.text.y = element_blank(),
    #panel.grid = element_blank(),
    legend.position = "none",
    axis.ticks = element_blank(),
    axis.title.x = element_blank())

#######################################

## Fraction by taxon -- Fig. 8b

plot.mag2 <- cazy %>%
  mutate(dbCAN = gsub("[A-z]","", dbCAN)) %>%
  filter(CAZY %in% c("GH","PL")) %>%
  mutate(
    type=factor(type, levels=c("surface","lower-photic")),
    CAZY=factor(CAZY, levels=c("GH","PL"))) %>%
  reshape2::melt(
    id.vars=c("type","CAZY","Order","MAG"), 
    measure.vars=c("frac")) %>% 
  group_by(type, CAZY, Order, MAG) %>%
  summarize(value=sum(value)) %>%
  #filter(!is.na(frac) & type!="unassigned") %>%
  group_by(type, CAZY, Order) %>%
  summarize(value=mean(value)) %>%
  #arrange(value) %>%
  ungroup %>%
  ggplot() +
  geom_bar(aes(
    x=reorder(Order, -value, FUN=sum), y=value, fill=CAZY), 
    position="stack", stat="identity") +
  scale_x_discrete(expand = c(0, 1)) +
  scale_y_continuous(
    n.breaks=4, expand = c(0, 0.01)) + 
  scale_fill_manual(values=c(
    "GH"="paleturquoise4","CBM"="sienna2",
    "PL"="goldenrod1")) +
  facet_grid(
    .~type, scales="free_x") +
  ylab("Fraction") +
  theme_classic2() + 
  theme(
    axis.text.x = element_text(
      angle=45, hjust=1, vjust=1.05),
    #strip.background.x = element_blank(),
    #strip.text.x = element_blank(),
    #panel.grid = element_blank(),
    legend.position = "bottom",
    axis.ticks = element_blank(),
    axis.title.x = element_blank())

cowplot::plot_grid(
  plot.mag1,
  plot.mag2,
  ncol=2,
  rel_heights = c(1,0.4),
  rel_widths = c(0.4,1),
  align="v",
  axis="lr")

