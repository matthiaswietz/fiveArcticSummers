
####################################################
   ### sPLS -- Fig. 7 		
####################################################

asv <- t(ASV.hel) %>% 
  as.data.frame()

meta <- ENV %>%
  remove_rownames %>%
  column_to_rownames("sample_title") %>%
  dplyr::select(c(
    "lon","temp","depth","sal","iceDist",
    "jday","chlConc","chlPast","O2_conc",
    "chl","chl_larger3um","chl_0.4-3um",
    "O2_sat","iceConc","IceConcPast","iceDist")) 

# match rows
asv <- asv[row.names(meta),] 

# Calculate sPLS
PLS <- spls(
  asv, meta, 
  ncomp = 3, 
  tol = 1e-06,
  max.iter = 100,
  near.zero.var = F,
  #keepX = c(400,400,100), 
  #keepY= c(9,9,9), 
  mode = "regression")

# Find significant ASVs
cim <- cim(
  PLS, 
  comp = 1:3, 
  cutoff = 0.4,
  dist.method = c(
    "correlation","correlation"), 
  clust.method=c(
    "complete", "complete"), 
  mapping = "XY")

# Subset significant ASVs
mat <- cim$mat
subset = row.names(mat)

# Subset ASVs to significant ASVs
asv <- asv[names(asv) %in% subset]

################################

# Final PLS with subset
PLS <- spls(
  asv, meta, 
  ncomp = 3, 
  tol = 1e-06,
  max.iter = 100,
  near.zero.var = F,
  mode = "regression")

cim <- cim(
  PLS,
  comp = 1:3, 
  #cutoff = 0.3,
  dist.method = c(
    "correlation","correlation"), 
  clust.method=c(
    "complete", "complete"), 
  mapping = "XY")

################################

# Extract correlations & clusters
corr = data.frame(cim$mat) %>%
  mutate(across(everything(), round, 2))

hc <- hclust(Dist(
  cim$mat, method="correlation"),"complete")	
cutree <- cutree(hc, k=11)

# Plot clusters 
#plot(as.phylo(hc), tip.color=col.cluster[cutree], cex=0.3)

# Reformat
clusters <- cbind(corr, cutree)
clusters <- data.frame(
  cbind(rownames(clusters), clusters[,"cutree"])) %>%
  dplyr::rename(asv=X1, cluster=X2) %>%
  mutate(cluster=paste0("C", cluster))

# Define colors
col.cluster <- c(
  "C1" = "plum4",  
  "C2" = "aquamarine2",  # mediumorchid3 darkcyan
  "C3" = "palevioletred1", #
  "C4" = "black", 
  "C5" = "gray55", 
  "C6" = "yellow2",
  "C7" = "dodgerblue4",
  "C8" = "greenyellow",
  "C9" = "orangered4",
  "C10" = "orange1")

## Combine cluster 10-11
tax <- TAX[match(
  clusters$asv, rownames(TAX)),] %>%
  add_column(cluster=clusters$cluster) %>%
  mutate(cluster=case_when(
  cluster %in% c("C10","C11")~"C10", 
    TRUE~cluster)) %>%
  rownames_to_column("asv")

# reassign to clustCol
clustCol <- tibble(
  cluster=names(col.cluster)) %>%
  left_join(tax) %>%  
  dplyr::slice(mixedorder(asv)) %>% 
  drop_na() 

# Add new colors for combined clusters
clustCol$col <- col.cluster[match(
  clustCol$cluster, names(col.cluster))]

# Plot
cim(
  PLS,  
  dist.method = c(
    "correlation","correlation"), 
  clust.method=c(
    "complete","complete"), 
  margins = c(7,20), 
  #col.cex = 0.2,
  row.names = F,
  mapping = "XY",
  #center = T,
  row.sideColors = clustCol$col,
  color = colorRampPalette(c(
    "khaki4","seashell3",
    "floralwhite","lightskyblue2",
    "skyblue4"))(20),
  save = "pdf",
  name.save = "../plots/Fig7_PLS") 

################################

## SUMMARIES

# ASVs per cluster
clustNum <- clustCol %>%
  group_by(cluster) %>%
  tally()

# Export cluster info for SI table
clustCol %>%
  dplyr::slice(mixedorder(cluster)) %>% 
  dplyr::select(-c("Kingdom","Species","col")) %>%
  write.table(file="../ASV_PLS-clust.txt",
  sep="\t", row.names=F, col.names=T, quote=F)

# Export correlations for SI table
corr %>% rownames_to_column("asv") %>%
  dplyr::slice(mixedorder(asv)) %>% 
  write.table(file="../ASV_PLS-corr.txt",
  sep="\t", row.names=F, col.names=T, quote=F)

# more unclassifieds in EGC/deep?
clustComp <- clustCol %>% 
  group_by(cluster) %>%
  summarize(count=sum(grepl("uc", Genus))) %>%
  left_join(clustNum) %>%
  mutate(frac=count/n*100)

# Export for SI table
clustComp %>%
  dplyr::slice(mixedorder(cluster)) %>% 
  relocate(cluster, n, count, frac) %>%
  dplyr::rename(
    `#ASV`=n,
    `#unclassified`=count,
    `%unclassified`=frac) %>%
  write.table(
    file="../ASV_PLS-comp.txt",
    sep="\t", row.names=F, col.names=T, quote=F)


################################################
   ###  SIGNATURE GENERA  ###
################################################

clustAsv = ASV.hel %>%
  rownames_to_column("asv")

# per ASV, select maximum value -- HELLINGER
# filter low-abundant ASVs
clustTax = ASV.hel %>%
  rownames_to_column("asv") %>% 
  reshape2::melt(id.vars=c("asv")) %>% 
  group_by(asv) %>% 
  summarise(value = max(value)) %>%
  left_join(clustCol) %>%
  drop_na() %>%
  #filter(value > 0.44) %>%  #for rel abd
  filter(value > 0.05) %>%
  left_join(clustAsv) %>%
  dplyr::select(-c(
    Kingdom, Phylum, Order, Family, Species)) %>%
  reshape2::melt() %>% 
  aggregate(
    value~cluster+Class+Genus+variable, 
    data=., FUN=sum) %>% 
  group_by(Genus) %>% 
  top_n(1, value) %>% 
  distinct(value, .keep_all=T) %>%
  mutate(cluster=factor(cluster, levels=c(
    "C1","C2","C3","C4","C5","C6",
    "C7","C8","C9","C10")))

# Sort genera by cluster
clustTax$Genus <- factor(
  clustTax$Genus, levels=rev(unique(
    clustTax$Genus[order(clustTax$cluster)])))

# filter most uc taxa / low-abundant
# filter ambiguties (e.g. Magnetospiraceae-uc/Magnetospira) 
# export size 9x4
ggplot(data=subset(clustTax, !Genus %in% c(
  "Flavobacteriales uc","Flavobacteriaceae uc",
  "Proteobacteria uc","Cellvibrionaceae uc","PeM15 uc",
  "SAR11 uc","Cryomorphaceae uc","Rhodobacteraceae uc",
  "Gammaproteobacteria uc","Cyclobacteriaceae uc","OM75",
  "JL-ETNP-F27","Marine_Benthic_Group_A uc","MB11C04",
  "Parvibaculaceae uc","KI89A uc",#"NS9 uc","NS7 uc",
  "Fluviicola","NS5",
  "Rickettsiales uc","SCGC_AAA164-E04","Porticoccus",
  "Marinoscillum","Gimesiaceae uc","Shewanella","OM27",
  "Magnetospira","Colwellia","Cyanobacteriia uc",
  "Winogradskyella","Vicingus","CL500-3"))) + 
  geom_point(aes(
    x=cluster, y=Genus,
    size=value, colour=cluster), 
    stat="identity", position="identity") + 
  scale_colour_manual(
    values=col.cluster) + 
  scale_size_continuous(
    range = c(1, 5), 
    breaks=c(0.3,0.6,0.9), 
    name="Maximum abundance") + 
  guides(fill=T, color="none") +
  scale_y_discrete(limits=rev) + 
  theme_bw() + 
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_text(
      size = 10, colour = "black"),
    legend.position="top")

# percent contribution of signature-pops -- layer/region
clustFrac = ASV.rel %>%
  rownames_to_column("asv") %>% 
  reshape2::melt(id.vars=c("asv")) %>% 
  left_join(clustCol) %>%
  mutate(type=case_when(
    is.na(cluster)~"unassigned", TRUE~"signature")) %>% 
  group_by(variable, type) %>% 
  summarize(sum=sum(value)) %>% 
  left_join(ENV, by=c("variable"="sample_title")) %>%
  group_by(type, region, layer) %>% 
  summarize(mean=mean(sum)) %>%
  ungroup()

# percent contribution of signature-pops -- year/region
clustYear = ASV.rel %>%
  rownames_to_column("asv") %>% 
  reshape2::melt(id.vars=c("asv")) %>% 
  left_join(clustCol) %>%
  group_by(variable,cluster) %>% 
  summarize(sum=sum(value)) %>% 
  #mutate(year = gsub("-.*","", variable)) %>% 
 left_join(ENV, by=c("variable"="sample_title")) %>%
  group_by(cluster, year, region) %>% 
  summarize(mean=mean(sum)) %>%
  ungroup()

# export for SI table
clustFrac %>%
  filter(type=="signature") %>%
  arrange(region) %>%
  dplyr::rename(
    `Relative abundance of signature clusters (%)`=mean) %>%
  dplyr::select(-c("type")) %>%
  mutate_if(is.numeric, round, 2) %>%
  write.table(
    file="../ASV_PLS-frac.txt",
    sep="\t", row.names=F, col.names=T, quote=F)
